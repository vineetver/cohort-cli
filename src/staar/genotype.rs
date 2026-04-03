//! Multi-sample VCF → genotype parquet (dense packed).
//!
//! Dosages are packed into a single FLOAT[] list column per variant — not one
//! column per sample. This gives DuckDB a single compression envelope to decode
//! instead of N_SAMPLES separate column chunks. For 3000 samples, that's 3000x
//! less column metadata overhead and sequential I/O instead of random seeks.
//!
//! Schema: chromosome, position, ref, alt, maf, dosages FLOAT[N_SAMPLES]
//! Sample order matches the VCF header (stored in a sidecar).

use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use arrow::array::{ArrayRef, Float32Builder, Int32Builder, ListBuilder, StringBuilder};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use serde::{Deserialize, Serialize};

use crate::error::FavorError;
use crate::ingest::vcf::{normalize_chrom, parsimony_normalize};
use crate::output::Output;

#[derive(Serialize, Deserialize)]
pub struct GenotypeMeta {
    pub version: u32,
    pub n_samples: usize,
    pub chromosomes: Vec<String>,
    pub source_vcf: String,
}

pub struct GenotypeResult {
    pub sample_names: Vec<String>,
    pub output_dir: PathBuf,
}

pub fn extract_genotypes(
    vcf_path: &Path,
    output_dir: &Path,
    available_memory: u64,
    output: &dyn Output,
) -> Result<GenotypeResult, FavorError> {
    let is_bgzf = vcf_path.extension().map(|e| e == "gz" || e == "bgz").unwrap_or(false);

    // Single header read
    let sample_names = read_sample_names(vcf_path, is_bgzf)?;
    let n_samples = sample_names.len();
    if n_samples == 0 {
        return Err(FavorError::Input(
            "VCF has no samples. STAAR requires a multi-sample VCF.".into(),
        ));
    }

    output.status(&format!("Extracting genotypes: {} samples", n_samples));

    // Batch size: each variant = 5 metadata fields + n_samples * 4 bytes dosage
    let bytes_per_variant = (n_samples as u64) * 4 + 200;
    let batch_size = ((available_memory / 4) / bytes_per_variant).max(1000).min(100_000) as usize;

    // Schema: 5 metadata columns + 1 packed dosage list
    let schema = Arc::new(packed_schema());
    let geno_dir = output_dir.join("genotypes");
    std::fs::create_dir_all(&geno_dir)?;

    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .set_max_row_group_size(batch_size)
        .build();

    let mut state = ExtractState {
        current_chrom: None,
        writer: None,
        batch: PackedBatchBuilder::new(n_samples, batch_size),
        total_variants: 0,
        chromosomes: Vec::new(),
    };

    // Parse VCF — macro deduplicates bgzf vs plain reader dispatch.
    macro_rules! extract_records {
        ($reader:expr) => {{
            let mut vcf_reader = noodles_vcf::io::Reader::new($reader);
            let _header = vcf_reader.read_header()?;
            for result in vcf_reader.records() {
                let record = result.map_err(|e| FavorError::Analysis(format!("VCF parse: {e}")))?;
                process_record(&record, n_samples, &mut state, &schema, &props, &geno_dir, output)?;
            }
        }};
    }

    let file = File::open(vcf_path)?;
    if is_bgzf {
        let bgzf = noodles_bgzf::Reader::new(file);
        extract_records!(BufReader::with_capacity(256 * 1024, bgzf));
    } else {
        extract_records!(BufReader::with_capacity(256 * 1024, file));
    }

    flush(&mut state, &schema)?;
    if let Some(w) = state.writer.take() {
        w.close().map_err(|e| FavorError::Resource(format!("Parquet close: {e}")))?;
    }

    // Write sample names sidecar (order matters for unpacking dosages)
    let sidecar = geno_dir.join("samples.txt");
    std::fs::write(&sidecar, sample_names.join("\n"))?;

    // Write metadata — this is the last step so partial extractions are detectable
    let meta = GenotypeMeta {
        version: 1,
        n_samples,
        chromosomes: state.chromosomes.clone(),
        source_vcf: vcf_path.display().to_string(),
    };
    std::fs::write(
        geno_dir.join("genotypes.json"),
        serde_json::to_string_pretty(&meta)
            .map_err(|e| FavorError::Resource(format!("JSON serialize: {e}")))?,
    )?;

    output.success(&format!(
        "Extracted {} variants × {} samples → packed dosage lists",
        state.total_variants, n_samples,
    ));

    Ok(GenotypeResult {
        sample_names,
        output_dir: geno_dir,
    })
}

struct ExtractState {
    current_chrom: Option<String>,
    writer: Option<ArrowWriter<File>>,
    batch: PackedBatchBuilder,
    total_variants: u64,
    chromosomes: Vec<String>,
}

fn flush(state: &mut ExtractState, schema: &Arc<Schema>) -> Result<(), FavorError> {
    if state.batch.count > 0 {
        if let Some(w) = state.writer.as_mut() {
            let rb = state.batch.finish(schema)?;
            w.write(&rb).map_err(|e| FavorError::Resource(format!("Parquet write: {e}")))?;
        }
    }
    Ok(())
}

fn switch_chrom(
    chrom: &str, state: &mut ExtractState, schema: &Arc<Schema>,
    props: &WriterProperties, geno_dir: &Path, output: &dyn Output,
) -> Result<(), FavorError> {
    flush(state, schema)?;
    if let Some(w) = state.writer.take() {
        w.close().map_err(|e| FavorError::Resource(format!("Parquet close: {e}")))?;
    }
    let chr_dir = geno_dir.join(format!("chromosome={chrom}"));
    std::fs::create_dir_all(&chr_dir)?;
    let f = File::create(chr_dir.join("data.parquet"))
        .map_err(|e| FavorError::Resource(format!("Create: {e}")))?;
    state.writer = Some(
        ArrowWriter::try_new(f, schema.clone(), Some(props.clone()))
            .map_err(|e| FavorError::Resource(format!("Writer init: {e}")))?,
    );
    state.current_chrom = Some(chrom.to_string());
    state.chromosomes.push(chrom.to_string());
    output.status(&format!("  chr{chrom}..."));
    Ok(())
}

fn process_record(
    record: &noodles_vcf::Record, n_samples: usize,
    state: &mut ExtractState, schema: &Arc<Schema>,
    props: &WriterProperties, geno_dir: &Path, output: &dyn Output,
) -> Result<(), FavorError> {
    let raw_chrom = record.reference_sequence_name();
    let chrom = match normalize_chrom(raw_chrom) {
        Some(c) => c,
        None => return Ok(()),
    };

    if state.current_chrom.as_deref() != Some(chrom) {
        switch_chrom(chrom, state, schema, props, geno_dir, output)?;
    }

    let pos = match record.variant_start() {
        Some(Ok(p)) => p.get() as i32,
        _ => return Ok(()),
    };

    let ref_allele = record.reference_bases().to_uppercase();
    let alt_bases = record.alternate_bases();
    let alt_str: &str = alt_bases.as_ref();
    let alts: Vec<&str> = alt_str.split(',').collect();

    let samples_raw = record.samples();
    let samples_str: &str = samples_raw.as_ref();
    let gt_index = 0; // GT is first FORMAT field by VCF spec

    let sample_fields: Vec<&str> = if samples_str.is_empty() || samples_str == "." {
        Vec::new()
    } else {
        samples_str.split('\t').collect()
    };

    for (alt_idx, alt) in alts.iter().enumerate() {
        let alt_upper = alt.trim().to_uppercase();
        if alt_upper == "*" || alt_upper == "." || alt_upper.is_empty() { continue; }

        let (norm_ref, norm_alt, norm_pos) = parsimony_normalize(&ref_allele, &alt_upper, pos);

        // Parse genotypes — stack buffer, no heap allocation for typical diploid GTs
        let mut dosages = vec![f32::NAN; n_samples];
        let mut ac: f64 = 0.0;
        let mut an: f64 = 0.0;

        for (i, sf) in sample_fields.iter().enumerate().take(n_samples) {
            let gt = extract_gt_field(sf, gt_index);
            let dose = gt_to_dosage(gt.as_bytes(), (alt_idx + 1) as u8);
            dosages[i] = dose;
            if dose.is_finite() {
                ac += dose as f64;
                an += 2.0;
            }
        }

        let af = if an > 0.0 { ac / an } else { 0.0 };
        let maf = af.min(1.0 - af) as f32;

        state.batch.push(chrom, norm_pos, &norm_ref, &norm_alt, maf, &dosages);
        state.total_variants += 1;

        if state.batch.is_full() {
            flush(state, schema)?;
            if state.total_variants % 500_000 == 0 {
                output.status(&format!("  {} variants...", state.total_variants));
            }
        }
    }

    Ok(())
}

fn packed_schema() -> Schema {
    Schema::new(vec![
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("position", DataType::Int32, false),
        Field::new("ref", DataType::Utf8, false),
        Field::new("alt", DataType::Utf8, false),
        Field::new("maf", DataType::Float32, true),
        Field::new("dosages",
            DataType::List(Arc::new(Field::new("item", DataType::Float32, true))),
            false),
    ])
}

struct PackedBatchBuilder {
    chromosome: StringBuilder,
    position: Int32Builder,
    ref_allele: StringBuilder,
    alt_allele: StringBuilder,
    maf: Float32Builder,
    dosages: ListBuilder<Float32Builder>,
    count: usize,
    capacity: usize,
}

impl PackedBatchBuilder {
    fn new(n_samples: usize, capacity: usize) -> Self {
        Self {
            chromosome: StringBuilder::with_capacity(capacity, capacity * 3),
            position: Int32Builder::with_capacity(capacity),
            ref_allele: StringBuilder::with_capacity(capacity, capacity * 4),
            alt_allele: StringBuilder::with_capacity(capacity, capacity * 4),
            maf: Float32Builder::with_capacity(capacity),
            dosages: ListBuilder::with_capacity(
                Float32Builder::with_capacity(capacity * n_samples),
                capacity,
            ),
            count: 0,
            capacity,
        }
    }

    fn push(&mut self, chrom: &str, pos: i32, ref_a: &str, alt_a: &str, maf: f32, doses: &[f32]) {
        self.chromosome.append_value(chrom);
        self.position.append_value(pos);
        self.ref_allele.append_value(ref_a);
        self.alt_allele.append_value(alt_a);
        self.maf.append_value(maf);

        let values = self.dosages.values();
        for &d in doses {
            if d.is_nan() { values.append_null(); } else { values.append_value(d); }
        }
        self.dosages.append(true);

        self.count += 1;
    }

    fn is_full(&self) -> bool { self.count >= self.capacity }

    fn finish(&mut self, schema: &Arc<Schema>) -> Result<RecordBatch, FavorError> {
        let columns: Vec<ArrayRef> = vec![
            Arc::new(self.chromosome.finish()),
            Arc::new(self.position.finish()),
            Arc::new(self.ref_allele.finish()),
            Arc::new(self.alt_allele.finish()),
            Arc::new(self.maf.finish()),
            Arc::new(self.dosages.finish()),
        ];
        self.count = 0;
        RecordBatch::try_new(schema.clone(), columns)
            .map_err(|e| FavorError::Resource(format!("Arrow batch: {e}")))
    }
}

/// Extract GT field from a colon-delimited sample field. Zero allocation.
#[inline]
fn extract_gt_field<'a>(sample_field: &'a str, gt_index: usize) -> &'a str {
    let mut start = 0;
    let mut field_idx = 0;
    let bytes = sample_field.as_bytes();
    for i in 0..bytes.len() {
        if bytes[i] == b':' {
            if field_idx == gt_index {
                return &sample_field[start..i];
            }
            field_idx += 1;
            start = i + 1;
        }
    }
    if field_idx == gt_index { &sample_field[start..] } else { "." }
}

/// Parse GT bytes to dosage. Branchless for common cases, no allocation.
/// "0/0" → 0.0, "0/1" → 1.0, "1/1" → 2.0, "./." → NaN
#[inline]
fn gt_to_dosage(gt: &[u8], alt_index: u8) -> f32 {
    // Fast path: most common genotypes are 3 bytes (0/0, 0/1, 1/1, ./.)
    if gt.len() == 3 {
        let a0 = gt[0];
        let a1 = gt[2];
        if a0 == b'.' || a1 == b'.' { return f32::NAN; }
        let d0 = if a0.wrapping_sub(b'0') == alt_index { 1.0f32 } else { 0.0 };
        let d1 = if a1.wrapping_sub(b'0') == alt_index { 1.0f32 } else { 0.0 };
        return d0 + d1;
    }
    // Slow path: multi-digit alleles or missing
    gt_to_dosage_slow(gt, alt_index)
}

fn gt_to_dosage_slow(gt: &[u8], alt_index: u8) -> f32 {
    let mut dose = 0.0f32;
    for allele in std::str::from_utf8(gt).unwrap_or(".").split(|c: char| c == '/' || c == '|') {
        if allele == "." { return f32::NAN; }
        if let Ok(idx) = allele.parse::<u8>() {
            if idx == alt_index { dose += 1.0; }
        }
    }
    dose
}

pub fn read_sample_names(vcf_path: &Path, is_bgzf: bool) -> Result<Vec<String>, FavorError> {
    macro_rules! read_header_samples {
        ($reader:expr) => {{
            let mut vcf_reader = noodles_vcf::io::Reader::new($reader);
            let header = vcf_reader.read_header().map_err(|e| FavorError::Input(format!("VCF header: {e}")))?;
            Ok(header.sample_names().iter().map(|s| s.to_string()).collect())
        }};
    }

    let file = File::open(vcf_path)?;
    if is_bgzf {
        let bgzf = noodles_bgzf::Reader::new(file);
        read_header_samples!(BufReader::with_capacity(256 * 1024, bgzf))
    } else {
        read_header_samples!(BufReader::with_capacity(256 * 1024, file))
    }
}

// ---------------------------------------------------------------------------
// Genotype loading from parquet (formerly geno_load.rs)
// ---------------------------------------------------------------------------

use faer::Mat;
use crate::db::DuckEngine;

pub fn dosage_columns(n_samples: usize) -> String {
    (1..=n_samples)
        .map(|i| format!("COALESCE(list_extract(g.dosages, {i}), 0)::DOUBLE"))
        .collect::<Vec<_>>()
        .join(", ")
}

/// Load ALL genotypes from a chromosome parquet file (no position filter).
/// Use when the entire chromosome fits in the genotype budget.
/// Returns flat column-major buffer: `buf[variant * n_samples + sample]`.
pub fn load_all(
    engine: &DuckEngine,
    geno_path: &str,
    n_samples: usize,
    extract_cols: &str,
) -> Result<(Vec<f64>, Vec<u32>), FavorError> {
    let sql = format!(
        "SELECT g.position, {extract_cols} FROM read_parquet('{geno_path}') g ORDER BY g.position"
    );
    let conn = engine.connection();
    let mut stmt = conn.prepare(&sql)
        .map_err(|e| FavorError::Analysis(format!("Genotype load: {e}")))?;
    let mut rows = stmt.query([])
        .map_err(|e| FavorError::Analysis(format!("Genotype load: {e}")))?;
    let mut flat = Vec::new();
    let mut positions = Vec::new();
    while let Ok(Some(row)) = rows.next() {
        positions.push(row.get::<_, i32>(0).unwrap_or(0) as u32);
        let base = flat.len();
        flat.resize(base + n_samples, 0.0);
        for si in 0..n_samples {
            flat[base + si] = row.get::<_, f64>(si + 1).unwrap_or(0.0);
        }
    }
    Ok((flat, positions))
}

/// Load genotypes for specific positions from a chromosome parquet file.
/// Returns flat column-major buffer: `buf[variant * n_samples + sample]`.
/// Positions are returned in sorted order (parquet scan order).
pub fn load(
    engine: &DuckEngine,
    geno_path: &str,
    positions: &[u32],
    n_samples: usize,
    extract_cols: &str,
) -> Result<Vec<f64>, FavorError> {
    let pos_str = positions.iter().map(|p| p.to_string()).collect::<Vec<_>>().join(",");
    let sql = format!(
        "SELECT {extract_cols} FROM read_parquet('{geno_path}') g \
         WHERE g.position IN ({pos_str}) ORDER BY g.position"
    );
    let n_variants = positions.len();
    let mut flat = vec![0.0; n_variants * n_samples];
    let conn = engine.connection();
    let mut stmt = conn.prepare(&sql)
        .map_err(|e| FavorError::Analysis(format!("Genotype load: {e}")))?;
    let mut rows = stmt.query([])
        .map_err(|e| FavorError::Analysis(format!("Genotype load: {e}")))?;
    let mut vi = 0;
    while let Ok(Some(row)) = rows.next() {
        if vi >= n_variants { break; }
        let base = vi * n_samples;
        for si in 0..n_samples {
            flat[base + si] = row.get::<_, f64>(si).unwrap_or(0.0);
        }
        vi += 1;
    }
    Ok(flat)
}

pub fn to_mat(flat: &[f64], n_samples: usize, n_variants: usize) -> Mat<f64> {
    Mat::from_fn(n_samples, n_variants, |row, col| flat[col * n_samples + row])
}
