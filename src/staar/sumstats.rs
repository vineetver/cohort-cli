use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, BooleanBuilder, Float64Builder, Int32Builder, ListBuilder, StringBuilder,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use faer::Mat;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use serde::{Deserialize, Serialize};

use super::genotype::GenotypeResult;
use super::null_model::NullModel;
use crate::db::DuckEngine;
use crate::db::query_strings;
use crate::error::FavorError;
use crate::output::Output;
use crate::types::{AnnotatedVariant, AnnotationWeights};

const SEGMENT_BP: u32 = 500_000;
const MAX_SEGMENT_VARIANTS: usize = 2000;

#[derive(Serialize, Deserialize)]
pub struct StudyMeta {
    pub favor_meta_version: u32,
    pub trait_type: String,
    pub trait_name: String,
    pub n_samples: usize,
    pub sigma2: f64,
    pub maf_cutoff: f64,
    pub covariates: Vec<String>,
    pub segment_size: u32,
}

/// Export MetaSTAAR summary statistics for all chromosomes.
///
/// Matches R MetaSTAAR scaling: both score vector and covariance are divided
/// by the null model dispersion sigma squared.
pub fn emit(
    engine: &DuckEngine,
    geno: &GenotypeResult,
    variants: &[AnnotatedVariant],
    null: &NullModel,
    n_samples: usize,
    output_dir: &Path,
    meta: &StudyMeta,
    out: &dyn Output,
) -> Result<(), FavorError> {
    out.status("MetaSTAAR: computing summary statistics...");

    let chromosomes = query_strings(engine, "SELECT DISTINCT chrom FROM _rare ORDER BY chrom")?;

    for chrom in &chromosomes {
        let indices: Vec<usize> = variants
            .iter()
            .enumerate()
            .filter(|(_, v)| v.chromosome.to_string() == format!("chr{chrom}") || v.chromosome.to_string() == *chrom)
            .map(|(i, _)| i)
            .collect();
        if indices.is_empty() {
            continue;
        }

        let dir = output_dir.join(format!("chromosome={chrom}"));
        std::fs::create_dir_all(&dir)?;
        emit_chromosome(engine, geno, variants, &indices, null, n_samples, chrom, &dir, out)?;
    }

    let meta_json = serde_json::to_string_pretty(meta)
        .map_err(|e| FavorError::Resource(format!("Failed to serialize meta_staar.json: {e}")))?;
    std::fs::write(output_dir.join("meta_staar.json"), meta_json)?;

    out.success(&format!("Summary statistics -> {}", output_dir.display()));
    Ok(())
}

fn emit_chromosome(
    engine: &DuckEngine,
    geno: &GenotypeResult,
    variants: &[AnnotatedVariant],
    chrom_indices: &[usize],
    null: &NullModel,
    n: usize,
    chrom: &str,
    dir: &Path,
    out: &dyn Output,
) -> Result<(), FavorError> {
    let mut coarse: std::collections::BTreeMap<i32, Vec<usize>> =
        std::collections::BTreeMap::new();
    for &gi in chrom_indices {
        coarse
            .entry((variants[gi].position / SEGMENT_BP) as i32)
            .or_default()
            .push(gi);
    }
    let mut segments: Vec<(i32, Vec<usize>)> = Vec::new();
    let mut next_id = 0i32;
    for (_, indices) in coarse {
        for chunk in indices.chunks(MAX_SEGMENT_VARIANTS) {
            segments.push((next_id, chunk.to_vec()));
            next_id += 1;
        }
    }
    out.status(&format!(
        "  chr{chrom}: {} variants, {} segments",
        chrom_indices.len(),
        segments.len()
    ));

    let geno_path = format!(
        "{}/chromosome={chrom}/data.parquet",
        geno.output_dir.display()
    );
    let inv_s2 = 1.0 / null.sigma2;

    let chrom_positions: Vec<u32> = chrom_indices.iter().map(|&i| variants[i].position).collect();
    let extract_cols = super::genotype::dosage_columns(n);
    let geno_flat =
        super::genotype::load(engine, &geno_path, &chrom_positions, n, &extract_cols)?;
    let pos_to_flat: std::collections::HashMap<u32, usize> = chrom_positions
        .iter()
        .enumerate()
        .map(|(local, &p)| (p, local))
        .collect();

    let total_variants = chrom_indices.len();
    let mut b_position = Int32Builder::with_capacity(total_variants);
    let mut b_ref = StringBuilder::with_capacity(total_variants, total_variants * 4);
    let mut b_alt = StringBuilder::with_capacity(total_variants, total_variants * 4);
    let mut b_maf = Float64Builder::with_capacity(total_variants);
    let mut b_mac = Int32Builder::with_capacity(total_variants);
    let mut b_n_obs = Int32Builder::with_capacity(total_variants);
    let mut b_u_stat = Float64Builder::with_capacity(total_variants);
    let mut b_v_stat = Float64Builder::with_capacity(total_variants);
    let mut b_segment_id = Int32Builder::with_capacity(total_variants);
    let mut b_gene = StringBuilder::with_capacity(total_variants, total_variants * 8);
    let mut b_region = StringBuilder::with_capacity(total_variants, total_variants * 8);
    let mut b_consequence = StringBuilder::with_capacity(total_variants, total_variants * 8);
    let mut b_cadd = Float64Builder::with_capacity(total_variants);
    let mut b_revel = Float64Builder::with_capacity(total_variants);
    let mut b_cage_prom = BooleanBuilder::with_capacity(total_variants);
    let mut b_cage_enh = BooleanBuilder::with_capacity(total_variants);
    let mut b_ccre_prom = BooleanBuilder::with_capacity(total_variants);
    let mut b_ccre_enh = BooleanBuilder::with_capacity(total_variants);
    let mut b_weights: [Float64Builder; 11] = std::array::from_fn(|_| Float64Builder::with_capacity(total_variants));

    let n_segments = segments.len();
    let mut s_segment_id = Int32Builder::with_capacity(n_segments);
    let mut s_n_variants = Int32Builder::with_capacity(n_segments);
    let mut s_positions = ListBuilder::new(Int32Builder::with_capacity(total_variants));
    let mut s_refs = ListBuilder::new(StringBuilder::with_capacity(total_variants, total_variants * 4));
    let mut s_alts = ListBuilder::new(StringBuilder::with_capacity(total_variants, total_variants * 4));
    let mut s_cov_lower = ListBuilder::new(Float64Builder::with_capacity(
        total_variants * (total_variants + 1) / 2,
    ));

    for (seg_id, seg_indices) in &segments {
        let seg_id = *seg_id;
        let m = seg_indices.len();

        let mut g = Mat::zeros(n, m);
        for (col, &gi) in seg_indices.iter().enumerate() {
            if let Some(&flat_idx) = pos_to_flat.get(&variants[gi].position) {
                let base = flat_idx * n;
                for row in 0..n {
                    g[(row, col)] = geno_flat[base + row];
                }
            }
        }

        let u = g.transpose() * &null.residuals;
        let pg = null.project(&g);
        let k = g.transpose() * &pg;

        for (j, &gi) in seg_indices.iter().enumerate() {
            let v = &variants[gi];
            b_position.append_value(v.position as i32);
            b_ref.append_value(&v.ref_allele);
            b_alt.append_value(&v.alt_allele);
            b_maf.append_value(v.maf);
            b_mac.append_value((2.0 * v.maf * n as f64).round() as i32);
            b_n_obs.append_value(n as i32);
            b_u_stat.append_value(u[(j, 0)] * inv_s2);
            b_v_stat.append_value(k[(j, j)] * inv_s2);
            b_segment_id.append_value(seg_id);
            b_gene.append_value(&v.gene_name);
            b_region.append_value(&v.annotation.region_type);
            b_consequence.append_value(&v.annotation.consequence);
            b_cadd.append_value(v.annotation.cadd_phred);
            b_revel.append_value(v.annotation.revel);
            b_cage_prom.append_value(v.annotation.regulatory.cage_promoter);
            b_cage_enh.append_value(v.annotation.regulatory.cage_enhancer);
            b_ccre_prom.append_value(v.annotation.regulatory.ccre_promoter);
            b_ccre_enh.append_value(v.annotation.regulatory.ccre_enhancer);
            for (i, builder) in b_weights.iter_mut().enumerate() {
                builder.append_value(v.annotation.weights.0[i]);
            }
        }

        s_segment_id.append_value(seg_id);
        s_n_variants.append_value(m as i32);

        let pos_builder = s_positions.values();
        for &gi in seg_indices {
            pos_builder.append_value(variants[gi].position as i32);
        }
        s_positions.append(true);

        let ref_builder = s_refs.values();
        for &gi in seg_indices {
            ref_builder.append_value(&variants[gi].ref_allele);
        }
        s_refs.append(true);

        let alt_builder = s_alts.values();
        for &gi in seg_indices {
            alt_builder.append_value(&variants[gi].alt_allele);
        }
        s_alts.append(true);

        let cov_builder = s_cov_lower.values();
        for i in 0..m {
            for j in 0..=i {
                cov_builder.append_value(k[(i, j)] * inv_s2);
            }
        }
        s_cov_lower.append(true);
    }

    write_variants_parquet(
        dir,
        &mut b_position, &mut b_ref, &mut b_alt, &mut b_maf, &mut b_mac,
        &mut b_n_obs, &mut b_u_stat, &mut b_v_stat, &mut b_segment_id,
        &mut b_gene, &mut b_region, &mut b_consequence,
        &mut b_cadd, &mut b_revel,
        &mut b_cage_prom, &mut b_cage_enh, &mut b_ccre_prom, &mut b_ccre_enh,
        &mut b_weights,
    )?;

    write_segments_parquet(
        dir,
        &mut s_segment_id,
        &mut s_n_variants,
        &mut s_positions,
        &mut s_refs,
        &mut s_alts,
        &mut s_cov_lower,
    )?;

    out.status(&format!("    chr{chrom} done"));
    Ok(())
}

fn variant_schema() -> Schema {
    let mut fields = vec![
        Field::new("position", DataType::Int32, false),
        Field::new("ref_allele", DataType::Utf8, false),
        Field::new("alt_allele", DataType::Utf8, false),
        Field::new("maf", DataType::Float64, false),
        Field::new("mac", DataType::Int32, false),
        Field::new("n_obs", DataType::Int32, false),
        Field::new("u_stat", DataType::Float64, false),
        Field::new("v_stat", DataType::Float64, false),
        Field::new("segment_id", DataType::Int32, false),
        Field::new("gene_name", DataType::Utf8, false),
        Field::new("region_type", DataType::Utf8, false),
        Field::new("consequence", DataType::Utf8, false),
        Field::new("cadd_phred", DataType::Float64, false),
        Field::new("revel", DataType::Float64, false),
        Field::new("is_cage_promoter", DataType::Boolean, false),
        Field::new("is_cage_enhancer", DataType::Boolean, false),
        Field::new("is_ccre_promoter", DataType::Boolean, false),
        Field::new("is_ccre_enhancer", DataType::Boolean, false),
    ];
    for name in AnnotationWeights::NAMES {
        fields.push(Field::new(name, DataType::Float64, false));
    }
    Schema::new(fields)
}

fn segment_schema() -> Schema {
    Schema::new(vec![
        Field::new("segment_id", DataType::Int32, false),
        Field::new("n_variants", DataType::Int32, false),
        Field::new(
            "positions",
            DataType::List(Arc::new(Field::new("item", DataType::Int32, true))),
            false,
        ),
        Field::new(
            "refs",
            DataType::List(Arc::new(Field::new("item", DataType::Utf8, true))),
            false,
        ),
        Field::new(
            "alts",
            DataType::List(Arc::new(Field::new("item", DataType::Utf8, true))),
            false,
        ),
        Field::new(
            "cov_lower",
            DataType::List(Arc::new(Field::new("item", DataType::Float64, true))),
            false,
        ),
    ])
}

fn parquet_props() -> WriterProperties {
    WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .set_max_row_group_size(5000)
        .build()
}

fn write_variants_parquet(
    dir: &Path,
    b_position: &mut Int32Builder,
    b_ref: &mut StringBuilder,
    b_alt: &mut StringBuilder,
    b_maf: &mut Float64Builder,
    b_mac: &mut Int32Builder,
    b_n_obs: &mut Int32Builder,
    b_u_stat: &mut Float64Builder,
    b_v_stat: &mut Float64Builder,
    b_segment_id: &mut Int32Builder,
    b_gene: &mut StringBuilder,
    b_region: &mut StringBuilder,
    b_consequence: &mut StringBuilder,
    b_cadd: &mut Float64Builder,
    b_revel: &mut Float64Builder,
    b_cage_prom: &mut BooleanBuilder,
    b_cage_enh: &mut BooleanBuilder,
    b_ccre_prom: &mut BooleanBuilder,
    b_ccre_enh: &mut BooleanBuilder,
    b_weights: &mut [Float64Builder; 11],
) -> Result<(), FavorError> {
    let schema = Arc::new(variant_schema());
    let mut columns: Vec<ArrayRef> = vec![
        Arc::new(b_position.finish()),
        Arc::new(b_ref.finish()),
        Arc::new(b_alt.finish()),
        Arc::new(b_maf.finish()),
        Arc::new(b_mac.finish()),
        Arc::new(b_n_obs.finish()),
        Arc::new(b_u_stat.finish()),
        Arc::new(b_v_stat.finish()),
        Arc::new(b_segment_id.finish()),
        Arc::new(b_gene.finish()),
        Arc::new(b_region.finish()),
        Arc::new(b_consequence.finish()),
        Arc::new(b_cadd.finish()),
        Arc::new(b_revel.finish()),
        Arc::new(b_cage_prom.finish()),
        Arc::new(b_cage_enh.finish()),
        Arc::new(b_ccre_prom.finish()),
        Arc::new(b_ccre_enh.finish()),
    ];
    for b in b_weights.iter_mut() {
        columns.push(Arc::new(b.finish()));
    }

    let batch = RecordBatch::try_new(schema.clone(), columns)
        .map_err(|e| FavorError::Resource(format!("Arrow batch: {e}")))?;

    let file = File::create(dir.join("variants.parquet"))
        .map_err(|e| FavorError::Resource(format!("Create variants.parquet: {e}")))?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(parquet_props()))
        .map_err(|e| FavorError::Resource(format!("Parquet writer init: {e}")))?;
    writer
        .write(&batch)
        .map_err(|e| FavorError::Resource(format!("Parquet write: {e}")))?;
    writer
        .close()
        .map_err(|e| FavorError::Resource(format!("Parquet close: {e}")))?;

    Ok(())
}

fn write_segments_parquet(
    dir: &Path,
    s_segment_id: &mut Int32Builder,
    s_n_variants: &mut Int32Builder,
    s_positions: &mut ListBuilder<Int32Builder>,
    s_refs: &mut ListBuilder<StringBuilder>,
    s_alts: &mut ListBuilder<StringBuilder>,
    s_cov_lower: &mut ListBuilder<Float64Builder>,
) -> Result<(), FavorError> {
    let schema = Arc::new(segment_schema());
    let columns: Vec<ArrayRef> = vec![
        Arc::new(s_segment_id.finish()),
        Arc::new(s_n_variants.finish()),
        Arc::new(s_positions.finish()),
        Arc::new(s_refs.finish()),
        Arc::new(s_alts.finish()),
        Arc::new(s_cov_lower.finish()),
    ];

    let batch = RecordBatch::try_new(schema.clone(), columns)
        .map_err(|e| FavorError::Resource(format!("Arrow batch: {e}")))?;

    let file = File::create(dir.join("segments.parquet"))
        .map_err(|e| FavorError::Resource(format!("Create segments.parquet: {e}")))?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(parquet_props()))
        .map_err(|e| FavorError::Resource(format!("Parquet writer init: {e}")))?;
    writer
        .write(&batch)
        .map_err(|e| FavorError::Resource(format!("Parquet write: {e}")))?;
    writer
        .close()
        .map_err(|e| FavorError::Resource(format!("Parquet close: {e}")))?;

    Ok(())
}
