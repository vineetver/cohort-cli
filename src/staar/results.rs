use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{ArrayRef, Float64Builder, Int32Builder, StringBuilder, UInt32Builder};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use serde_json::json;

use crate::db::DuckEngine;
use crate::db::query_scalar;
use crate::error::FavorError;
use crate::output::Output;
use crate::staar::{self, GeneResult, MaskType};
use crate::types::{AnnotatedVariant, AnnotationWeights};

pub fn write_individual_results(
    _engine: &DuckEngine,
    pvals: &[(usize, f64)],
    variants: &[AnnotatedVariant],
    output_dir: &Path,
    out: &dyn Output,
) -> Result<(), FavorError> {
    let out_path = output_dir.join("individual.parquet");
    let n = pvals.len();

    // Sort by p-value for output ordering
    let mut sorted: Vec<(usize, f64)> = pvals.to_vec();
    sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut b_chrom = StringBuilder::with_capacity(n, n * 2);
    let mut b_pos = Int32Builder::with_capacity(n);
    let mut b_ref = StringBuilder::with_capacity(n, n * 4);
    let mut b_alt = StringBuilder::with_capacity(n, n * 4);
    let mut b_maf = Float64Builder::with_capacity(n);
    let mut b_gene = StringBuilder::with_capacity(n, n * 8);
    let mut b_region = StringBuilder::with_capacity(n, n * 8);
    let mut b_consequence = StringBuilder::with_capacity(n, n * 8);
    let mut b_cadd = Float64Builder::with_capacity(n);
    let mut b_pvalue = Float64Builder::with_capacity(n);

    for &(idx, pval) in &sorted {
        let v = &variants[idx];
        b_chrom.append_value(v.chromosome.to_string());
        b_pos.append_value(v.position as i32);
        b_ref.append_value(&v.ref_allele);
        b_alt.append_value(&v.alt_allele);
        b_maf.append_value(v.maf);
        b_gene.append_value(&v.gene_name);
        b_region.append_value(&v.annotation.region_type);
        b_consequence.append_value(&v.annotation.consequence);
        b_cadd.append_value(v.annotation.cadd_phred);
        b_pvalue.append_value(pval);
    }

    let schema = Arc::new(Schema::new(vec![
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("position", DataType::Int32, false),
        Field::new("ref_allele", DataType::Utf8, false),
        Field::new("alt_allele", DataType::Utf8, false),
        Field::new("maf", DataType::Float64, false),
        Field::new("gene_name", DataType::Utf8, false),
        Field::new("region_type", DataType::Utf8, false),
        Field::new("consequence", DataType::Utf8, false),
        Field::new("cadd_phred", DataType::Float64, false),
        Field::new("pvalue", DataType::Float64, false),
    ]));
    let columns: Vec<ArrayRef> = vec![
        Arc::new(b_chrom.finish()), Arc::new(b_pos.finish()),
        Arc::new(b_ref.finish()), Arc::new(b_alt.finish()),
        Arc::new(b_maf.finish()), Arc::new(b_gene.finish()),
        Arc::new(b_region.finish()), Arc::new(b_consequence.finish()),
        Arc::new(b_cadd.finish()), Arc::new(b_pvalue.finish()),
    ];
    let batch = RecordBatch::try_new(schema.clone(), columns)
        .map_err(|e| FavorError::Resource(format!("Arrow batch: {e}")))?;

    let file = File::create(&out_path)
        .map_err(|e| FavorError::Resource(format!("Create individual.parquet: {e}")))?;
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .build();
    let mut writer = ArrowWriter::try_new(file, schema, Some(props))
        .map_err(|e| FavorError::Resource(format!("Parquet writer: {e}")))?;
    writer.write(&batch).map_err(|e| FavorError::Resource(format!("Parquet write: {e}")))?;
    writer.close().map_err(|e| FavorError::Resource(format!("Parquet close: {e}")))?;

    let n_sig = pvals.iter().filter(|(_, p)| *p < 5e-8).count();
    out.success(&format!("  individual -> {} variants, {} genome-wide significant",
        pvals.len(), n_sig));
    Ok(())
}

pub fn write_results(
    engine: &DuckEngine,
    all_mask_results: &[(MaskType, Vec<GeneResult>)],
    trait_names: &[String],
    maf_cutoff: f64,
    output_dir: &Path,
    null_model: &staar::null_model::NullModel,
    trait_type: staar::TraitType,
    n: usize,
    out: &dyn Output,
) -> Result<(), FavorError> {
    out.status("Step 6/6: Writing results...");
    let mut significant_genes: Vec<serde_json::Value> = Vec::new();

    let n_rare: i64 = query_scalar(engine, "SELECT COUNT(*) FROM _rare")?;

    let channels = AnnotationWeights::DISPLAY_NAMES;
    let n_channels = channels.len();

    for (mask_type, results) in all_mask_results {
        if results.is_empty() { continue; }
        let out_path = output_dir.join(format!("{}.parquet", mask_type.file_stem()));

        // Sort by STAAR-O p-value
        let mut sorted_results: Vec<&GeneResult> = results.iter().collect();
        sorted_results.sort_by(|a, b| a.staar.staar_o.partial_cmp(&b.staar.staar_o).unwrap_or(std::cmp::Ordering::Equal));

        let nr = sorted_results.len();
        let mut b_ensembl = StringBuilder::with_capacity(nr, nr * 16);
        let mut b_symbol = StringBuilder::with_capacity(nr, nr * 12);
        let mut b_chrom = StringBuilder::with_capacity(nr, nr * 2);
        let mut b_start = UInt32Builder::with_capacity(nr);
        let mut b_end = UInt32Builder::with_capacity(nr);
        let mut b_nvariants = UInt32Builder::with_capacity(nr);
        let mut b_cmac = UInt32Builder::with_capacity(nr);

        // 6 base + 6*n_channels per-annotation + 6 per-test omnibus + acat_o + staar_o = 6 + 6*11 + 6 + 2 = 80
        let n_pval_cols = 6 + 6 * n_channels + 6 + 2;
        let mut pval_builders: Vec<Float64Builder> = (0..n_pval_cols)
            .map(|_| Float64Builder::with_capacity(nr))
            .collect();

        for r in &sorted_results {
            let s = &r.staar;
            b_ensembl.append_value(&r.ensembl_id);
            b_symbol.append_value(&r.gene_symbol);
            b_chrom.append_value(&r.chromosome);
            b_start.append_value(r.start);
            b_end.append_value(r.end);
            b_nvariants.append_value(r.n_variants);
            b_cmac.append_value(r.cumulative_mac);

            let mut pi = 0;
            for p in [s.burden_1_25, s.burden_1_1, s.skat_1_25, s.skat_1_1, s.acat_v_1_25, s.acat_v_1_1] {
                pval_builders[pi].append_value(p);
                pi += 1;
            }
            for ann_p in &s.per_annotation {
                for &v in ann_p {
                    pval_builders[pi].append_value(v);
                    pi += 1;
                }
            }
            // Pad if fewer annotation channels than expected
            while pi < 6 + 6 * n_channels {
                pval_builders[pi].append_value(f64::NAN);
                pi += 1;
            }
            for p in [s.staar_b_1_25, s.staar_b_1_1, s.staar_s_1_25, s.staar_s_1_1, s.staar_a_1_25, s.staar_a_1_1, s.acat_o, s.staar_o] {
                pval_builders[pi].append_value(p);
                pi += 1;
            }
        }

        // Build schema
        let test_names = ["Burden(1,25)", "Burden(1,1)", "SKAT(1,25)", "SKAT(1,1)", "ACAT-V(1,25)", "ACAT-V(1,1)"];
        let mut fields = vec![
            Field::new("ensembl_id", DataType::Utf8, false),
            Field::new("gene_symbol", DataType::Utf8, false),
            Field::new("chromosome", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("end", DataType::UInt32, false),
            Field::new("n_variants", DataType::UInt32, false),
            Field::new("cMAC", DataType::UInt32, false),
        ];
        for test in &test_names {
            fields.push(Field::new(*test, DataType::Float64, true));
        }
        for ch in channels {
            for test in &test_names {
                fields.push(Field::new(format!("{test}-{ch}"), DataType::Float64, true));
            }
        }
        for name in ["STAAR-B(1,25)", "STAAR-B(1,1)", "STAAR-S(1,25)", "STAAR-S(1,1)", "STAAR-A(1,25)", "STAAR-A(1,1)", "ACAT-O", "STAAR-O"] {
            fields.push(Field::new(name, DataType::Float64, true));
        }
        let schema = Arc::new(Schema::new(fields));

        let mut columns: Vec<ArrayRef> = vec![
            Arc::new(b_ensembl.finish()), Arc::new(b_symbol.finish()),
            Arc::new(b_chrom.finish()), Arc::new(b_start.finish()),
            Arc::new(b_end.finish()), Arc::new(b_nvariants.finish()),
            Arc::new(b_cmac.finish()),
        ];
        for b in &mut pval_builders {
            columns.push(Arc::new(b.finish()));
        }

        let batch = RecordBatch::try_new(schema.clone(), columns)
            .map_err(|e| FavorError::Resource(format!("Arrow batch: {e}")))?;
        let file = File::create(&out_path)
            .map_err(|e| FavorError::Resource(format!("Create {}: {e}", out_path.display())))?;
        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(Default::default()))
            .build();
        let mut writer = ArrowWriter::try_new(file, schema, Some(props))
            .map_err(|e| FavorError::Resource(format!("Parquet writer: {e}")))?;
        writer.write(&batch).map_err(|e| FavorError::Resource(format!("Parquet write: {e}")))?;
        writer.close().map_err(|e| FavorError::Resource(format!("Parquet close: {e}")))?;

        let n_sig = results.iter().filter(|r| r.staar.staar_o < 2.5e-6).count();
        for r in results.iter().filter(|r| r.staar.staar_o < 2.5e-6) {
            significant_genes.push(json!({
                "gene": r.ensembl_id, "mask": mask_type.file_stem(),
                "STAAR-O": r.staar.staar_o, "n_variants": r.n_variants,
            }));
        }
        out.success(&format!("  {} -> {} genes, {} significant",
            mask_type.file_stem(), results.len(), n_sig));
    }

    let meta = json!({
        "favor_staar_version": 1,
        "traits": trait_names, "trait_type": format!("{:?}", trait_type),
        "n_samples": n, "n_rare_variants": n_rare, "maf_cutoff": maf_cutoff,
        "sigma2": null_model.sigma2, "significant_genes": significant_genes,
    });
    let _ = std::fs::write(output_dir.join("staar.meta.json"),
        serde_json::to_string_pretty(&meta).unwrap_or_default());

    match staar::summary::generate_report(all_mask_results, trait_names, n, n_rare, output_dir, "STAAR Rare Variant Association") {
        Ok(()) => out.success(&format!("  summary.html -> {}", output_dir.join("summary.html").display())),
        Err(e) => out.warn(&format!("  Summary report failed: {e}")),
    }

    out.success(&format!("STAAR complete -> {}", output_dir.display()));
    out.result_json(&meta);
    Ok(())
}
