use std::path::PathBuf;

use serde_json::json;

use crate::data::AnnotationDb;
use crate::commands;
use crate::config::{Config, Tier};
use crate::db::DuckEngine;
use crate::db::{query_scalar, query_strings};
use crate::error::FavorError;
use crate::ingest::JoinKey;
use crate::output::Output;
use crate::resource::Resources;
use crate::data::{VariantSet, VariantSetKind, VariantSetWriter};

const JOIN_KEY_COLUMNS: &[&str] = &["chromosome", "position", "ref", "alt"];

pub fn run(
    input: PathBuf,
    output_path: Option<PathBuf>,
    full: bool,
    out: &dyn Output,
    dry_run: bool,
) -> Result<(), FavorError> {
    if !input.exists() {
        return Err(FavorError::Input(format!("Not found: {}", input.display())));
    }

    let input_vs = VariantSet::open(&input)?;
    let join_key = input_vs.join_key();

    let config = Config::load_configured()?;
    let tier = if full { Tier::Full } else { config.data.tier };
    let ann_db = AnnotationDb::open_tier(&config, tier)?;

    let output_path = output_path.unwrap_or_else(|| {
        let name = input.file_name().unwrap_or_default().to_string_lossy();
        let stem = name.strip_suffix(".ingested").or_else(|| name.strip_suffix("/"))
            .unwrap_or(&name);
        input.parent().unwrap_or(&input).join(format!("{stem}.annotated"))
    });

    out.status(&format!("Join key: {:?}", join_key));

    let input_count = input_vs.count() as i64;
    out.status(&format!("Input: {} variants", input_count));

    if input_count == 0 {
        return Err(FavorError::Input("Input has 0 variants.".into()));
    }

    let input_struct_cols: Vec<&str> = input_vs.columns().iter()
        .map(|s| s.as_str())
        .filter(|c| !JOIN_KEY_COLUMNS.contains(c))
        .collect();

    out.status(&format!("Input: {} columns ({} in input struct)",
        input_vs.columns().len(), input_struct_cols.len()));

    if dry_run {
        let plan = commands::DryRunPlan {
            command: "annotate".into(),
            inputs: json!({
                "file": input.to_string_lossy(),
                "variant_count": input_count,
                "join_key": format!("{:?}", join_key),
                "tier": tier.as_str(),
            }),
            memory: commands::MemoryEstimate::duckdb_only(),
            output_path: output_path.to_string_lossy().into(),
        };
        commands::emit(&plan, out);
        return Ok(());
    }

    let resources = Resources::detect_with_config(&config.resources);
    let engine = DuckEngine::new(&resources)?;

    let input_struct_expr = if input_struct_cols.is_empty() {
        String::new()
    } else {
        let fields: Vec<String> = input_struct_cols.iter()
            .map(|c| format!("\"{c}\": u.\"{c}\""))
            .collect();
        format!("{{{}}}", fields.join(", "))
    };

    engine.execute("SET preserve_insertion_order = false;")?;
    let gb: u64 = 1024 * 1024 * 1024;
    let temp_limit_gb = (resources.memory_bytes * 3 / gb).max(50);
    engine.execute(&format!("SET max_temp_directory_size = '{}GiB';", temp_limit_gb))?;

    let input_chroms = input_vs.chromosomes();
    out.status(&format!("Chromosomes: {}", input_chroms.join(", ")));
    out.status(&format!("Annotating against favor-{tier} ({})...", tier.size_human()));

    let mut out_writer = VariantSetWriter::new(&output_path, join_key,
        &input.display().to_string())?;
    out_writer.set_kind(VariantSetKind::Annotated { tier: tier.as_str().to_string() });

    match join_key {
        JoinKey::ChromPosRefAlt | JoinKey::ChromPos => {
            let exact_allele = join_key == JoinKey::ChromPosRefAlt;
            let total_chroms = input_chroms.len();

            for (i, chrom) in input_chroms.iter().enumerate() {
                let chrom_parquet = match ann_db.chrom_parquet(chrom) {
                    Some(p) => p,
                    None => {
                        out.warn(&format!("  chr{chrom}: skipping (no annotation file)"));
                        continue;
                    }
                };

                out.status(&format!("  chr{chrom} ({}/{})", i + 1, total_chroms));

                let out_path = out_writer.chrom_path(chrom)?;
                let input_read = input_vs.read_chrom(chrom);
                let sql = build_chrom_query(
                    &input_read, &chrom_parquet, &out_path,
                    &input_struct_expr, chrom, exact_allele,
                );
                engine.execute(&sql)?;

                let count = query_scalar(&engine, &format!(
                    "SELECT COUNT(*) FROM read_parquet('{}')", out_path.display()
                )).unwrap_or(0) as u64;
                let size = std::fs::metadata(&out_path).map_or(0, |m| m.len());
                if count > 0 {
                    out_writer.register_chrom(chrom, count, size);
                } else {
                    let _ = std::fs::remove_file(&out_path);
                }
            }
        }
        JoinKey::Rsid => {
            out.status("  Resolving rsids across all chromosomes...");
            // Rsid is cross-chromosome — load full input into temp table
            engine.execute(&format!(
                "CREATE TEMP TABLE _user AS SELECT * FROM {}",
                input_vs.read_all(),
            ))?;

            let annotation_glob = format!(
                "{}/chromosome=*/sorted.parquet", ann_db.root().display()
            );
            let rsid_sql = build_rsid_query(
                &annotation_glob, &output_path, &config, &input_struct_expr,
            );
            engine.execute(&rsid_sql)?;
            let _ = engine.execute("DROP TABLE IF EXISTS _user");

            // rsid query wrote via PARTITION_BY — scan to register
            out_writer.scan_and_register(&engine)?;
        }
    }

    let output_vs = out_writer.finish()?;
    let output_count = output_vs.count() as i64;

    let match_rate = if input_count > 0 {
        output_count as f64 / input_count as f64
    } else {
        0.0
    };

    if output_count == 0 {
        let diagnostic = diagnose_zero_matches(&engine, &input_vs, ann_db.root())?;
        out.warn("0 variants annotated.");
        out.warn(&diagnostic.message);
        out.result_json(&json!({
            "status": "no_matches",
            "input_count": input_count,
            "output_count": 0,
            "diagnostic": diagnostic.message,
            "suggested_fix": diagnostic.suggested_fix,
        }));
        return Err(FavorError::Analysis(
            format!("0 variants matched annotations. {}", diagnostic.message)
        ));
    }

    if match_rate < 0.05 {
        let diagnostic = diagnose_zero_matches(&engine, &input_vs, ann_db.root())?;
        out.warn(&format!(
            "Low match rate: {}/{} ({:.1}%). {}",
            output_count, input_count, match_rate * 100.0, diagnostic.message,
        ));
        out.result_json(&json!({
            "status": "low_match_rate",
            "input_count": input_count,
            "output_count": output_count,
            "match_rate": match_rate,
            "diagnostic": diagnostic.message,
            "suggested_fix": diagnostic.suggested_fix,
            "output": output_vs.root().to_string_lossy(),
        }));
        return Ok(());
    }

    let unmatched = input_count - output_count;
    out.success(&format!(
        "Annotated {output_count}/{input_count} variants ({:.1}%) → {}",
        match_rate * 100.0,
        output_vs.root().display(),
    ));
    if unmatched > 0 {
        out.status(&format!("  {unmatched} variants not found in annotations"));
    }

    out.result_json(&json!({
        "status": "ok",
        "output": output_vs.root().to_string_lossy(),
        "input_count": input_count,
        "annotated_count": output_count,
        "match_rate": match_rate,
        "tier": tier.as_str(),
    }));

    Ok(())
}

/// Per-chromosome annotation query. Reads directly from the input partition
/// (no global temp table). Two-pass approach for memory efficiency:
///   Pass 1: SEMI JOIN on position + alleles to find matching vids
///   Pass 2: Read only matched rows from the annotation parquet
fn build_chrom_query(
    input_read_expr: &str,
    chrom_parquet: &std::path::Path,
    output: &std::path::Path,
    input_struct_expr: &str,
    chrom: &str,
    exact_allele: bool,
) -> String {
    let input_select = if input_struct_expr.is_empty() {
        String::new()
    } else {
        format!("{input_struct_expr} AS input, ")
    };

    let allele_condition = if exact_allele {
        "AND u.ref = a.ref_vcf AND u.alt = a.alt_vcf "
    } else {
        ""
    };

    format!(
        "CREATE OR REPLACE TEMP TABLE _input_{chrom} AS \
         SELECT * FROM {input_read_expr}; \
         \
         CREATE OR REPLACE TEMP TABLE _matched_{chrom} AS \
         SELECT a.vid \
         FROM _input_{chrom} u \
         INNER JOIN (\
             SELECT vid, position, ref_vcf, alt_vcf \
             FROM read_parquet('{chrom_parquet}')\
         ) a ON u.position = a.position {allele_condition}; \
         \
         COPY (\
             SELECT \
                 {input_select}\
                 a.* \
             FROM _input_{chrom} u \
             INNER JOIN (\
                 SELECT * FROM read_parquet('{chrom_parquet}') \
                 WHERE vid IN (SELECT vid FROM _matched_{chrom})\
             ) a ON u.position = a.position {allele_condition}\
         ) TO '{output}' (FORMAT PARQUET, COMPRESSION ZSTD); \
         \
         DROP TABLE IF EXISTS _matched_{chrom}; \
         DROP TABLE IF EXISTS _input_{chrom};",
        input_read_expr = input_read_expr,
        input_select = input_select,
        chrom_parquet = chrom_parquet.display(),
        chrom = chrom,
        allele_condition = allele_condition,
        output = output.display(),
    )
}

/// Rsid join: resolve rsids via lookup, then write partitioned output.
fn build_rsid_query(
    annotation_glob: &str,
    output_dir: &std::path::Path,
    config: &Config,
    input_struct_expr: &str,
) -> String {
    let lookup_glob = format!(
        "{}/lookup/rsid_lookup/chromosome=*/*.parquet",
        config.root_dir().display()
    );

    let input_select = if input_struct_expr.is_empty() {
        String::new()
    } else {
        format!("{input_struct_expr} AS input, ")
    };

    format!(
        "CREATE TEMP TABLE _rsid_resolved AS \
         SELECT u.*, lk.chromosome AS _chr, lk.position AS _pos, \
                lk.ref_vcf AS _ref, lk.alt_vcf AS _alt \
         FROM _user u \
         INNER JOIN read_parquet('{lookup_glob}', hive_partitioning=true) lk \
             ON u.rsid = lk.rsid; \
         \
         COPY (\
             SELECT \
                 {input_select}\
                 a.* \
             FROM _rsid_resolved r \
             INNER JOIN read_parquet('{annotation_glob}', hive_partitioning=true) a \
                 ON r._chr = a.chromosome \
                 AND r._pos = a.position \
                 AND r._ref = a.ref_vcf \
                 AND r._alt = a.alt_vcf\
         ) TO '{output_dir}' (FORMAT PARQUET, PARTITION_BY (chromosome), COMPRESSION ZSTD); \
         \
         DROP TABLE IF EXISTS _rsid_resolved;",
        input_select = input_select,
        lookup_glob = lookup_glob,
        annotation_glob = annotation_glob,
        output_dir = output_dir.display(),
    )
}

struct Diagnostic {
    message: String,
    suggested_fix: String,
}

fn diagnose_zero_matches(
    engine: &DuckEngine,
    input_vs: &VariantSet,
    annotations_dir: &std::path::Path,
) -> Result<Diagnostic, FavorError> {
    let chr1_path = annotations_dir.join("chromosome=1/sorted.parquet");
    if !chr1_path.exists() {
        return Ok(Diagnostic {
            message: "Annotations incomplete — chromosome=1 not found.".into(),
            suggested_fix: "favor data pull".into(),
        });
    }

    // Load chr1 input for probing
    let chr1_input = input_vs.read_chrom("1");
    let chroms = query_strings(engine, &format!(
        "SELECT DISTINCT chromosome FROM {} LIMIT 5", input_vs.read_all()
    )).unwrap_or_default();

    if chroms.iter().any(|c| c.starts_with("chr")) {
        return Ok(Diagnostic {
            message: "Chromosome names have 'chr' prefix — annotations use bare numbers. Re-ingest.".into(),
            suggested_fix: "favor ingest <original_file>".into(),
        });
    }

    let position_hits = query_scalar(engine, &format!(
        "SELECT COUNT(*) FROM (\
            SELECT position FROM {chr1_input} LIMIT 10\
        ) s SEMI JOIN read_parquet('{}') a ON s.position = a.position",
        chr1_path.display()
    )).unwrap_or(0);

    if position_hits == 0 {
        return Ok(Diagnostic {
            message: "No positions match hg38 annotations. Input may be hg19 or different build.".into(),
            suggested_fix: "favor ingest <original_file> --build hg19 --emit-sql".into(),
        });
    }

    Ok(Diagnostic {
        message: "Positions exist but alleles don't match. Check ref/alt normalization.".into(),
        suggested_fix: "favor ingest <original_file> --emit-sql".into(),
    })
}
