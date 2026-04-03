use std::path::{Path, PathBuf};

use serde_json::json;

use crate::commands::dry_run;
use crate::config::Config;
use crate::db::engine::DuckEngine;
use crate::db::query::{query_scalar, query_strings};
use crate::error::FavorError;
use crate::output::Output;
use crate::resource::Resources;
use crate::variant_set::VariantSet;

struct EnrichTable {
    dir: &'static str,
    name: &'static str,
    has_tissue: bool,
}

const VARIANT_TABLES: &[EnrichTable] = &[
    EnrichTable { dir: "variant_eqtl",               name: "eqtl",                has_tissue: true },
    EnrichTable { dir: "variant_sqtl",               name: "sqtl",                has_tissue: true },
    EnrichTable { dir: "variant_apaqtl",             name: "apaqtl",              has_tissue: true },
    EnrichTable { dir: "variant_eqtl_susie",         name: "eqtl_susie",          has_tissue: true },
    EnrichTable { dir: "variant_eqtl_catalogue",     name: "eqtl_catalogue",      has_tissue: true },
    EnrichTable { dir: "variant_sc_eqtl",            name: "sc_eqtl",             has_tissue: true },
    EnrichTable { dir: "variant_sc_eqtl_dice",       name: "sc_eqtl_dice",        has_tissue: true },
    EnrichTable { dir: "variant_sc_eqtl_psychencode",name: "sc_eqtl_psychencode", has_tissue: true },
    EnrichTable { dir: "variant_tissue_scores",      name: "tissue_scores",       has_tissue: true },
    EnrichTable { dir: "variant_chrombpnet",         name: "chrombpnet",           has_tissue: true },
    EnrichTable { dir: "variant_allelic_imbalance",  name: "allelic_imbalance",    has_tissue: true },
    EnrichTable { dir: "variant_allelic_methylation", name: "allelic_methylation",  has_tissue: true },
    EnrichTable { dir: "variant_eqtl_ccre",          name: "eqtl_ccre",           has_tissue: true },
];

pub fn run(
    input: PathBuf,
    tissue: String,
    output_path: Option<PathBuf>,
    out: &dyn Output,
    dry_run: bool,
) -> Result<(), FavorError> {
    if !input.exists() {
        return Err(FavorError::Input(format!("Not found: {}", input.display())));
    }

    let input_vs = VariantSet::open(&input)?;
    if !input_vs.has_column("vid") {
        return Err(FavorError::Input(
            "Input has no 'vid' column. Run `favor annotate` first.".into()
        ));
    }
    let input_count = input_vs.count() as i64;

    let config = Config::load_configured()?;
    let tissue_dir = config.tissue_dir();

    if !tissue_dir.exists() {
        return Err(FavorError::DataMissing(format!(
            "Tissue data not found at {}. Run `favor data pull --pack eqtl` first.",
            tissue_dir.display()
        )));
    }

    let output_dir = output_path.unwrap_or_else(|| {
        let name = input.file_name().unwrap_or_default().to_string_lossy();
        let stem = name.strip_suffix(".annotated").or_else(|| name.strip_suffix("/"))
            .unwrap_or(&name);
        input.parent().unwrap_or(&input).join(format!("{stem}.enriched"))
    });
    std::fs::create_dir_all(&output_dir)?;

    let resources = Resources::detect_with_config(&config.resources);
    let engine = DuckEngine::new(&resources)?;

    out.status(&format!("Input: {} variants", input_count));

    let tissue_vocab_path = tissue_dir.join("reference/tissue_vocab.parquet");
    let resolved_tissues = resolve_tissue(&engine, &tissue_vocab_path, &tissue)?;

    if resolved_tissues.is_empty() {
        let groups = list_tissue_groups(&engine, &tissue_vocab_path)?;
        return Err(FavorError::Input(format!(
            "Unknown tissue '{}'. Available groups: {}",
            tissue, groups.join(", ")
        )));
    }

    out.status(&format!("Tissue '{}' → {} subtissues", tissue, resolved_tissues.len()));

    if dry_run {
        let available_tables: Vec<&str> = VARIANT_TABLES.iter()
            .filter(|t| tissue_dir.join(t.dir).is_dir())
            .map(|t| t.name)
            .collect();
        let plan = dry_run::DryRunPlan {
            command: "enrich".into(),
            inputs: json!({
                "file": input.to_string_lossy(),
                "variant_count": input_count,
                "tissue": tissue,
                "resolved_tissues": resolved_tissues,
                "available_tables": available_tables,
            }),
            memory: dry_run::MemoryEstimate::duckdb_only(),
            output_path: output_dir.to_string_lossy().into(),
        };
        dry_run::emit(&plan, out);
        return Ok(());
    }

    let tissue_filter: String = resolved_tissues.iter()
        .map(|t| format!("'{}'", t.replace('\'', "''")))
        .collect::<Vec<_>>()
        .join(", ");

    let read_all = input_vs.read_all();

    let annotated_out = output_dir.join("annotated.parquet");
    engine.execute(&format!(
        "COPY (SELECT * FROM {read_all}) TO '{}' (FORMAT PARQUET, COMPRESSION ZSTD);",
        annotated_out.display(),
    ))?;
    out.status(&format!("  annotated.parquet ({} variants)", input_count));

    engine.execute(&format!(
        "CREATE TEMP TABLE _input_vids AS SELECT DISTINCT vid FROM {read_all}",
    ))?;

    let mut tables_written: Vec<(String, i64)> = Vec::new();

    for table in VARIANT_TABLES {
        let table_path = tissue_dir.join(table.dir);
        if !table_path.is_dir() { continue; }

        let glob = format!("{}/chrom_id=*/data_0.parquet", table_path.display());
        let out_path = output_dir.join(format!("{}.parquet", table.name));

        let where_clause = if table.has_tissue {
            format!("WHERE t.vid IN (SELECT vid FROM _input_vids) \
                     AND t.tissue_name IN ({tissue_filter})")
        } else {
            "WHERE t.vid IN (SELECT vid FROM _input_vids)".to_string()
        };

        let sql = format!(
            "COPY (\
                SELECT t.* EXCLUDE (chrom_id) \
                FROM read_parquet('{glob}', hive_partitioning=true) t \
                {where_clause}\
            ) TO '{out_path}' (FORMAT PARQUET, COMPRESSION ZSTD);",
            glob = glob,
            where_clause = where_clause,
            out_path = out_path.display(),
        );

        out.status(&format!("  {}: joining...", table.name));
        if let Err(e) = engine.execute(&sql) {
            out.warn(&format!("  {}: skipped ({e})", table.name));
            continue;
        }

        let row_count = query_scalar(&engine, &format!(
            "SELECT COUNT(*) FROM read_parquet('{}')", out_path.display()
        )).unwrap_or(0);

        if row_count == 0 {
            let _ = std::fs::remove_file(&out_path);
        } else {
            out.status(&format!("  {}.parquet ({} rows)", table.name, row_count));
            tables_written.push((table.name.to_string(), row_count));
        }
    }

    let _ = engine.execute("DROP TABLE IF EXISTS _input_vids");

    if tables_written.is_empty() {
        out.warn("No enrichment data found for these variants in this tissue.");
    } else {
        out.success(&format!(
            "Enriched → {} ({} tables)",
            output_dir.display(), tables_written.len(),
        ));
    }

    let meta = json!({
        "favor_enrich_version": 2,
        "source": input.to_string_lossy(),
        "tissue": tissue,
        "resolved_tissues": resolved_tissues,
        "output_dir": output_dir.to_string_lossy(),
        "tables": tables_written.iter()
            .map(|(name, rows)| json!({"name": name, "rows": rows}))
            .collect::<Vec<_>>(),
        "input_count": input_count,
        "join_key": "vid",
        "usage": "SELECT a.*, e.* FROM 'annotated.parquet' a INNER JOIN 'eqtl.parquet' e ON a.vid = e.vid",
    });
    let meta_path = output_dir.join("enriched.meta.json");
    if let Ok(json_str) = serde_json::to_string_pretty(&meta) {
        let _ = std::fs::write(&meta_path, json_str);
    }

    out.result_json(&json!({
        "status": "ok",
        "output_dir": output_dir.to_string_lossy(),
        "tables": tables_written.iter()
            .map(|(name, rows)| json!({"name": name, "rows": rows}))
            .collect::<Vec<_>>(),
        "tissue": tissue,
        "input_count": input_count,
    }));

    Ok(())
}

fn resolve_tissue(
    engine: &DuckEngine,
    tissue_vocab_path: &Path,
    tissue_query: &str,
) -> Result<Vec<String>, FavorError> {
    if !tissue_vocab_path.exists() {
        return Ok(vec![tissue_query.to_string()]);
    }

    let escaped = tissue_query.replace('\'', "''");

    let sql = format!(
        "SELECT DISTINCT tissue_norm FROM read_parquet('{}') \
         WHERE tissue_group ILIKE '%{escaped}%' \
            OR tissue_norm ILIKE '%{escaped}%' \
            OR tissue_raw ILIKE '%{escaped}%'",
        tissue_vocab_path.display(),
    );

    query_strings(engine, &sql)
}

fn list_tissue_groups(
    engine: &DuckEngine,
    tissue_vocab_path: &Path,
) -> Result<Vec<String>, FavorError> {
    if !tissue_vocab_path.exists() {
        return Ok(vec!["(tissue_vocab.parquet not found)".into()]);
    }

    query_strings(engine, &format!(
        "SELECT DISTINCT tissue_group FROM read_parquet('{}') ORDER BY tissue_group",
        tissue_vocab_path.display(),
    ))
}
