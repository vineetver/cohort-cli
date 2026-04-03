use std::path::PathBuf;

use serde_json::json;

use crate::cli::GenomeBuild;
use crate::commands;
use crate::config::Config;
use crate::db::DuckEngine;
use crate::error::FavorError;
use crate::ingest::{self, BuildGuess, InputFormat};
use crate::output::Output;
use crate::resource::Resources;
use crate::data::{VariantSetKind, VariantSetWriter};

pub fn run(
    input: PathBuf,
    output: Option<PathBuf>,
    emit_sql: bool,
    build_override: Option<GenomeBuild>,
    out: &dyn Output,
    dry_run: bool,
) -> Result<(), FavorError> {
    if !input.exists() {
        return Err(FavorError::Input(format!(
            "File not found: {}",
            input.display()
        )));
    }

    let output_path = output.unwrap_or_else(|| {
        let stem = input.file_stem().unwrap_or_default().to_string_lossy();
        let stem = stem.strip_suffix(".vcf").or_else(|| stem.strip_suffix(".tsv"))
            .or_else(|| stem.strip_suffix(".csv"))
            .unwrap_or(&stem);
        input.with_file_name(format!("{stem}.ingested"))
    });

    out.status(&format!("Analyzing {}...", input.display()));
    let mut analysis = ingest::analyze(&input)?;

    out.status(&format!("  Format: {:?}", analysis.format));
    out.status(&format!("  Columns: {}", analysis.raw_columns.join(", ")));
    out.status(&format!("  Join key: {:?}", analysis.join_key));

    for mapping in &analysis.columns {
        if mapping.canonical != mapping.input_name.to_lowercase() {
            out.status(&format!("  {} → {}", mapping.input_name, mapping.canonical));
        }
    }
    for amb in &analysis.ambiguous {
        out.warn(&format!("  {} — {}", amb.column, amb.reason));
    }

    match build_override {
        Some(GenomeBuild::Hg19) => {
            analysis.build_guess = BuildGuess::Hg19 {
                match_rate_hg38: 0.0,
                match_rate_hg19: 1.0,
            };
        }
        Some(GenomeBuild::Hg38) => {
            analysis.build_guess = BuildGuess::Hg38;
        }
        None if analysis.format != InputFormat::Vcf => {
            if let Ok(config) = Config::load() {
                if !config.data.root_dir.is_empty() {
                    let resources = Resources::detect_with_config(&config.resources);
                    if let Ok(engine) = DuckEngine::new(&resources) {
                        let _ = ingest::detect::detect_build_and_coords(
                            &mut analysis, &input, &engine, &config,
                        );
                    }
                }
            }
        }
        None => {}
    }

    match &analysis.build_guess {
        BuildGuess::Hg38 => out.status("  Build: hg38"),
        BuildGuess::Hg19 { .. } => out.warn("  Build: likely hg19 — needs liftover"),
        BuildGuess::Unknown => out.status("  Build: unknown (no annotations to probe)"),
    }

    match analysis.coord_base {
        ingest::CoordBase::OneBased => out.status("  Coordinates: 1-based"),
        ingest::CoordBase::ZeroBased => out.warn("  Coordinates: 0-based (will convert to 1-based)"),
        ingest::CoordBase::Unknown => {}
    }

    if dry_run {
        let plan = commands::DryRunPlan {
            command: "ingest".into(),
            inputs: json!({
                "file": input.to_string_lossy(),
                "file_size": commands::file_size(&input),
                "format": format!("{:?}", analysis.format),
                "join_key": format!("{:?}", analysis.join_key),
                "build": format!("{:?}", analysis.build_guess),
                "columns_mapped": analysis.columns.len(),
                "columns_ambiguous": analysis.ambiguous.len(),
                "needs_intervention": analysis.needs_intervention(),
            }),
            memory: commands::MemoryEstimate::duckdb_only(),
            output_path: output_path.to_string_lossy().into(),
        };
        commands::emit(&plan, out);
        return Ok(());
    }

    let config_resources = Config::load()
        .map(|c| c.resources)
        .unwrap_or_default();

    if analysis.format == InputFormat::Vcf {
        if emit_sql {
            let resources = Resources::detect_with_config(&config_resources);
            let sql = ingest::sql::generate_sql(&analysis, &input, &output_path, &resources);
            let script_path = output_path.with_extension("ingest.sql");
            std::fs::write(&script_path, &sql)?;
            out.status(&format!("VCF hint script: {}", script_path.display()));
            return Ok(());
        }

        let resources = Resources::detect_with_config(&config_resources);
        out.status("Ingesting VCF (streaming, per-chromosome writers)...");

        let mut vs_writer = VariantSetWriter::new(
            &output_path, ingest::JoinKey::ChromPosRefAlt,
            &input.display().to_string(),
        )?;
        vs_writer.set_kind(VariantSetKind::Ingested);
        let result = ingest::vcf::ingest_vcf(&input, &mut vs_writer, resources.memory_bytes, out)?;
        let vs = vs_writer.finish()?;

        out.success(&format!("Ingested {} variants → {}", result.variant_count, vs.root().display()));
        if result.filtered_contigs > 0 {
            out.status(&format!("  {} records on non-standard contigs filtered", result.filtered_contigs));
        }
        if result.multiallelic_split > 0 {
            out.status(&format!("  {} multi-allelic sites split to biallelic", result.multiallelic_split));
        }

        out.result_json(&json!({
            "status": "ok",
            "output": vs.root().to_string_lossy(),
            "variant_count": result.variant_count,
            "join_key": "chrom_pos_ref_alt",
        }));

        return Ok(());
    }

    let resources = Resources::detect_with_config(&config_resources);
    let sql = ingest::sql::generate_sql(&analysis, &input, &output_path, &resources);

    if emit_sql || analysis.needs_intervention() {
        let script_path = output_path.with_extension("ingest.sql");
        std::fs::write(&script_path, &sql)?;

        if analysis.needs_intervention() {
            out.warn("Cannot auto-ingest — ambiguities detected. Edit the script and run:");
            out.warn(&format!("  duckdb < {}", script_path.display()));
        } else {
            out.status(&format!("SQL script written to: {}", script_path.display()));
            out.status(&format!("  Run: duckdb < {}", script_path.display()));
        }

        out.result_json(&json!({
            "status": analysis.status(),
            "format": analysis.format,
            "join_key": analysis.join_key,
            "build": analysis.build_guess,
            "coord_base": analysis.coord_base,
            "columns_mapped": analysis.columns.len(),
            "columns_ambiguous": analysis.ambiguous.len(),
            "columns_unmapped": analysis.unmapped.len(),
            "script_path": script_path.to_string_lossy(),
            "suggested_command": if analysis.needs_intervention() {
                format!("duckdb < {}", script_path.display())
            } else {
                format!("favor ingest {}", input.display())
            },
        }));

        return Ok(());
    }

    out.status("Ingesting...");

    let mut vs_writer = VariantSetWriter::new(
        &output_path, analysis.join_key,
        &input.display().to_string(),
    )?;
    vs_writer.set_kind(VariantSetKind::Ingested);

    // SQL with PARTITION_BY writes directly to the output directory
    let engine = DuckEngine::new(&resources)?;
    engine.execute(&sql)?;

    vs_writer.scan_and_register(&engine)?;
    let vs = vs_writer.finish()?;

    out.success(&format!("Ingested → {}", vs.root().display()));
    out.status(&format!("  Join key: {:?}", analysis.join_key));

    out.result_json(&json!({
        "status": "ok",
        "output": vs.root().to_string_lossy(),
        "join_key": analysis.join_key,
        "build": analysis.build_guess,
    }));

    Ok(())
}
