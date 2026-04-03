//! SQL generator for ingest scripts.
//!
//! Takes an Analysis and produces a complete, runnable DuckDB SQL string.
//! Ambiguities become `-- FIXME:` comments that agents or humans can edit.
//! Rules enforced: no row counts, all_varchar=true, explicit CAST, streaming COPY.

use std::fmt::Write;
use std::path::Path;

use super::{Analysis, BuildGuess, CoordBase, Delimiter, InputFormat};
use super::canonical_type;
use crate::resource::Resources;

/// Generate a complete DuckDB SQL script from an Analysis.
/// If the analysis has ambiguities, FIXME comments are injected.
pub fn generate_sql(
    analysis: &Analysis,
    input_path: &Path,
    output_path: &Path,
    resources: &Resources,
) -> String {
    let mut sql = String::with_capacity(4096);

    write_header(&mut sql, analysis, input_path, output_path);
    write_pragmas(&mut sql, resources);

    match analysis.format {
        InputFormat::Tabular => write_tabular_query(&mut sql, analysis, input_path, output_path),
        InputFormat::Parquet => write_parquet_query(&mut sql, analysis, input_path, output_path),
        InputFormat::Vcf => write_vcf_stub(&mut sql, input_path, output_path),
    }

    sql
}

fn write_header(sql: &mut String, analysis: &Analysis, input_path: &Path, output_path: &Path) {
    let _ = writeln!(sql, "-- FAVOR ingest: {} -> {}", input_path.display(), output_path.display());
    let _ = writeln!(sql, "-- Format: {:?} | Join key: {:?}", analysis.format, analysis.join_key);

    if analysis.needs_intervention() {
        let _ = writeln!(sql, "-- STATUS: needs_edit — resolve FIXME comments before running");
    } else {
        let _ = writeln!(sql, "-- STATUS: ok — ready to run as-is");
    }

    let _ = writeln!(sql, "--");
    let _ = writeln!(sql, "-- Rules (do not violate):");
    let _ = writeln!(sql, "--   1. NEVER COUNT(*) or aggregate the input file (expensive on CSV)");
    let _ = writeln!(sql, "--   2. ALWAYS read_csv with all_varchar=true (no type guessing)");
    let _ = writeln!(sql, "--   3. ALWAYS explicit CAST on every output column");
    let _ = writeln!(sql, "--   4. STREAM via COPY (SELECT FROM read_csv) TO — never CREATE TABLE");
    let _ = writeln!(sql, "--   5. ZSTD compression, always");
    let _ = writeln!(sql, "--   6. Filter non-standard chromosomes in WHERE");
    let _ = writeln!(sql, "--");
    let _ = writeln!(sql);
}

fn write_pragmas(sql: &mut String, resources: &Resources) {
    let _ = writeln!(sql, "PRAGMA temp_directory = '{}';", resources.temp_dir.display());
    let _ = writeln!(sql, "PRAGMA threads = {};", resources.threads);
    let _ = writeln!(sql, "PRAGMA memory_limit = '{}';", resources.duckdb_memory());
    let _ = writeln!(sql);
}

fn write_tabular_query(
    sql: &mut String,
    analysis: &Analysis,
    input_path: &Path,
    output_path: &Path,
) {
    let delim = analysis.delimiter.unwrap_or(Delimiter::Tab);

    let _ = writeln!(sql, "COPY (");
    let _ = writeln!(sql, "    SELECT");

    if let Some(chr_col) = &analysis.chr_col {
        let _ = writeln!(sql, "        -- Chromosome: strip 'chr' prefix, normalize M→MT");
        let _ = writeln!(sql, "        CASE UPPER(REGEXP_REPLACE(CAST(\"{chr_col}\" AS VARCHAR), '^chr', '', 'i'))");
        let _ = writeln!(sql, "            WHEN 'M' THEN 'MT'");
        let _ = writeln!(sql, "            WHEN '23' THEN 'X'");
        let _ = writeln!(sql, "            WHEN '24' THEN 'Y'");
        let _ = writeln!(sql, "            ELSE UPPER(REGEXP_REPLACE(CAST(\"{chr_col}\" AS VARCHAR), '^chr', '', 'i'))");
        let _ = writeln!(sql, "        END AS chromosome,");
    } else {
        let _ = writeln!(sql, "        -- FIXME: no chromosome column detected. Add: CAST(\"<your_chr_col>\" AS VARCHAR) AS chromosome,");
    }

    if let Some(pos_col) = &analysis.pos_col {
        match analysis.coord_base {
            CoordBase::ZeroBased => {
                let _ = writeln!(sql, "        -- Position: detected 0-based, converting to 1-based");
                let _ = writeln!(sql, "        CAST(\"{pos_col}\" AS INTEGER) + 1 AS position,");
            }
            _ => {
                let _ = writeln!(sql, "        CAST(\"{pos_col}\" AS INTEGER) AS position,");
            }
        }
    } else {
        let _ = writeln!(sql, "        -- FIXME: no position column detected. Add: CAST(\"<your_pos_col>\" AS INTEGER) AS position,");
    }

    if let Some(ref_col) = &analysis.ref_col {
        let _ = writeln!(sql, "        UPPER(TRIM(CAST(\"{ref_col}\" AS VARCHAR))) AS ref,");
    }
    if let Some(alt_col) = &analysis.alt_col {
        let _ = writeln!(sql, "        UPPER(TRIM(CAST(\"{alt_col}\" AS VARCHAR))) AS alt,");
    }

    // Handle ambiguous A1/A2 columns
    for amb in &analysis.ambiguous {
        let _ = writeln!(sql, "        -- FIXME: '{}' is ambiguous — {}", amb.column, amb.reason);
        let _ = writeln!(sql, "        -- Most common convention: A1=effect(alt), A2=other(ref).");
        let _ = writeln!(sql, "        -- If your study reverses this, swap the AS labels below:");
        if analysis.ref_col.is_none() && analysis.alt_col.is_none() {
            // No ref/alt mapped yet — fill from ambiguous columns
            if amb.column.to_lowercase() == "a1" || amb.column.to_lowercase() == "allele1" {
                let _ = writeln!(sql, "        UPPER(TRIM(CAST(\"{}\" AS VARCHAR))) AS alt,  -- assumed effect allele", amb.column);
            } else {
                let _ = writeln!(sql, "        UPPER(TRIM(CAST(\"{}\" AS VARCHAR))) AS ref,  -- assumed other allele", amb.column);
            }
        }
    }

    if analysis.ref_col.is_none() && analysis.alt_col.is_none() && analysis.ambiguous.is_empty() {
        let _ = writeln!(sql, "        -- FIXME: no ref/alt columns detected.");
        let _ = writeln!(sql, "        -- Add: UPPER(TRIM(CAST(\"<ref_col>\" AS VARCHAR))) AS ref,");
        let _ = writeln!(sql, "        -- Add: UPPER(TRIM(CAST(\"<alt_col>\" AS VARCHAR))) AS alt,");
    }

    if let Some(rsid_col) = &analysis.rsid_col {
        let _ = writeln!(sql, "        CAST(\"{rsid_col}\" AS VARCHAR) AS rsid,");
    }

    for mapping in &analysis.columns {
        // Skip the core columns we already handled
        if ["chromosome", "position", "ref", "alt", "rsid"].contains(&mapping.canonical) {
            continue;
        }

        let typ = canonical_type(mapping.canonical);

        // -log10(p): keep as-is, don't convert. Downstream can handle either representation.
        if mapping.canonical == "neglog10p" {
            let _ = writeln!(sql, "        CAST(\"{}\" AS DOUBLE) AS neglog10p,", mapping.input_name);
            continue;
        }

        let _ = writeln!(sql, "        CAST(\"{}\" AS {typ}) AS {},", mapping.input_name, mapping.canonical);
    }

    for col in &analysis.unmapped {
        let _ = writeln!(sql, "        \"{col}\",");
    }

    // Remove trailing comma from last column
    // (DuckDB tolerates trailing commas in SELECT, so this is cosmetic)

    let _ = writeln!(sql, "    FROM read_csv(");
    let _ = writeln!(sql, "        '{}',", input_path.display());
    let _ = writeln!(sql, "        header = true,");
    let _ = writeln!(sql, "        all_varchar = true,");
    let _ = writeln!(sql, "        delim = {}", delim.sql_literal());
    let _ = writeln!(sql, "    )");

    if let Some(chr_col) = &analysis.chr_col {
        let _ = writeln!(sql, "    WHERE UPPER(REGEXP_REPLACE(CAST(\"{}\" AS VARCHAR), '^chr', '', 'i')) IN (",
            chr_col);
        let _ = writeln!(sql, "        '1','2','3','4','5','6','7','8','9','10',");
        let _ = writeln!(sql, "        '11','12','13','14','15','16','17','18','19','20',");
        let _ = writeln!(sql, "        '21','22','X','Y','MT','M','23','24'");
        let _ = writeln!(sql, "    )");
    }

    let _ = writeln!(sql, "    ORDER BY position");

    let _ = writeln!(sql, ") TO '{}' (FORMAT PARQUET, PARTITION_BY (chromosome), COMPRESSION ZSTD);", output_path.display());

    match &analysis.build_guess {
        BuildGuess::Hg19 { match_rate_hg38, match_rate_hg19 } => {
            let _ = writeln!(sql);
            let _ = writeln!(sql, "-- WARNING: Input appears to be hg19 (hg38 match: {:.0}%, hg19 match: {:.0}%)",
                match_rate_hg38 * 100.0, match_rate_hg19 * 100.0);
            let _ = writeln!(sql, "-- FIXME: This data needs liftover to hg38 before annotation.");
            let _ = writeln!(sql, "-- Use UCSC liftOver: liftOver input.bed hg19ToHg38.over.chain.gz output.bed unmapped.bed");
        }
        _ => {}
    }
}

fn write_parquet_query(
    sql: &mut String,
    _analysis: &Analysis,
    input_path: &Path,
    output_path: &Path,
) {
    let _ = writeln!(sql, "-- Parquet re-normalization: inspect schema and rename columns");
    let _ = writeln!(sql, "-- FIXME: verify column names match your parquet schema");
    let _ = writeln!(sql, "COPY (");
    let _ = writeln!(sql, "    SELECT *");
    let _ = writeln!(sql, "    FROM read_parquet('{}')  ", input_path.display());
    let _ = writeln!(sql, ") TO '{}' (FORMAT PARQUET, PARTITION_BY (chromosome), COMPRESSION ZSTD);", output_path.display());
}

fn write_vcf_stub(
    sql: &mut String,
    input_path: &Path,
    output_path: &Path,
) {
    let _ = writeln!(sql, "-- VCF ingestion is not yet supported via SQL generation.");
    let _ = writeln!(sql, "-- FIXME: Convert VCF to TSV first, then re-run ingest:");
    let _ = writeln!(sql, "--   bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n' {} > variants.tsv", input_path.display());
    let _ = writeln!(sql, "--   favor ingest variants.tsv -o {}", output_path.display());
}
