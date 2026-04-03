//! Cheap build and coordinate-base detection.
//!
//! Samples 100 rows from the input, probes annotation parquets to determine:
//! - hg38 vs hg19 (position match rate against known annotations)
//! - 0-based vs 1-based (position vs position+1 match rate)
//!
//! NEVER counts all rows. NEVER scans the full file. Reads exactly 100 rows.

use std::path::Path;

use super::{Analysis, BuildGuess, CoordBase};
use crate::annotation_db::AnnotationDb;
use crate::config::Config;
use crate::db::engine::DuckEngine;
use crate::error::FavorError;

/// Chromosomes to try probing, in priority order.
/// Falls back through the list until one has data in both input and annotations.
const PROBE_CHROMS: &[&str] = &["1", "2", "22", "10", "X"];

/// Run build and coordinate detection on an analysis.
/// Mutates analysis.build_guess and analysis.coord_base in place.
/// Requires DuckDB engine and configured annotation paths.
pub fn detect_build_and_coords(
    analysis: &mut Analysis,
    input_path: &Path,
    engine: &DuckEngine,
    config: &Config,
) -> Result<(), FavorError> {
    let chr_col = match &analysis.chr_col {
        Some(c) => c.clone(),
        None => return Ok(()),
    };
    let pos_col = match &analysis.pos_col {
        Some(c) => c.clone(),
        None => return Ok(()),
    };

    let ann_db = match AnnotationDb::open(config) {
        Ok(db) => db,
        Err(_) => return Ok(()),
    };

    let read_csv_expr = match analysis.delimiter {
        Some(delim) => format!(
            "read_csv('{}', header=true, all_varchar=true, delim={})",
            input_path.display(), delim.sql_literal()
        ),
        None => format!(
            "read_parquet('{}')",
            input_path.display()
        ),
    };

    // Try each chromosome until we get a sample with hits
    for probe_chrom in PROBE_CHROMS {
        let probe_parquet = match ann_db.chrom_parquet(probe_chrom) {
            Some(p) => p,
            None => continue,
        };

        let result = probe_chromosome(
            engine, &read_csv_expr, &chr_col, &pos_col,
            probe_chrom, &probe_parquet.to_string_lossy(),
        );

        match result {
            Ok(Some((rate_1based, rate_0based))) => {
                // Coordinate base
                if rate_1based > 0.5 {
                    analysis.coord_base = CoordBase::OneBased;
                } else if rate_0based > 0.5 {
                    analysis.coord_base = CoordBase::ZeroBased;
                }

                // Build
                let best_rate = rate_1based.max(rate_0based);
                if best_rate > 0.5 {
                    analysis.build_guess = BuildGuess::Hg38;
                } else if best_rate < 0.2 {
                    analysis.build_guess = BuildGuess::Hg19 {
                        match_rate_hg38: best_rate,
                        match_rate_hg19: 1.0 - best_rate,
                    };
                }
                return Ok(());
            }
            Ok(None) => continue, // no data for this chrom, try next
            Err(_) => continue,   // query failed, try next
        }
    }

    Ok(()) // couldn't detect from any chromosome
}

/// Probe a single chromosome. Returns Some((rate_1based, rate_0based)) or None if no data.
fn probe_chromosome(
    engine: &DuckEngine,
    read_csv_expr: &str,
    chr_col: &str,
    pos_col: &str,
    chrom: &str,
    annotation_path: &str,
) -> Result<Option<(f64, f64)>, FavorError> {
    // Sample 100 positions from input for this chromosome
    let sample_sql = format!(
        "CREATE OR REPLACE TEMP TABLE _ingest_probe AS
         SELECT CAST(\"{pos_col}\" AS INTEGER) AS pos
         FROM {read_csv_expr}
         WHERE UPPER(REGEXP_REPLACE(CAST(\"{chr_col}\" AS VARCHAR), '^chr', '', 'i')) = '{chrom}'
         LIMIT 100"
    );

    engine.execute(&sample_sql)?;

    // Check how many samples we got
    let count = query_single_i64(engine, "SELECT COUNT(*) FROM _ingest_probe")?;
    if count < 5 {
        let _ = engine.execute("DROP TABLE IF EXISTS _ingest_probe");
        return Ok(None); // not enough data
    }

    // 1-based match
    let hits_1 = query_single_i64(engine, &format!(
        "SELECT COUNT(*) FROM _ingest_probe s
         SEMI JOIN read_parquet('{annotation_path}') a ON s.pos = a.position"
    ))?;

    // 0-based match (pos + 1)
    let hits_0 = query_single_i64(engine, &format!(
        "SELECT COUNT(*) FROM _ingest_probe s
         SEMI JOIN read_parquet('{annotation_path}') a ON s.pos + 1 = a.position"
    ))?;

    let _ = engine.execute("DROP TABLE IF EXISTS _ingest_probe");

    let rate_1 = hits_1 as f64 / count as f64;
    let rate_0 = hits_0 as f64 / count as f64;

    Ok(Some((rate_1, rate_0)))
}

fn query_single_i64(engine: &DuckEngine, sql: &str) -> Result<i64, FavorError> {
    let conn = engine.connection();
    let mut stmt = conn.prepare(sql)
        .map_err(|e| FavorError::Analysis(format!("Probe query failed: {e}")))?;
    let mut rows = stmt.query([])
        .map_err(|e| FavorError::Analysis(format!("Probe query failed: {e}")))?;

    if let Some(row) = rows.next()
        .map_err(|e| FavorError::Analysis(format!("Probe read failed: {e}")))? {
        row.get(0).map_err(|e| FavorError::Analysis(format!("Probe parse failed: {e}")))
    } else {
        Ok(0)
    }
}
