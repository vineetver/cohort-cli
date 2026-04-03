use faer::Mat;

use crate::db::engine::DuckEngine;
use crate::error::FavorError;

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
