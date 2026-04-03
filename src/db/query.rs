use crate::db::engine::DuckEngine;
use crate::error::FavorError;

/// Execute SQL returning a single i64 scalar (COUNT, MAX, etc). Returns 0 if no rows.
pub fn query_scalar(engine: &DuckEngine, sql: &str) -> Result<i64, FavorError> {
    let conn = engine.connection();
    let mut stmt = conn.prepare(sql)
        .map_err(|e| FavorError::Analysis(format!("{e}")))?;
    let mut rows = stmt.query([])
        .map_err(|e| FavorError::Analysis(format!("{e}")))?;
    match rows.next() {
        Ok(Some(row)) => row.get(0).map_err(|e| FavorError::Analysis(format!("{e}"))),
        _ => Ok(0),
    }
}

/// Execute SQL returning a Vec<String> from the first column.
pub fn query_strings(engine: &DuckEngine, sql: &str) -> Result<Vec<String>, FavorError> {
    let conn = engine.connection();
    let mut stmt = conn.prepare(sql)
        .map_err(|e| FavorError::Analysis(format!("{e}")))?;
    let mut rows = stmt.query([])
        .map_err(|e| FavorError::Analysis(format!("{e}")))?;
    let mut result = Vec::new();
    while let Some(row) = rows.next()
        .map_err(|e| FavorError::Analysis(format!("{e}")))? {
        let v: String = row.get(0).map_err(|e| FavorError::Analysis(format!("{e}")))?;
        result.push(v);
    }
    Ok(result)
}
