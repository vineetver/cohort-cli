//! DuckDB wrapper: engine + query helpers.

use crate::error::FavorError;
use crate::resource::Resources;

// ---------------------------------------------------------------------------
// DuckEngine (was db/engine.rs)
// ---------------------------------------------------------------------------

pub struct DuckEngine {
    conn: duckdb::Connection,
}

impl DuckEngine {
    pub fn new(resources: &Resources) -> Result<Self, FavorError> {
        let conn = duckdb::Connection::open_in_memory()
            .map_err(|e| FavorError::Resource(format!("DuckDB init failed: {e}")))?;

        conn.execute_batch(&format!(
            "SET memory_limit = '{}';
             SET threads TO {};
             SET temp_directory = '{}';
             SET preserve_insertion_order = true;",
            resources.duckdb_memory(),
            resources.threads,
            resources.temp_dir.display(),
        ))
        .map_err(|e| FavorError::Resource(format!("DuckDB config failed: {e}")))?;

        Ok(Self { conn })
    }

    pub fn execute(&self, sql: &str) -> Result<(), FavorError> {
        self.conn
            .execute_batch(sql)
            .map_err(|e| FavorError::Analysis(format!("Query failed: {e}")))?;
        Ok(())
    }

    pub fn connection(&self) -> &duckdb::Connection {
        &self.conn
    }
}

// ---------------------------------------------------------------------------
// Query helpers (was db/query.rs)
// ---------------------------------------------------------------------------

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
