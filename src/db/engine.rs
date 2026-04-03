use crate::error::FavorError;
use crate::resource::Resources;

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
