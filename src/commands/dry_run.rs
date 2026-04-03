use std::path::Path;

use serde::Serialize;

use crate::output::Output;

const GB: u64 = 1024 * 1024 * 1024;

#[derive(Serialize)]
pub struct DryRunPlan {
    pub command: String,
    pub inputs: serde_json::Value,
    pub memory: MemoryEstimate,
    pub output_path: String,
}

#[derive(Serialize)]
pub struct MemoryEstimate {
    pub minimum: String,
    pub recommended: String,
    pub minimum_bytes: u64,
    pub recommended_bytes: u64,
}

impl MemoryEstimate {
    /// For DuckDB-only commands: any memory works (DuckDB spills to disk).
    pub fn duckdb_only() -> Self {
        Self {
            minimum: "4G".into(),
            recommended: "16G".into(),
            minimum_bytes: 4 * GB,
            recommended_bytes: 16 * GB,
        }
    }
}

pub fn emit(plan: &DryRunPlan, out: &dyn Output) {
    out.result_json(&serde_json::to_value(plan).unwrap_or_default());
}

pub fn file_size(path: &Path) -> u64 {
    std::fs::metadata(path).map(|m| m.len()).unwrap_or(0)
}

pub use super::data::human_size as human_bytes;
