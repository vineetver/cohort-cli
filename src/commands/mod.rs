pub mod annotate;
pub mod enrich;
pub mod ingest;
pub mod inspect;
pub mod interpret;
pub mod meta_staar;
pub mod staar;

use std::path::{Path, PathBuf};

use serde::Serialize;

use crate::config::Tier;
use crate::error::CohortError;
use crate::output::Output;
use crate::staar::MaskCategory;

pub fn parse_mask_categories(masks: &[String]) -> Result<Vec<MaskCategory>, CohortError> {
    masks
        .iter()
        .map(|s| {
            s.parse::<MaskCategory>().map_err(|_| {
                CohortError::Input(format!(
                    "Unknown mask '{s}'. Available: coding, noncoding, sliding-window, scang, custom"
                ))
            })
        })
        .collect()
}

pub struct IngestConfig {
    pub inputs: Vec<PathBuf>,
    pub output: PathBuf,
    pub emit_sql: bool,
    pub build_override: Option<crate::cli::GenomeBuild>,
    pub annotations: Option<PathBuf>,
    pub cohort_id: Option<String>,
    pub rebuild: bool,
}

pub fn derive_cohort_id(genotypes: &Path) -> String {
    let stem = genotypes
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("cohort");
    let stem = stem
        .strip_suffix(".vcf.gz")
        .or_else(|| stem.strip_suffix(".vcf.bgz"))
        .or_else(|| stem.strip_suffix(".vcf"))
        .or_else(|| stem.strip_suffix(".bcf"))
        .unwrap_or(stem);
    let mut out = String::with_capacity(stem.len());
    let mut last_underscore = false;
    for ch in stem.chars() {
        if ch.is_ascii_alphanumeric() {
            out.push(ch);
            last_underscore = false;
        } else if !last_underscore && !out.is_empty() {
            out.push('_');
            last_underscore = true;
        }
    }
    while out.ends_with('_') {
        out.pop();
    }
    if out.is_empty() {
        "cohort".into()
    } else {
        out
    }
}

pub struct AnnotateConfig {
    pub input: PathBuf,
    pub output: PathBuf,
    pub tier: Tier,
}

pub struct EnrichConfig {
    pub input: PathBuf,
    pub output: PathBuf,
    pub tissue_name: String,
}

pub struct MetaStaarConfig {
    pub study_dirs: Vec<PathBuf>,
    pub mask_categories: Vec<MaskCategory>,
    pub maf_cutoff: f64,
    pub window_size: u32,
    pub output_dir: PathBuf,
}

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
    pub fn default_estimate() -> Self {
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

pub use crate::data::transfer::human_size as human_bytes;
