//! Closed-loop learning: accumulate signal across STAAR runs.
//!
//! Phase 1: type definitions and load/save stubs only.
//! The feedback loop implementation is Phase 9.

use std::path::Path;

use serde::{Deserialize, Serialize};

const HISTORY_FILE: &str = ".favor_history.json";

/// What happened in one pipeline run.
#[derive(Serialize, Deserialize)]
pub struct RunRecord {
    pub run_id: String,
    pub timestamp: String,
    pub trait_name: String,
    pub n_samples: usize,
    pub maf_cutoff: f64,
    pub masks: Vec<String>,

    /// Per-channel: fraction of significant genes where this channel's
    /// weighted test produced the best (lowest) p-value among the 11.
    pub channel_lift: [f64; 11],

    /// Per-mask: fraction of tested genes that passed significance.
    pub mask_yield: Vec<(String, f64)>,

    /// MAC distribution of significant hits.
    pub signal_mac_bins: Vec<(u32, u32, usize)>,

    pub n_significant: usize,
    pub top_genes: Vec<String>,
}

/// Accumulated intelligence across runs.
/// Persisted to `.favor_history.json` in the output directory.
#[derive(Serialize, Deserialize, Default)]
pub struct PipelineHistory {
    pub version: u32,
    pub runs: Vec<RunRecord>,
}

impl PipelineHistory {
    /// Load history from the output directory. Returns empty history if the
    /// file doesn't exist or can't be parsed.
    pub fn load(output_dir: &Path) -> Self {
        let path = output_dir.join(HISTORY_FILE);
        std::fs::read_to_string(&path)
            .ok()
            .and_then(|s| serde_json::from_str(&s).ok())
            .unwrap_or_default()
    }

    /// Persist history to the output directory.
    pub fn save(&self, output_dir: &Path) -> std::io::Result<()> {
        let path = output_dir.join(HISTORY_FILE);
        let json = serde_json::to_string_pretty(self)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        std::fs::write(&path, json)
    }

    /// Append a run record.
    pub fn append(&mut self, record: RunRecord) {
        self.runs.push(record);
    }

    /// Learned channel priors from accumulated signal.
    /// Returns `[1.0; 11]` (uniform) until Phase 9 implements EMA.
    pub fn channel_priors(&self) -> [f64; 11] {
        [1.0; 11]
    }

    /// Masks that have produced >0 yield in any prior run.
    pub fn productive_masks(&self) -> Vec<String> {
        Vec::new()
    }

    /// Suggest a tighter MAF cutoff if signal concentrates at ultra-rare MAC.
    pub fn recommended_maf(&self) -> Option<f64> {
        None
    }

    /// Genes significant in 2+ runs.
    pub fn recurrent_genes(&self) -> Vec<String> {
        Vec::new()
    }
}
