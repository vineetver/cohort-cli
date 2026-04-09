//! Derived caches keyed by cohort id. Today: per-phenotype score cache.

pub mod score_cache;

use std::path::PathBuf;

use crate::error::CohortError;
use crate::store::ids::{CacheKey, CohortId};
use crate::store::layout::Layout;

pub struct CacheStore<'a> {
    layout: &'a Layout,
}

impl<'a> CacheStore<'a> {
    pub(crate) fn new(layout: &'a Layout) -> Self {
        Self { layout }
    }

    /// Directory for one (cohort, cache_key) pair. The on-disk path is
    /// `<store_root>/cache/score_cache/<cohort>/<key>/`.
    pub fn score_cache_dir(&self, cohort: &CohortId, key: &CacheKey) -> PathBuf {
        self.layout.score_cache_dir(cohort, key)
    }

    /// Delete the score-cache directory owned by `cohort`. Called when
    /// `--rebuild-store` flips the content fingerprint and the existing
    /// cache no longer matches.
    pub fn prune_cohort(&self, cohort: &CohortId) -> Result<PruneSummary, CohortError> {
        let mut summary = PruneSummary::default();
        let score_dir = self
            .layout
            .cache_root()
            .join("score_cache")
            .join(cohort.as_str());
        if score_dir.is_dir() {
            summary.bytes_freed += dir_size(&score_dir);
            std::fs::remove_dir_all(&score_dir).map_err(|e| {
                CohortError::Resource(format!("rm {}: {e}", score_dir.display()))
            })?;
            summary.removed_score_caches = 1;
        }
        Ok(summary)
    }
}

#[derive(Debug, Default)]
pub struct PruneSummary {
    pub removed_score_caches: usize,
    pub bytes_freed: u64,
}

fn dir_size(path: &std::path::Path) -> u64 {
    let mut total = 0u64;
    if let Ok(entries) = std::fs::read_dir(path) {
        for entry in entries.flatten() {
            let p = entry.path();
            if p.is_file() {
                total += std::fs::metadata(&p).map(|m| m.len()).unwrap_or(0);
            } else if p.is_dir() {
                total += dir_size(&p);
            }
        }
    }
    total
}
