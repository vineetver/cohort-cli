use std::path::{Path, PathBuf};

use crate::config::{Config, Tier};
use crate::db::engine::DuckEngine;
use crate::db::query::query_strings;
use crate::error::FavorError;

/// The FAVOR annotation database on disk: hive-partitioned parquets at
/// `{root_dir}/{tier}/chromosome=*/sorted.parquet`.
pub struct AnnotationDb {
    root: PathBuf,
    tier: Tier,
}

impl AnnotationDb {
    pub fn open(config: &Config) -> Result<Self, FavorError> {
        Self::open_tier(config, config.data.tier)
    }

    pub fn open_tier(config: &Config, tier: Tier) -> Result<Self, FavorError> {
        let root = config.root_dir().join(tier.as_str());
        if !root.exists() {
            return Err(FavorError::DataMissing(format!(
                "Annotations not found at {}. Run `favor data pull` first.",
                root.display()
            )));
        }
        Ok(Self { root, tier })
    }

    pub fn chrom_parquet(&self, chrom: &str) -> Option<PathBuf> {
        let p = self.root.join(format!("chromosome={chrom}/sorted.parquet"));
        p.exists().then_some(p)
    }

    pub fn read_chrom(&self, chrom: &str) -> String {
        format!(
            "read_parquet('{}/chromosome={chrom}/sorted.parquet')",
            self.root.display()
        )
    }

    pub fn read_all(&self) -> String {
        format!(
            "read_parquet('{}/chromosome=*/sorted.parquet', hive_partitioning=true)",
            self.root.display()
        )
    }

    pub fn root(&self) -> &Path {
        &self.root
    }

    pub fn tier(&self) -> Tier {
        self.tier
    }

    /// Validate that the annotation parquet has the top-level columns STAAR requires.
    /// Catches base-vs-full tier mismatches before a multi-hour run.
    pub fn validate_staar_columns(&self, engine: &DuckEngine) -> Result<(), FavorError> {
        let required = ["gencode", "main", "cage", "apc"];
        let cols = query_strings(engine, &format!(
            "SELECT column_name FROM (DESCRIBE SELECT * FROM {} LIMIT 0)",
            self.read_all()
        ))?;
        for r in required {
            if !cols.iter().any(|c| c == r) {
                return Err(FavorError::DataMissing(format!(
                    "Annotation column '{}' missing. You may need favor-full tier. \
                     Current tier: {}. Run `favor data pull` to update.",
                    r, self.tier
                )));
            }
        }
        Ok(())
    }
}
