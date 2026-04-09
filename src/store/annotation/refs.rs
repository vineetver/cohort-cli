//! Registry of attached external annotation databases.
//!
//! On-disk shape: `<store_root>/annotations/refs.toml`. Each entry names
//! an annotation source by alias and points at the directory holding the
//! data. The registry never owns or copies the bytes — it is only the
//! resolver from alias → typed reader.
//!
//! Default registry seeds from the user's `Config`: `favor-base` and
//! `favor-full` map to the configured FAVOR root + tier; `favor-tissue`
//! maps to the tissue add-on dir. `attach()` writes a refs.toml entry
//! that overrides or extends the defaults; the on-disk file wins when an
//! alias appears in both.

use std::collections::BTreeMap;
use std::path::PathBuf;

use serde::{Deserialize, Serialize};

use crate::config::{Config, Tier};
use crate::error::CohortError;
use crate::store::annotation::{AnnotationDb, TissueDb};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", tag = "kind")]
pub enum AnnotationKind {
    /// FAVOR base or full tier — chromosome-sharded `sorted.parquet`.
    Favor { tier: Tier },
    /// FAVOR tissue add-on — per-table parquet directories.
    Tissue,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotationRef {
    pub name: String,
    #[serde(flatten)]
    pub kind: AnnotationKind,
    pub path: PathBuf,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
struct RefsFile {
    #[serde(default)]
    refs: Vec<AnnotationRef>,
}

pub struct AnnotationRegistry {
    entries: BTreeMap<String, AnnotationRef>,
}

impl AnnotationRegistry {
    /// Build a registry by merging the on-disk `refs.toml` (if present)
    /// over the defaults derived from the user's `Config`. Entries in
    /// the file shadow defaults of the same alias.
    pub fn load(refs_path: PathBuf, config: &Config) -> Result<Self, CohortError> {
        let mut entries: BTreeMap<String, AnnotationRef> = BTreeMap::new();
        for entry in default_entries(config) {
            entries.insert(entry.name.clone(), entry);
        }
        if refs_path.exists() {
            let body = std::fs::read_to_string(&refs_path).map_err(|e| {
                CohortError::Resource(format!("read {}: {e}", refs_path.display()))
            })?;
            let parsed: RefsFile = toml::from_str(&body).map_err(|e| {
                CohortError::Input(format!("parse {}: {e}", refs_path.display()))
            })?;
            for entry in parsed.refs {
                entries.insert(entry.name.clone(), entry);
            }
        }
        Ok(Self { entries })
    }

    /// Open the FAVOR-style annotation DB at `name`. Errors if the
    /// alias is unknown or points at something that is not a Favor ref.
    pub fn open_db(&self, name: &str) -> Result<AnnotationDb, CohortError> {
        let r = self.require(name)?;
        match r.kind {
            AnnotationKind::Favor { tier } => {
                let parent = r.path.parent().ok_or_else(|| {
                    CohortError::Input(format!(
                        "annotation ref {name}: path {} has no parent",
                        r.path.display()
                    ))
                })?;
                AnnotationDb::open_tier(tier, parent)
            }
            AnnotationKind::Tissue => Err(CohortError::Input(format!(
                "annotation ref {name} is a tissue ref; use open_tissue"
            ))),
        }
    }

    /// Open the tissue add-on DB at `name`.
    pub fn open_tissue(&self, name: &str) -> Result<TissueDb, CohortError> {
        let r = self.require(name)?;
        match r.kind {
            AnnotationKind::Tissue => TissueDb::open(&r.path),
            AnnotationKind::Favor { .. } => Err(CohortError::Input(format!(
                "annotation ref {name} is a favor ref; use open_db"
            ))),
        }
    }

    fn require(&self, name: &str) -> Result<&AnnotationRef, CohortError> {
        self.entries.get(name).ok_or_else(|| {
            CohortError::Input(format!(
                "annotation ref {name} not registered; known: {}",
                self.entries
                    .keys()
                    .cloned()
                    .collect::<Vec<_>>()
                    .join(", ")
            ))
        })
    }
}

fn default_entries(config: &Config) -> Vec<AnnotationRef> {
    let mut out = Vec::new();
    let root = config.root_dir();
    if !root.as_os_str().is_empty() {
        out.push(AnnotationRef {
            name: "favor-base".into(),
            kind: AnnotationKind::Favor { tier: Tier::Base },
            path: root.join(Tier::Base.as_str()),
        });
        out.push(AnnotationRef {
            name: "favor-full".into(),
            kind: AnnotationKind::Favor { tier: Tier::Full },
            path: root.join(Tier::Full.as_str()),
        });
        out.push(AnnotationRef {
            name: "favor-tissue".into(),
            kind: AnnotationKind::Tissue,
            path: config.tissue_dir(),
        });
    }
    out
}

/// Convenience: pick the default alias for a tier.
pub fn favor_alias_for(tier: Tier) -> &'static str {
    match tier {
        Tier::Base => "favor-base",
        Tier::Full => "favor-full",
    }
}

