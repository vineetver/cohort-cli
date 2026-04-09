pub mod annotation;
pub mod backend;
pub mod cache;
pub mod cohort;
pub mod config;
pub mod ids;
pub mod layout;
pub mod list;
pub mod manifest;

use std::fs;

use crate::config::Config;
use crate::error::CohortError;

use annotation::AnnotationRegistry;
use backend::{Backend, LocalFs};
use cache::CacheStore;
use cohort::{CohortHandle, CohortId};
use config::StoreConfig;
use layout::Layout;

pub struct Store {
    layout: Layout,
    backend: Box<dyn Backend>,
}

impl Store {
    pub fn open(config: StoreConfig) -> Result<Self, CohortError> {
        let layout = Layout::new(config.root);
        fs::create_dir_all(layout.root()).map_err(|e| {
            CohortError::Resource(format!(
                "create store root {}: {e}",
                layout.root().display()
            ))
        })?;
        for sub in layout.known_subdirs() {
            fs::create_dir_all(&sub).map_err(|e| {
                CohortError::Resource(format!("create {}: {e}", sub.display()))
            })?;
        }
        Ok(Self {
            layout,
            backend: Box::new(LocalFs::new()),
        })
    }

    pub(crate) fn backend(&self) -> &dyn Backend {
        &*self.backend
    }

    /// Open a cohort handle by id. Lazy: does not read the manifest until
    /// a method on the returned handle is called.
    pub fn cohort(&self, id: &CohortId) -> CohortHandle<'_> {
        let dir = self.layout.cohort_dir(id);
        CohortHandle::new(self, id.clone(), dir)
    }

    /// Build the annotation registry from `<store_root>/annotations/refs.toml`
    /// merged over the defaults derived from `config`. Reads from disk
    /// every call — registries are cheap and the on-disk file may have
    /// changed between commands.
    pub fn annotations(&self, config: &Config) -> Result<AnnotationRegistry, CohortError> {
        let refs_path = self.layout.annotations_refs();
        if let Some(parent) = refs_path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| {
                CohortError::Resource(format!("create {}: {e}", parent.display()))
            })?;
        }
        AnnotationRegistry::load(refs_path, config)
    }

    /// Cache subsystem — score caches keyed by cohort id.
    pub fn cache(&self) -> CacheStore<'_> {
        CacheStore::new(&self.layout)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn open_creates_root_and_subdirs() {
        let tmp = tempfile::tempdir().unwrap();
        let root = tmp.path().join("store_root");
        let store = Store::open(StoreConfig { root: root.clone() }).unwrap();
        assert!(root.is_dir());
        assert!(store.layout.cohorts_root().is_dir());
        assert!(store.layout.lists_root().is_dir());
        assert!(store.layout.cache_root().is_dir());
    }

    #[test]
    fn backend_open_parquet_reads_testdata() {
        let path = std::path::Path::new(
            "testdata/ukb_chr22/ingested/chromosome=22/data.parquet",
        );
        if !path.exists() {
            eprintln!("skipping: {} not present", path.display());
            return;
        }
        let tmp = tempfile::tempdir().unwrap();
        let store = Store::open(StoreConfig { root: tmp.path().join("s") }).unwrap();
        let mut reader = store.backend().open_parquet(path).unwrap();
        let batch = reader
            .next()
            .expect("at least one batch")
            .expect("decode ok");
        assert!(batch.num_rows() > 0);
        assert!(batch.num_columns() > 0);
    }
}
