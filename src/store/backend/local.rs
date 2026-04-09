//! Local filesystem backend.

use std::fs::File;
use std::path::Path;

use memmap2::MmapOptions;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

use crate::error::CohortError;

use super::{Backend, BoxedBatchReader, MappedBytes};

pub struct LocalFs;

impl LocalFs {
    pub fn new() -> Self {
        Self
    }
}

impl Backend for LocalFs {
    fn open_parquet(&self, path: &Path) -> Result<BoxedBatchReader, CohortError> {
        let file = File::open(path)
            .map_err(|e| CohortError::Resource(format!("open {}: {e}", path.display())))?;
        let reader = ParquetRecordBatchReaderBuilder::try_new(file)
            .map_err(|e| {
                CohortError::Resource(format!("parquet header {}: {e}", path.display()))
            })?
            .build()
            .map_err(|e| {
                CohortError::Resource(format!("parquet reader {}: {e}", path.display()))
            })?;
        Ok(Box::new(reader))
    }

    fn mmap(&self, path: &Path) -> Result<MappedBytes, CohortError> {
        let file = File::open(path)
            .map_err(|e| CohortError::Resource(format!("open {}: {e}", path.display())))?;
        let mmap = unsafe { MmapOptions::new().map(&file) }
            .map_err(|e| CohortError::Resource(format!("mmap {}: {e}", path.display())))?;
        Ok(MappedBytes::new(mmap))
    }
}
