//! Engine-facing newtypes for the cohort layer.
//!
//! `VariantVcf` is the dense per-cohort variant index. The newtype is
//! `repr(transparent)` so `&[VariantVcf]` can be cast to `&[u32]` for the
//! sparse_g.bin reader without copying.

#[repr(transparent)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct VariantVcf(pub u32);

impl From<u32> for VariantVcf {
    #[inline]
    fn from(v: u32) -> Self {
        Self(v)
    }
}

/// `&[VariantVcf]` ↔ `&[u32]` zero-cost cast (both `repr(transparent)` over `u32`).
#[inline]
pub(crate) fn as_u32_slice(vcfs: &[VariantVcf]) -> &[u32] {
    // SAFETY: VariantVcf is `#[repr(transparent)]` over `u32`.
    unsafe { std::slice::from_raw_parts(vcfs.as_ptr() as *const u32, vcfs.len()) }
}

#[inline]
pub(crate) fn from_u32_slice(vcfs: &[u32]) -> &[VariantVcf] {
    // SAFETY: VariantVcf is `#[repr(transparent)]` over `u32`.
    unsafe { std::slice::from_raw_parts(vcfs.as_ptr() as *const VariantVcf, vcfs.len()) }
}

/// Sparse carrier output for a batch of variants. `entries[i]` aligns
/// with the i-th input vcf in the slice that produced this batch.
pub struct CarrierBatch {
    pub entries: Vec<super::variants::CarrierList>,
}
