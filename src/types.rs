//! Canonical object model for FAVOR.
//!
//! One struct per concept. Downstream code composes these types — no parallel
//! structs, no field duplication, no `Vec<f64>` where `[f64; 11]` suffices.

use std::fmt;
use std::str::FromStr;

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Chromosome
// ---------------------------------------------------------------------------

/// Strongly-typed chromosome identifier.
///
/// `Ord` gives natural sort order: 1..22 < X < Y < MT.
/// Replaces all string-based chromosome fields and `chrom_sort_key()`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Chromosome {
    Autosome(u8), // 1–22
    X,
    Y,
    MT,
}

impl Chromosome {
    fn sort_key(self) -> (u8, u8) {
        match self {
            Chromosome::Autosome(n) => (0, n),
            Chromosome::X => (1, 0),
            Chromosome::Y => (1, 1),
            Chromosome::MT => (1, 2),
        }
    }
}

impl PartialOrd for Chromosome {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Chromosome {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.sort_key().cmp(&other.sort_key())
    }
}

impl fmt::Display for Chromosome {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Chromosome::Autosome(n) => write!(f, "chr{n}"),
            Chromosome::X => write!(f, "chrX"),
            Chromosome::Y => write!(f, "chrY"),
            Chromosome::MT => write!(f, "chrMT"),
        }
    }
}

impl FromStr for Chromosome {
    type Err = String;

    /// Parse chromosome from common representations:
    /// "chr1", "1", "chrX", "X", "chrM", "MT", "chrMT", etc.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let raw = s.strip_prefix("chr").unwrap_or(s);
        match raw {
            "X" | "x" => Ok(Chromosome::X),
            "Y" | "y" => Ok(Chromosome::Y),
            "M" | "MT" | "Mt" | "m" | "mt" => Ok(Chromosome::MT),
            _ => {
                let n: u8 = raw
                    .parse()
                    .map_err(|_| format!("invalid chromosome: {s}"))?;
                if (1..=22).contains(&n) {
                    Ok(Chromosome::Autosome(n))
                } else {
                    Err(format!("autosome out of range 1–22: {n}"))
                }
            }
        }
    }
}

impl Serialize for Chromosome {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for Chromosome {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let s = String::deserialize(deserializer)?;
        Chromosome::from_str(&s).map_err(serde::de::Error::custom)
    }
}

// ---------------------------------------------------------------------------
// AnnotationWeights
// ---------------------------------------------------------------------------

/// Fixed-size annotation weights for the 11 STAAR channels.
///
/// Replaces `Vec<f64>` annotation_weights, `WEIGHT_COL_START` / `WEIGHT_COL_COUNT`
/// magic constants, and `ANNOTATION_CHANNELS` (now `DISPLAY_NAMES`).
#[derive(Clone, Copy, Debug)]
pub struct AnnotationWeights(pub [f64; 11]);

impl Default for AnnotationWeights {
    fn default() -> Self {
        Self([0.0; 11])
    }
}

impl AnnotationWeights {
    /// Column names as written to parquet / read from DuckDB.
    pub const NAMES: [&str; 11] = [
        "w_cadd",
        "w_linsight",
        "w_fathmm_xf",
        "w_apc_epi_active",
        "w_apc_epi_repressed",
        "w_apc_epi_transcription",
        "w_apc_conservation",
        "w_apc_protein_function",
        "w_apc_local_nd",
        "w_apc_mutation_density",
        "w_apc_tf",
    ];

    /// Display names matching STAARpipeline R output.
    pub const DISPLAY_NAMES: [&str; 11] = [
        "cadd_phred",
        "linsight",
        "fathmm_xf",
        "apc_epigenetics_active",
        "apc_epigenetics_repressed",
        "apc_epigenetics_transcription",
        "apc_conservation",
        "apc_protein_function",
        "apc_local_nucleotide_diversity",
        "apc_mutation_density",
        "apc_transcription_factor",
    ];

    /// Apply learned priors from the feedback loop.
    ///
    /// Multiplies each channel weight by the corresponding prior.
    /// Priors default to `[1.0; 11]` (uniform — no adjustment).
    pub fn with_priors(&self, priors: &[f64; 11]) -> Self {
        let mut out = self.0;
        for i in 0..11 {
            out[i] *= priors[i];
        }
        Self(out)
    }

    /// Iterator over (column_name, value) pairs — for parquet / SQL generation.
    pub fn named_values(&self) -> impl Iterator<Item = (&'static str, f64)> + '_ {
        Self::NAMES.iter().zip(self.0.iter()).map(|(n, &v)| (*n, v))
    }
}

// ---------------------------------------------------------------------------
// RegulatoryFlags
// ---------------------------------------------------------------------------

/// Regulatory region flags used by mask predicates.
///
/// Replaces 4 loose `bool` fields repeated in the old AnnotatedVariant and
/// MetaVariant structs.
#[derive(Clone, Copy, Debug, Default)]
pub struct RegulatoryFlags {
    pub cage_promoter: bool,
    pub cage_enhancer: bool,
    pub ccre_promoter: bool,
    pub ccre_enhancer: bool,
}

// ---------------------------------------------------------------------------
// FunctionalAnnotation
// ---------------------------------------------------------------------------

/// Everything STAAR needs to classify (masks) and weight (score tests) a variant.
///
/// Mask predicates access `consequence` / `region_type`.
/// Score tests access `weights`.
#[derive(Clone, Debug)]
pub struct FunctionalAnnotation {
    pub region_type: String,
    pub consequence: String,
    pub cadd_phred: f64,
    pub revel: f64,
    pub regulatory: RegulatoryFlags,
    pub weights: AnnotationWeights,
}

// ---------------------------------------------------------------------------
// AnnotatedVariant (canonical)
// ---------------------------------------------------------------------------

/// Single canonical variant type used by masks, score tests, sumstats export,
/// result writing, and MetaSTAAR merge.
///
/// Replaces the parallel structs in `masks.rs` and `meta.rs`.
#[derive(Clone, Debug)]
pub struct AnnotatedVariant {
    pub chromosome: Chromosome,
    pub position: u32,
    pub ref_allele: String,
    pub alt_allele: String,
    pub maf: f64,
    pub gene_name: String,
    pub annotation: FunctionalAnnotation,
}

// ---------------------------------------------------------------------------
// MetaVariant
// ---------------------------------------------------------------------------

/// A variant in a meta-analysis context. Composes `AnnotatedVariant` instead
/// of duplicating its fields.
///
/// Where masks need `&AnnotatedVariant`, pass `&mv.variant`. Zero cost,
/// zero duplication. Deletes the need for `meta_to_annotated()`.
pub struct MetaVariant {
    pub variant: AnnotatedVariant,
    pub u_meta: f64,
    pub mac_total: i64,
    pub n_total: i64,
    pub study_segments: Vec<(usize, i32)>,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -- Chromosome --------------------------------------------------------

    #[test]
    fn parse_autosomes() {
        assert_eq!("chr1".parse::<Chromosome>().unwrap(), Chromosome::Autosome(1));
        assert_eq!("22".parse::<Chromosome>().unwrap(), Chromosome::Autosome(22));
        assert_eq!("chr10".parse::<Chromosome>().unwrap(), Chromosome::Autosome(10));
    }

    #[test]
    fn parse_sex_and_mt() {
        assert_eq!("chrX".parse::<Chromosome>().unwrap(), Chromosome::X);
        assert_eq!("Y".parse::<Chromosome>().unwrap(), Chromosome::Y);
        assert_eq!("chrM".parse::<Chromosome>().unwrap(), Chromosome::MT);
        assert_eq!("MT".parse::<Chromosome>().unwrap(), Chromosome::MT);
        assert_eq!("chrMT".parse::<Chromosome>().unwrap(), Chromosome::MT);
    }

    #[test]
    fn parse_invalid() {
        assert!("chr0".parse::<Chromosome>().is_err());
        assert!("chr23".parse::<Chromosome>().is_err());
        assert!("banana".parse::<Chromosome>().is_err());
    }

    #[test]
    fn display_roundtrip() {
        for s in ["chr1", "chr22", "chrX", "chrY", "chrMT"] {
            let c: Chromosome = s.parse().unwrap();
            assert_eq!(c.to_string(), s);
        }
    }

    #[test]
    fn ordering() {
        let mut chroms: Vec<Chromosome> = vec![
            Chromosome::MT,
            Chromosome::Autosome(22),
            Chromosome::X,
            Chromosome::Autosome(1),
            Chromosome::Y,
            Chromosome::Autosome(3),
        ];
        chroms.sort();
        assert_eq!(
            chroms,
            vec![
                Chromosome::Autosome(1),
                Chromosome::Autosome(3),
                Chromosome::Autosome(22),
                Chromosome::X,
                Chromosome::Y,
                Chromosome::MT,
            ]
        );
    }

    #[test]
    fn serde_roundtrip() {
        let c = Chromosome::Autosome(7);
        let json = serde_json::to_string(&c).unwrap();
        assert_eq!(json, "\"chr7\"");
        let back: Chromosome = serde_json::from_str(&json).unwrap();
        assert_eq!(back, c);
    }

    // -- AnnotationWeights -------------------------------------------------

    #[test]
    fn channel_counts() {
        assert_eq!(AnnotationWeights::NAMES.len(), 11);
        assert_eq!(AnnotationWeights::DISPLAY_NAMES.len(), 11);
    }

    #[test]
    fn with_priors_uniform() {
        let w = AnnotationWeights([1.0; 11]);
        let out = w.with_priors(&[1.0; 11]);
        assert_eq!(out.0, [1.0; 11]);
    }

    #[test]
    fn with_priors_scales() {
        let w = AnnotationWeights([2.0; 11]);
        let mut priors = [1.0; 11];
        priors[0] = 0.5;
        let out = w.with_priors(&priors);
        assert!((out.0[0] - 1.0).abs() < 1e-10);
        assert!((out.0[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn named_values_matches() {
        let w = AnnotationWeights([0.0; 11]);
        let pairs: Vec<_> = w.named_values().collect();
        assert_eq!(pairs.len(), 11);
        assert_eq!(pairs[0].0, "w_cadd");
        assert_eq!(pairs[10].0, "w_apc_tf");
    }

    // -- Display names match existing ANNOTATION_CHANNELS -------------------

    #[test]
    fn display_names_match_annotation_channels() {
        // These must stay in sync with weights.rs::ANNOTATION_CHANNELS
        let expected = [
            "cadd_phred",
            "linsight",
            "fathmm_xf",
            "apc_epigenetics_active",
            "apc_epigenetics_repressed",
            "apc_epigenetics_transcription",
            "apc_conservation",
            "apc_protein_function",
            "apc_local_nucleotide_diversity",
            "apc_mutation_density",
            "apc_transcription_factor",
        ];
        assert_eq!(AnnotationWeights::DISPLAY_NAMES, expected);
    }
}
