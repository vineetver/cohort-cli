pub mod ancestry;
pub mod cauchy;
pub mod geno_load;
pub mod genotype;
pub mod masks;
pub mod meta;
pub mod multi;
pub mod null_model;
pub mod scang;
pub mod score_test;
pub mod skat_pval;
pub mod spa;
pub mod summary;
pub mod sumstats;
pub mod weights;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TraitType {
    Continuous,
    Binary,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MaskCategory {
    Coding,
    Noncoding,
    SlidingWindow,
    Scang,
    Custom,
}

impl std::str::FromStr for MaskCategory {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "coding" => Ok(Self::Coding),
            "noncoding" => Ok(Self::Noncoding),
            "sliding-window" | "window" => Ok(Self::SlidingWindow),
            "scang" => Ok(Self::Scang),
            "custom" => Ok(Self::Custom),
            _ => Err(format!("unknown mask: {s}")),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum MaskType {
    PLof,
    Missense,
    DisruptiveMissense,
    PLofMissense,
    Synonymous,
    Ptv,
    PtvDs,
    Upstream,
    Downstream,
    Utr,
    PromoterCage,
    PromoterDhs,
    EnhancerCage,
    EnhancerDhs,
    Ncrna,
    SlidingWindow,
    Scang,
    Custom { name: String },
}

impl MaskType {
    pub fn file_stem(&self) -> String {
        match self {
            Self::PLof => "coding_pLoF".into(),
            Self::Missense => "coding_missense".into(),
            Self::DisruptiveMissense => "coding_disruptive_missense".into(),
            Self::PLofMissense => "coding_pLoF_missense".into(),
            Self::Synonymous => "coding_synonymous".into(),
            Self::Ptv => "coding_ptv".into(),
            Self::PtvDs => "coding_ptv_ds".into(),
            Self::Upstream => "noncoding_upstream".into(),
            Self::Downstream => "noncoding_downstream".into(),
            Self::Utr => "noncoding_utr".into(),
            Self::PromoterCage => "noncoding_promoter_CAGE".into(),
            Self::PromoterDhs => "noncoding_promoter_DHS".into(),
            Self::EnhancerCage => "noncoding_enhancer_CAGE".into(),
            Self::EnhancerDhs => "noncoding_enhancer_DHS".into(),
            Self::Ncrna => "noncoding_ncRNA".into(),
            Self::SlidingWindow => "sliding_window".into(),
            Self::Scang => "scang".into(),
            Self::Custom { name } => format!("custom_{name}"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct GeneResult {
    pub ensembl_id: String,
    pub gene_symbol: String,
    pub chromosome: String,
    pub start: u32,
    pub end: u32,
    pub n_variants: u32,
    pub cumulative_mac: u32,
    pub staar: score_test::StaarResult,
}

