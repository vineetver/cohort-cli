//! Column alias map — THE single source of truth for GWAS column name normalization.

use super::{Ambiguity, ColumnMapping};

/// (lowercased_input_name, canonical_output_name)
///
/// When a raw column name (lowercased) matches the left side, it maps to the right side.
/// Ambiguities (e.g., A1 could be effect or other allele) are detected separately.
static ALIASES: &[(&str, &str)] = &[
    ("chromosome", "chromosome"),
    ("chr", "chromosome"),
    ("chrom", "chromosome"),
    ("#chrom", "chromosome"),
    ("#chr", "chromosome"),
    ("chr_name", "chromosome"),
    ("chrom_id", "chromosome"),
    ("hg38chrc", "chromosome"),

    ("position", "position"),
    ("pos", "position"),
    ("bp", "position"),
    ("base_pair_location", "position"),
    ("bp_hg38", "position"),
    ("bp_grch38", "position"),
    ("pos_hg38", "position"),
    ("start", "position"),
    ("chromstart", "position"),
    ("beg", "position"),

    ("ref", "ref"),
    ("reference", "ref"),
    ("ref_allele", "ref"),
    ("other_allele", "ref"),
    ("nea", "ref"),
    ("non_effect_allele", "ref"),

    ("alt", "alt"),
    ("alternate", "alt"),
    ("alt_allele", "alt"),
    ("effect_allele", "alt"),
    ("ea", "alt"),
    ("tested_allele", "alt"),
    ("risk_allele", "alt"),

    ("rsid", "rsid"),
    ("rsids", "rsid"),
    ("snp", "rsid"),
    ("snpid", "rsid"),
    ("variant_id", "rsid"),
    ("markername", "rsid"),
    ("rs", "rsid"),
    ("rs_id", "rsid"),

    ("beta", "beta"),
    ("effect", "beta"),
    ("effect_size", "beta"),
    ("b", "beta"),
    ("log_or", "beta"),

    ("se", "se"),
    ("stderr", "se"),
    ("standard_error", "se"),
    ("sebeta", "se"),

    ("pvalue", "pvalue"),
    ("p", "pvalue"),
    ("pval", "pvalue"),
    ("p_value", "pvalue"),
    ("p.value", "pvalue"),
    ("p_val", "pvalue"),

    ("log10p", "neglog10p"),
    ("neglog10p", "neglog10p"),
    ("mlogp", "neglog10p"),
    ("log10_p", "neglog10p"),
    ("nlog10p", "neglog10p"),

    ("z", "zscore"),
    ("zscore", "zscore"),
    ("z_score", "zscore"),
    ("z_stat", "zscore"),

    ("n", "n"),
    ("neff", "n"),
    ("n_eff", "n"),
    ("sample_size", "n"),
    ("n_total", "n"),

    ("pip", "pip"),
    ("posterior_prob", "pip"),
    ("posterior_inclusion_prob", "pip"),
    ("pp", "pip"),

    ("cs", "cs_id"),
    ("cs_id", "cs_id"),
    ("credible_set", "cs_id"),
    ("cs_index", "cs_id"),

    ("or", "odds_ratio"),
    ("odds_ratio", "odds_ratio"),

    ("af", "allele_freq"),
    ("freq", "allele_freq"),
    ("eaf", "allele_freq"),
    ("maf", "allele_freq"),
    ("allele_frequency", "allele_freq"),
    ("effect_allele_frequency", "allele_freq"),
];

/// Columns that are AMBIGUOUS — A1/A2 could be effect/other or ref/alt
/// depending on the study convention. We flag these instead of guessing.
static AMBIGUOUS_ALLELE_COLS: &[&str] = &["a1", "a2", "allele1", "allele2"];

/// Apply the alias map to raw input column names.
/// Returns (mapped, ambiguous, unmapped) — pure function, no side effects (Cmd III).
pub fn map_columns(raw_columns: &[String]) -> (Vec<ColumnMapping>, Vec<Ambiguity>, Vec<String>) {
    let mut mapped = Vec::new();
    let mut ambiguous = Vec::new();
    let mut unmapped = Vec::new();

    for col in raw_columns {
        let lower = col.to_lowercase().trim().to_string();

        // Check ambiguous allele columns first
        if AMBIGUOUS_ALLELE_COLS.contains(&lower.as_str()) {
            ambiguous.push(Ambiguity {
                column: col.clone(),
                candidates: vec!["ref", "alt"],
                reason: "A1/A2 convention varies across studies. \
                         Most common: A1=effect(alt), A2=other(ref). \
                         Check your study's README.",
            });
            continue;
        }

        // Look up in alias map
        if let Some((_, canonical)) = ALIASES.iter().find(|(alias, _)| *alias == lower) {
            mapped.push(ColumnMapping {
                input_name: col.clone(),
                canonical,
            });
        } else {
            unmapped.push(col.clone());
        }
    }

    (mapped, ambiguous, unmapped)
}

/// Get the DuckDB type to CAST a canonical column to.
pub fn canonical_type(canonical: &str) -> &'static str {
    match canonical {
        "chromosome" => "VARCHAR",
        "position" => "INTEGER",
        "ref" | "alt" | "rsid" => "VARCHAR",
        "beta" | "se" | "pvalue" | "zscore" | "neglog10p"
        | "pip" | "odds_ratio" | "allele_freq" => "DOUBLE",
        "n" | "cs_id" => "INTEGER",
        _ => "VARCHAR",
    }
}
