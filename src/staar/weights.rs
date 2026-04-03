use statrs::distribution::{Beta, Continuous};

/// Beta density weight: dbeta(maf, a1, a2).
///
/// For beta(1,25): upweights ultra-rare variants (MAF near 0).
/// For beta(1,1): uniform weight = 1.0 for all MAFs in (0,1).
///
/// Returns 0.0 for maf=0 or maf≥0.5 (monomorphic or not minor).
pub fn beta_density_weight(maf: f64, a1: f64, a2: f64) -> f64 {
    if maf <= 0.0 || maf >= 0.5 || !maf.is_finite() {
        return 0.0;
    }
    // Beta(1,1) is uniform — just return 1.0 without creating a distribution
    if (a1 - 1.0).abs() < 1e-10 && (a2 - 1.0).abs() < 1e-10 {
        return 1.0;
    }
    match Beta::new(a1, a2) {
        Ok(dist) => dist.pdf(maf),
        Err(_) => 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Convert a CADD PHRED score to an annotation weight in [0, 1].
    ///
    /// Matches the SQL in staar.rs: `1.0 - POW(10.0, -a.main.cadd.phred / 10.0)`.
    /// Non-finite or non-positive PHRED values yield 0.0 (no weight).
    fn phred_to_weight(phred: f64) -> f64 {
        if !phred.is_finite() || phred <= 0.0 {
            return 0.0;
        }
        1.0 - 10.0_f64.powf(-phred / 10.0)
    }

    #[test]
    fn test_beta_1_25_rare() {
        let w = beta_density_weight(0.001, 1.0, 25.0);
        assert!(w > 20.0, "Very rare variant should have high beta(1,25) weight: {w}");
    }

    #[test]
    fn test_beta_1_1_uniform() {
        let w = beta_density_weight(0.01, 1.0, 1.0);
        assert!((w - 1.0).abs() < 1e-10, "Beta(1,1) should give 1.0: {w}");
    }

    #[test]
    fn test_phred_conversion() {
        assert!((phred_to_weight(0.0)).abs() < 1e-10);
        assert!((phred_to_weight(10.0) - 0.9).abs() < 1e-10);
        assert!((phred_to_weight(20.0) - 0.99).abs() < 1e-10);
        assert_eq!(phred_to_weight(f64::NAN), 0.0);
    }
}
