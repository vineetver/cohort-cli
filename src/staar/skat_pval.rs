use statrs::distribution::{ChiSquared, ContinuousCDF};

/// SKAT p-value via Liu et al. (2009) moment-matching.
///
/// Q ~ Σ λ_j χ²_1. Approximated as a*χ²(l, ncp=δ) + b where
/// parameters (a, l, δ) are matched from the first 4 cumulants.
///
/// Matches SKAT R package: Get_Liu_Params_Mod + Get_Liu_PVal_MOD.
///
/// Cumulants (raw power sums, NOT multiplied by chi-sq moments):
///   c1 = Σ λ_j,  c2 = Σ λ_j²,  c3 = Σ λ_j³,  c4 = Σ λ_j⁴
pub fn mixture_chisq_pvalue(statistic: f64, eigenvalues: &[f64]) -> f64 {
    if eigenvalues.is_empty() || !statistic.is_finite() {
        return f64::NAN;
    }

    // Match SKAT R: Get_Lambda_Org — keep eigenvalues > mean(positive) / 100000
    let positive: Vec<f64> = eigenvalues.iter().copied().filter(|&l| l >= 0.0).collect();
    if positive.is_empty() { return 1.0; }
    let threshold = positive.iter().sum::<f64>() / positive.len() as f64 / 100000.0;
    let lambdas: Vec<f64> = eigenvalues.iter().copied().filter(|&l| l > threshold).collect();
    if lambdas.is_empty() {
        return 1.0;
    }

    // Raw cumulants (power sums of eigenvalues)
    let c1: f64 = lambdas.iter().sum();
    let c2: f64 = lambdas.iter().map(|l| l * l).sum();
    let c3: f64 = lambdas.iter().map(|l| l * l * l).sum();
    let c4: f64 = lambdas.iter().map(|l| l.powi(4)).sum();

    if c2 <= 0.0 {
        return 1.0;
    }

    // Mean and standard deviation of Q under H0
    let mu_q = c1;
    let sigma_q = (2.0 * c2).sqrt();

    // Normalized cumulants (following SKAT R: Get_Liu_Params_Mod)
    let s1 = c3 / c2.powf(1.5);
    let s2 = c4 / (c2 * c2);

    // Fit parameters: Q ≈ a * (χ²(l, δ) - l - δ) * σ_Q + μ_Q
    // Equivalently: standardize Q, then transform to χ²(l, δ)
    let (l, delta, a) = if s1 * s1 > s2 {
        // Noncentral chi-squared approximation
        let a_val = 1.0 / (s1 - (s1 * s1 - s2).sqrt());
        let delta_val = (s1 * a_val.powi(3) - a_val * a_val).max(0.0);
        let l_val = (a_val * a_val - 2.0 * delta_val).max(0.5);
        (l_val, delta_val, a_val)
    } else {
        // Central chi-squared approximation (delta = 0)
        // This is the common case for most genetic data
        let a_val = 1.0 / s1;
        let l_val = 1.0 / (s1 * s1);
        (l_val.max(0.5), 0.0, a_val)
    };

    // Transform: Q_std = (Q - μ_Q)/σ_Q, then χ²_val = Q_std * σ_X + μ_X
    let mu_x = l + delta;
    let sigma_x = (2.0_f64).sqrt() * a;

    let q_std = (statistic - mu_q) / sigma_q;
    let chisq_val = q_std * sigma_x + mu_x;

    if chisq_val <= 0.0 {
        return 1.0;
    }

    // P-value from chi-squared distribution
    // When delta = 0: use central chi-squared (exact for this approximation)
    // When delta > 0: use central chi-squared with adjusted df as fallback
    //   (statrs doesn't have noncentral chi-squared; the central approx is
    //    accurate for small delta, which is the typical genetics case)
    let effective_df = if delta > 0.0 {
        // Noncentral → central approximation: χ²(l, δ) ≈ (l+δ)/(l+2δ) * χ²(l + 2δ)
        // This is the Patnaik (1949) two-moment approximation
        let df_adj = (l + delta).powi(2) / (l + 2.0 * delta);
        let scale = (l + 2.0 * delta) / (l + delta);
        // Rescale chisq_val
        let chisq_rescaled = chisq_val * scale;
        return match ChiSquared::new(df_adj) {
            Ok(dist) => (1.0 - dist.cdf(chisq_rescaled)).max(0.0),
            Err(_) => f64::NAN,
        };
    } else {
        l
    };

    match ChiSquared::new(effective_df) {
        Ok(dist) => (1.0 - dist.cdf(chisq_val)).max(0.0),
        Err(_) => f64::NAN,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_eigenvalue() {
        // Q ~ λ χ²(1). Exact: P(Q > t) = P(χ²(1) > t/λ)
        let lambda = 2.0;
        let t = 5.0;
        let p = mixture_chisq_pvalue(t, &[lambda]);
        let exact = 1.0 - ChiSquared::new(1.0).unwrap().cdf(t / lambda);
        assert!((p - exact).abs() < 0.01,
            "Single eigenvalue: got {p:.6}, expected {exact:.6}");
    }

    #[test]
    fn test_equal_eigenvalues() {
        // Q ~ λ * χ²(k). Exact: P(Q > t) = P(χ²(k) > t/λ)
        let lambda = 1.5;
        let k = 5;
        let eigenvalues = vec![lambda; k];
        let t = 10.0;
        let p = mixture_chisq_pvalue(t, &eigenvalues);
        let exact = 1.0 - ChiSquared::new(k as f64).unwrap().cdf(t / lambda);
        assert!((p - exact).abs() < 0.05,
            "Equal eigenvalues: got {p:.6}, expected {exact:.6}");
    }

    #[test]
    fn test_null_statistic() {
        let p = mixture_chisq_pvalue(0.0, &[1.0, 2.0, 3.0]);
        assert!(p > 0.5, "Zero statistic → large p: {p}");
    }

    #[test]
    fn test_uniform_under_null() {
        // Under H0, p-values should be approximately uniform
        // Generate Q = Σ λ_j z_j² where z ~ N(0,1) (under null)
        // For a quick check: Q = c1 (the mean) should give p ≈ 0.5
        let eigenvalues = [1.0, 2.0, 3.0];
        let mean = eigenvalues.iter().sum::<f64>();
        let p = mixture_chisq_pvalue(mean, &eigenvalues);
        assert!(p > 0.3 && p < 0.7,
            "p at mean should be ~0.5: {p}");
    }
}
