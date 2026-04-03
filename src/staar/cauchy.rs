use std::f64::consts::PI;

/// Cauchy combination test (CCT) — equal weights.
///
/// T = (1/K) * sum(tan((0.5 - p_i) * pi))
/// p = 0.5 - arctan(T) / pi
///
/// Handles edge cases: p near 0 or 1 via clamping. Skips NaN.
pub fn cauchy_combine(p_values: &[f64]) -> f64 {
    cauchy_combine_weighted(p_values, &[])
}

/// Weighted Cauchy combination.
///
/// T = sum(w_i * tan((0.5 - p_i) * pi)) / sum(w_i)
///
/// If weights is empty, uses equal weights. Weights need not sum to 1.
/// Used by ACAT-V where annotation weights modulate each variant's contribution.
pub fn cauchy_combine_weighted(p_values: &[f64], weights: &[f64]) -> f64 {
    let use_weights = !weights.is_empty();
    if use_weights {
        assert_eq!(p_values.len(), weights.len());
    }

    let mut t_sum = 0.0;
    let mut w_sum = 0.0;
    let mut count = 0;

    for (i, &p) in p_values.iter().enumerate() {
        if !p.is_finite() || p <= 0.0 || p >= 1.0 {
            // Clamp valid p-values; skip truly invalid ones
            if p.is_nan() {
                continue;
            }
        }
        let p_clamped = p.clamp(1e-300, 1.0 - 1e-15);
        let w = if use_weights { weights[i].max(0.0) } else { 1.0 };
        if w == 0.0 {
            continue;
        }
        t_sum += w * ((0.5 - p_clamped) * PI).tan();
        w_sum += w;
        count += 1;
    }

    if count == 0 || w_sum == 0.0 {
        return f64::NAN;
    }

    let t = t_sum / w_sum;
    let p = 0.5 - t.atan() / PI;
    p.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cauchy_combine_basic() {
        let p = cauchy_combine(&[0.01, 0.05, 0.5]);
        assert!(p > 0.0 && p < 0.05, "Combined p should be small: {p}");
    }

    #[test]
    fn test_cauchy_combine_all_large() {
        let p = cauchy_combine(&[0.5, 0.5, 0.5]);
        assert!((p - 0.5).abs() < 0.01, "All 0.5 should give ~0.5: {p}");
    }

    #[test]
    fn test_cauchy_combine_extreme() {
        let p = cauchy_combine(&[1e-10, 0.99]);
        assert!(p < 0.01, "One very small p should dominate: {p}");
    }

    #[test]
    fn test_cauchy_combine_empty() {
        assert!(cauchy_combine(&[]).is_nan());
    }

    #[test]
    fn test_weighted_vs_equal() {
        let pvals = [0.01, 0.05, 0.5];
        let equal = cauchy_combine(&pvals);
        let weighted = cauchy_combine_weighted(&pvals, &[1.0, 1.0, 1.0]);
        assert!((equal - weighted).abs() < 1e-12);
    }

    #[test]
    fn test_weighted_emphasis() {
        let pvals = [0.01, 0.99];
        // High weight on the small p-value
        let p_high = cauchy_combine_weighted(&pvals, &[10.0, 1.0]);
        // High weight on the large p-value
        let p_low = cauchy_combine_weighted(&pvals, &[1.0, 10.0]);
        assert!(p_high < p_low, "Weighting small p should give smaller combined p");
    }
}
