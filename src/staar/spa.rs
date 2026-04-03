use statrs::distribution::{Continuous, ContinuousCDF, Normal};

/// Saddlepoint approximation for binary trait score tests.
///
/// Provides calibrated p-values for imbalanced case-control designs where
/// the chi-squared approximation inflates Type I error. Falls back to the
/// normal approximation when p > 0.05 (root-finding is expensive).
///
/// Reference: Dey et al. (2017), "A Fast and Accurate Algorithm to Test
/// for Binary Phenotypes and Its Application to PheWAS"

const P_FILTER: f64 = 0.05;
const ROOT_TOL: f64 = 1e-9;
const ROOT_MAX_ITER: usize = 200;
const EXP_CLAMP: f64 = 500.0;

/// Two-sided SPA p-value for Burden or per-variant score tests.
///
/// `score`: the observed score statistic (w'G'(Y-μ) for Burden, G_j'(Y-μ) for single variant).
/// `mu`: fitted probabilities from logistic null model (length n).
/// `g`: (weighted) genotype vector (length n).
pub fn spa_pvalue(score: f64, mu: &[f64], g: &[f64]) -> f64 {
    debug_assert_eq!(mu.len(), g.len());
    if !score.is_finite() {
        return 1.0;
    }

    let norm = Normal::new(0.0, 1.0).unwrap();

    // Mean and variance of the score under H0
    let mean: f64 = mu.iter().zip(g).map(|(&m, &gi)| m * gi).sum();
    let var: f64 = mu
        .iter()
        .zip(g)
        .map(|(&m, &gi)| m * (1.0 - m) * gi * gi)
        .sum();
    if var <= 0.0 {
        return 1.0;
    }

    // Normal pre-filter: skip expensive root-finding for non-significant results
    let z = (score - mean) / var.sqrt();
    let p_normal = 2.0 * norm.cdf(-z.abs());
    if p_normal > P_FILTER {
        return p_normal;
    }

    // Two-sided via saddlepoint on both tails
    let p_upper = tail_prob(score, mu, g, mean, &norm);
    let p_lower = 1.0 - tail_prob(2.0 * mean - score, mu, g, mean, &norm);

    (2.0 * p_upper.min(p_lower)).clamp(0.0, 1.0)
}

/// Upper tail probability P(S >= q) via Lugannani-Rice formula.
fn tail_prob(q: f64, mu: &[f64], g: &[f64], mean: f64, norm: &Normal) -> f64 {
    let t_hat = match find_root(q, mu, g) {
        Some(t) if t.abs() > 1e-12 => t,
        _ => {
            let var: f64 = mu
                .iter()
                .zip(g)
                .map(|(&m, &gi)| m * (1.0 - m) * gi * gi)
                .sum();
            if var <= 0.0 {
                return if q > mean { 0.0 } else { 1.0 };
            }
            return norm.cdf(-(q - mean) / var.sqrt());
        }
    };

    let k = cgf(t_hat, mu, g);
    let k2 = cgf_d2(t_hat, mu, g);
    if k2 <= 0.0 {
        return 0.5;
    }

    let w_sq = 2.0 * (t_hat * q - k);
    if w_sq < 0.0 {
        return 0.5;
    }

    let w = t_hat.signum() * w_sq.sqrt();
    let v = t_hat * k2.sqrt();
    if w.abs() < 1e-12 {
        return 0.5;
    }

    (norm.cdf(-w) + norm.pdf(w) * (1.0 / w - 1.0 / v)).clamp(0.0, 1.0)
}

// ---------------------------------------------------------------------------
// Cumulant generating function and derivatives
//
// K(t)  = Σ log(1 - μ_i + μ_i exp(t g_i))
// K'(t) = Σ μ_i g_i exp(t g_i) / (1 - μ_i + μ_i exp(t g_i))
// K''(t)= Σ μ_i(1-μ_i) g_i² exp(t g_i) / (1 - μ_i + μ_i exp(t g_i))²
// ---------------------------------------------------------------------------

fn cgf(t: f64, mu: &[f64], g: &[f64]) -> f64 {
    mu.iter()
        .zip(g)
        .map(|(&m, &gi)| (1.0 - m + m * (t * gi).clamp(-EXP_CLAMP, EXP_CLAMP).exp()).ln())
        .sum()
}

fn cgf_d1(t: f64, mu: &[f64], g: &[f64]) -> f64 {
    mu.iter()
        .zip(g)
        .map(|(&m, &gi)| {
            let e = (t * gi).clamp(-EXP_CLAMP, EXP_CLAMP).exp();
            m * gi * e / (1.0 - m + m * e)
        })
        .sum()
}

fn cgf_d2(t: f64, mu: &[f64], g: &[f64]) -> f64 {
    mu.iter()
        .zip(g)
        .map(|(&m, &gi)| {
            let e = (t * gi).clamp(-EXP_CLAMP, EXP_CLAMP).exp();
            let d = 1.0 - m + m * e;
            m * (1.0 - m) * gi * gi * e / (d * d)
        })
        .sum()
}

/// Bisection root-finding for K'(t*) = q.
///
/// K'(t) is monotonically increasing so bisection always converges.
fn find_root(q: f64, mu: &[f64], g: &[f64]) -> Option<f64> {
    let k1_0 = cgf_d1(0.0, mu, g);
    let diff = q - k1_0;
    if diff.abs() < ROOT_TOL {
        return None;
    }

    let going_right = diff > 0.0;
    let (mut lo, mut hi) = if going_right { (0.0, 1.0) } else { (-1.0, 0.0) };

    for _ in 0..50 {
        if going_right {
            if cgf_d1(hi, mu, g) >= q {
                break;
            }
            hi *= 2.0;
        } else {
            if cgf_d1(lo, mu, g) <= q {
                break;
            }
            lo *= 2.0;
        }
    }

    for _ in 0..ROOT_MAX_ITER {
        let mid = 0.5 * (lo + hi);
        if (hi - lo).abs() < ROOT_TOL {
            return Some(mid);
        }
        if cgf_d1(mid, mu, g) < q {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cgf_at_zero_is_zero() {
        let mu = vec![0.5, 0.3, 0.8];
        let g = vec![1.0, 0.0, 2.0];
        assert!(cgf(0.0, &mu, &g).abs() < 1e-12);
    }

    #[test]
    fn k1_at_zero_equals_mean() {
        let mu = vec![0.5, 0.3, 0.8];
        let g = vec![1.0, 0.0, 2.0];
        let expected: f64 = mu.iter().zip(&g).map(|(m, gi)| m * gi).sum();
        assert!((cgf_d1(0.0, &mu, &g) - expected).abs() < 1e-12);
    }

    #[test]
    fn k2_at_zero_equals_variance() {
        let mu = vec![0.5, 0.3, 0.8];
        let g = vec![1.0, 0.0, 2.0];
        let expected: f64 = mu
            .iter()
            .zip(&g)
            .map(|(m, gi)| m * (1.0 - m) * gi * gi)
            .sum();
        assert!((cgf_d2(0.0, &mu, &g) - expected).abs() < 1e-12);
    }

    #[test]
    fn balanced_matches_normal() {
        let mu = vec![0.5; 100];
        let g: Vec<f64> = (0..100).map(|i| if i < 10 { 1.0 } else { 0.0 }).collect();
        let score = 3.0;
        let p_spa = spa_pvalue(score, &mu, &g);

        let norm = Normal::new(0.0, 1.0).unwrap();
        let mean: f64 = mu.iter().zip(&g).map(|(m, gi)| m * gi).sum();
        let var: f64 = mu
            .iter()
            .zip(&g)
            .map(|(m, gi)| m * (1.0 - m) * gi * gi)
            .sum();
        let z = (score - mean) / var.sqrt();
        let p_normal = 2.0 * norm.cdf(-z.abs());

        assert!(
            (p_spa - p_normal).abs() < 0.01,
            "Balanced: SPA={p_spa:.6}, normal={p_normal:.6}"
        );
    }

    #[test]
    fn score_at_mean_gives_large_p() {
        let mu = vec![0.1; 50];
        let g = vec![1.0; 50];
        let mean: f64 = mu.iter().zip(&g).map(|(m, gi)| m * gi).sum();
        let p = spa_pvalue(mean, &mu, &g);
        assert!(p > 0.9, "Score at mean should give large p: {p}");
    }

    #[test]
    fn extreme_imbalance_does_not_panic() {
        let mu = vec![0.001; 1000];
        let g: Vec<f64> = (0..1000).map(|i| if i < 5 { 1.0 } else { 0.0 }).collect();
        let score = 4.0;
        let p = spa_pvalue(score, &mu, &g);
        assert!(p >= 0.0 && p <= 1.0, "p must be valid: {p}");
    }

    #[test]
    fn zero_genotypes_gives_one() {
        let mu = vec![0.3; 50];
        let g = vec![0.0; 50];
        assert!((spa_pvalue(0.0, &mu, &g) - 1.0).abs() < 1e-10);
    }
}
