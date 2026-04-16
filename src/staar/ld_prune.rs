//! LD pruning via sequential conditional analysis.
//!
//! Forward selection on conditional p-values, not r² correlations. Starts
//! with the most significant marginal variant, then at each step refits
//! the null on `[X, G_known]` and picks the candidate with the smallest
//! conditional p-value. Stops when no remaining candidate beats
//! `cond_p_thresh`.
//!
//! Mirrors STAARpipeline R/LD_pruning.R lines 146–185 on the gaussian,
//! unrelated, single-trait path.

use faer::Mat;

use crate::error::CohortError;
use crate::staar::carrier::sparse_score::{carriers_to_dense_compact, individual_score_test};
use crate::staar::carrier::AnalysisVectors;
use crate::staar::model::{augment_covariates, fit_glm};
use crate::store::cohort::types::VariantVcf;
use crate::store::cohort::variants::CarrierList;
use crate::store::cohort::ChromosomeView;
use crate::types::Chromosome;

#[derive(Clone, Debug)]
pub struct Candidate {
    pub position: u32,
    pub ref_allele: String,
    pub alt_allele: String,
}

#[derive(Clone, Debug)]
pub struct KeptVariant {
    pub chromosome: Chromosome,
    pub position: u32,
    pub ref_allele: String,
    pub alt_allele: String,
    /// −log10 of the p-value under which this variant entered the pruned
    /// set. The first pick carries its marginal p; later picks carry their
    /// conditional p at the time of selection.
    pub entry_log10p: f64,
}

pub struct LdPruneParams {
    pub maf_cutoff: f64,
    pub cond_p_thresh: f64,
}

/// Resolve candidates to cohort-indexed carriers, MAF-filter, and hand off
/// to `ld_prune_from_carriers`. Keeps disk IO and math isolated so tests
/// can exercise the loop with synthetic carrier lists.
pub fn ld_prune_chromosome(
    view: &ChromosomeView<'_>,
    chromosome: Chromosome,
    y: &Mat<f64>,
    x_base: &Mat<f64>,
    pheno_mask: &[bool],
    candidates: &[Candidate],
    params: &LdPruneParams,
) -> Result<Vec<KeptVariant>, CohortError> {
    if candidates.is_empty() {
        return Ok(Vec::new());
    }

    let index = view.index()?;
    let all_entries = index.all_entries();

    let mut matched: Vec<(u32, Candidate)> = Vec::with_capacity(candidates.len());
    for c in candidates {
        for (i, e) in all_entries.iter().enumerate() {
            if e.position == c.position
                && e.ref_allele.as_ref() == c.ref_allele.as_str()
                && e.alt_allele.as_ref() == c.alt_allele.as_str()
                && e.maf > params.maf_cutoff
            {
                matched.push((i as u32, c.clone()));
                break;
            }
        }
    }
    if matched.is_empty() {
        return Ok(Vec::new());
    }

    matched.sort_by_key(|(v, _)| *v);
    let sorted_vcfs: Vec<VariantVcf> = matched.iter().map(|(v, _)| VariantVcf(*v)).collect();
    let sorted_candidates: Vec<Candidate> = matched.into_iter().map(|(_, c)| c).collect();
    let carriers: Vec<CarrierList> = view.carriers_batch(&sorted_vcfs)?.entries;

    ld_prune_from_carriers(
        chromosome,
        &sorted_candidates,
        &carriers,
        y,
        x_base,
        pheno_mask,
        params.cond_p_thresh,
    )
}

/// Core forward-selection loop over pre-loaded carriers. `candidates[i]`
/// identifies the variant that produced `carriers[i]`; vectors must have
/// the same length.
pub fn ld_prune_from_carriers(
    chromosome: Chromosome,
    candidates: &[Candidate],
    carriers: &[CarrierList],
    y: &Mat<f64>,
    x_base: &Mat<f64>,
    pheno_mask: &[bool],
    cond_p_thresh: f64,
) -> Result<Vec<KeptVariant>, CohortError> {
    assert_eq!(candidates.len(), carriers.len());
    let m = carriers.len();

    if m == 0 {
        return Ok(Vec::new());
    }

    let null0 = fit_glm(y, x_base);
    let analysis0 = AnalysisVectors::from_null_model(&null0, pheno_mask)?;

    if m == 1 {
        let s = individual_score_test(&carriers[0], &analysis0);
        let c = &candidates[0];
        return Ok(vec![KeptVariant {
            chromosome,
            position: c.position,
            ref_allele: c.ref_allele.clone(),
            alt_allele: c.alt_allele.clone(),
            entry_log10p: safe_log10p(s.pvalue),
        }]);
    }

    let mut best = 0usize;
    let mut best_log10p = f64::NEG_INFINITY;
    for (i, carrier) in carriers.iter().enumerate() {
        let s = individual_score_test(carrier, &analysis0);
        let lp = safe_log10p(s.pvalue);
        if lp > best_log10p {
            best_log10p = lp;
            best = i;
        }
    }

    let mut known: Vec<usize> = vec![best];
    let mut kept: Vec<KeptVariant> = Vec::with_capacity(4);
    kept.push(kept_variant(chromosome, &candidates[best], best_log10p));

    let cond_log10_thresh = -cond_p_thresh.log10();

    loop {
        if known.len() == m {
            break;
        }

        let known_carriers: Vec<CarrierList> =
            known.iter().map(|&i| carriers[i].clone()).collect();
        let g_known = carriers_to_dense_compact(
            &known_carriers,
            &analysis0.vcf_to_pheno,
            analysis0.n_pheno,
        );
        let x_cond = augment_covariates(x_base, &g_known);
        let null_cond = fit_glm(y, &x_cond);
        let analysis_cond = AnalysisVectors::from_null_model(&null_cond, pheno_mask)?;

        let mut pick: Option<usize> = None;
        let mut pick_log10p = cond_log10_thresh;
        for (i, carrier) in carriers.iter().enumerate() {
            if known.contains(&i) {
                continue;
            }
            let s = individual_score_test(carrier, &analysis_cond);
            let lp = safe_log10p(s.pvalue);
            if lp > pick_log10p {
                pick_log10p = lp;
                pick = Some(i);
            }
        }

        let Some(idx) = pick else { break };
        known.push(idx);
        kept.push(kept_variant(chromosome, &candidates[idx], pick_log10p));
    }

    Ok(kept)
}

fn kept_variant(chromosome: Chromosome, c: &Candidate, entry_log10p: f64) -> KeptVariant {
    KeptVariant {
        chromosome,
        position: c.position,
        ref_allele: c.ref_allele.clone(),
        alt_allele: c.alt_allele.clone(),
        entry_log10p,
    }
}

#[inline]
fn safe_log10p(p: f64) -> f64 {
    if p > 0.0 {
        -p.log10()
    } else {
        f64::INFINITY
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::store::cohort::variants::{CarrierEntry, CarrierList};

    fn dense_carrier(n: usize, dosages: impl Fn(usize) -> u8) -> CarrierList {
        let entries: Vec<CarrierEntry> = (0..n)
            .map(|i| CarrierEntry {
                sample_idx: i as u32,
                dosage: dosages(i),
            })
            .filter(|e| e.dosage != 0)
            .collect();
        CarrierList { entries }
    }

    fn fake_candidate(pos: u32) -> Candidate {
        Candidate {
            position: pos,
            ref_allele: "A".into(),
            alt_allele: "T".into(),
        }
    }

    #[test]
    fn single_variant_returns_one_row() {
        let n = 200;
        let pheno_mask = vec![true; n];
        let y = Mat::from_fn(n, 1, |i, _| (i as f64) * 0.01);
        let x = Mat::from_fn(n, 1, |_, _| 1.0);
        let carriers = vec![dense_carrier(n, |i| if i < 10 { 1 } else { 0 })];
        let cands = vec![fake_candidate(100)];

        let kept = ld_prune_from_carriers(
            Chromosome::Autosome(1),
            &cands,
            &carriers,
            &y,
            &x,
            &pheno_mask,
            1e-4,
        )
        .unwrap();
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].position, 100);
    }

    #[test]
    fn perfectly_collinear_variant_is_pruned() {
        // Two variants carry identical genotype; the partner should drop
        // out on the conditional pass because its signal is already
        // absorbed by the first pick.
        let n = 400;
        let pheno_mask = vec![true; n];
        // y correlates with carrier pattern so marginal p is tiny.
        let y = Mat::from_fn(n, 1, |i, _| if i < 40 { 2.0 } else { 0.0 });
        let x = Mat::from_fn(n, 1, |_, _| 1.0);
        let pattern = |i: usize| if i < 40 { 1u8 } else { 0u8 };
        let carriers = vec![
            dense_carrier(n, pattern),
            dense_carrier(n, pattern),
        ];
        let cands = vec![fake_candidate(100), fake_candidate(200)];

        let kept = ld_prune_from_carriers(
            Chromosome::Autosome(1),
            &cands,
            &carriers,
            &y,
            &x,
            &pheno_mask,
            1e-4,
        )
        .unwrap();
        assert_eq!(kept.len(), 1, "collinear partner must not be kept");
    }

    #[test]
    fn independent_signals_both_kept() {
        // Two orthogonal carrier patterns on two disjoint halves of y.
        let n = 400;
        let pheno_mask = vec![true; n];
        let mut y_vals = vec![0.0_f64; n];
        for v in y_vals.iter_mut().take(20) {
            *v = 3.0;
        }
        for v in y_vals.iter_mut().take(220).skip(200) {
            *v = 3.0;
        }
        let y = Mat::from_fn(n, 1, |i, _| y_vals[i]);
        let x = Mat::from_fn(n, 1, |_, _| 1.0);
        let carriers = vec![
            dense_carrier(n, |i| if i < 20 { 1 } else { 0 }),
            dense_carrier(n, |i| if (200..220).contains(&i) { 1 } else { 0 }),
        ];
        let cands = vec![fake_candidate(100), fake_candidate(500)];

        let kept = ld_prune_from_carriers(
            Chromosome::Autosome(1),
            &cands,
            &carriers,
            &y,
            &x,
            &pheno_mask,
            1e-4,
        )
        .unwrap();
        assert_eq!(kept.len(), 2, "orthogonal signals must both survive");
    }

    #[test]
    fn empty_input_returns_empty() {
        let n = 10;
        let y = Mat::<f64>::zeros(n, 1);
        let x = Mat::from_fn(n, 1, |_, _| 1.0);
        let kept = ld_prune_from_carriers(
            Chromosome::Autosome(1),
            &[],
            &[],
            &y,
            &x,
            &vec![true; n],
            1e-4,
        )
        .unwrap();
        assert!(kept.is_empty());
    }
}
