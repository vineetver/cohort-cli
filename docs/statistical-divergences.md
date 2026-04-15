# Statistical Divergences from R

> [Back to README](../README.md) · [STAAR](staar.md) · [Validation](validation.md)

Where favor-cli differs from the R reference packages and why.

## SKAT p-value: moment-matching vs eigenvalue saddlepoint

R STAAR uses eigenvalue-based saddlepoint approximation for the mixture-of-chi-squared distribution. We use Liu et al. (2009) moment-matching (match first 4 cumulants to a scaled chi-squared). Validated within 2e-3 tolerance against R output.

Moment-matching is O(m) where m = number of eigenvalues. Saddlepoint is iterative per query. For genome-wide scans across thousands of genes the difference adds up.

Ref: `stats.rs:85-182`, `ground_truth_test.rs:12` (`TOL_SKAT = 2e-3`)

## ACAT-V per-variant test: chi-squared(1) vs t-distribution

R STAAR uses the t-distribution (df = n - p) for per-variant score tests on continuous traits. We use chi-squared(1) everywhere. For n > 100 the t and normal are indistinguishable; for n < 30 there could be a small difference in tail p-values.

Validated within 1e-4 tolerance against R output across all test genes.

Ref: `score.rs:260-274`, `ground_truth_test.rs:11` (`TOL = 1e-4`)

## Eigenvalue floor: 1e-8 fixed vs SKAT-R adaptive

R SKAT drops eigenvalues below `mean(positive_eigenvalues) / 100000`. We drop below 1e-8, matching the STAAR C++ code (`STAAR_O.cpp`, `MetaSTAAR_O_SMMAT.cpp`). The fixed floor is simpler and matches the package we're reimplementing. Could diverge from standalone SKAT-R results on genes with very small eigenvalue spread.

Ref: `stats.rs:100-103`

## CCT edge cases

Match R `CCT.R` exactly: small-p asymptotic `1/(p*pi)` below 1e-16, large-T asymptotic `1/(pi*T)` above 1e15. P-value floor at 5e-324 (smallest positive f64) to prevent underflow poisoning Cauchy combination downstream.

Ref: `stats.rs:18-83`
