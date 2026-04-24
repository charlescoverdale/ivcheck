# Frandsen-Lefgren-Leslie (2023) test for instrument validity in judge-fixed-effects designs

Jointly tests the local exclusion and monotonicity assumptions when the
instruments are a set of mutually exclusive dummy variables (the
leniency-of-assigned-judge design). Supports binary and multivalued
discrete treatments. Under the joint null, the per-judge mean outcome
`mu_j = E[Y | J = j]` must be a linear function of the per-judge
treatment propensities `P(D = d | J = j)`. Rejection is evidence that at
least one of exclusion or monotonicity fails.

## Usage

``` r
iv_testjfe(object, ...)

# Default S3 method
iv_testjfe(
  object,
  d,
  z,
  x = NULL,
  n_boot = 1000,
  alpha = 0.05,
  method = c("asymptotic", "bootstrap"),
  weights = NULL,
  basis_order = 1L,
  parallel = TRUE,
  ...
)

# S3 method for class 'fixest'
iv_testjfe(
  object,
  x = NULL,
  n_boot = 1000,
  alpha = 0.05,
  method = c("asymptotic", "bootstrap"),
  weights = NULL,
  basis_order = 1L,
  parallel = TRUE,
  ...
)

# S3 method for class 'ivreg'
iv_testjfe(
  object,
  x = NULL,
  n_boot = 1000,
  alpha = 0.05,
  method = c("asymptotic", "bootstrap"),
  weights = NULL,
  basis_order = 1L,
  parallel = TRUE,
  ...
)
```

## Arguments

- object:

  For the default method: a numeric outcome vector. For the `fixest` and
  `ivreg` methods: a fitted instrumental variable model from
  [fixest::feols](https://rdrr.io/pkg/fixest/man/feols.html) or
  [`ivreg::ivreg()`](https://rdrr.io/pkg/ivreg/man/ivreg.html).

- ...:

  Further arguments passed to methods.

- d:

  Binary 0/1 treatment vector (default method only).

- z:

  Factor, integer, or matrix of mutually exclusive dummy variables
  identifying the judge (or other random-assignment unit).

- x:

  Optional numeric vector, matrix, or data frame of covariates. If
  supplied, `y` and `d` are residualised on `x` before the per- judge
  means are computed.

- n_boot:

  Number of multiplier-bootstrap replications. Default 1000.

- alpha:

  Significance level for the returned verdict. Default 0.05.

- method:

  Reference distribution for the p-value. `"asymptotic"` (default) uses
  the chi-squared with `K - (basis_order + 1)` degrees of freedom.
  `"bootstrap"` uses the multiplier bootstrap of the restricted-model
  residual process. Asymptotic is fast and accurate for moderate `K`;
  bootstrap is preferred for small `K` or if errors are far from normal.

- weights:

  Optional survey weights. A non-negative numeric vector of length equal
  to the sample size. Scaled internally so the mean weight is 1.0
  (preserving effective sample-size interpretation). Applied to the
  empirical CDFs, the bootstrap multiplier process, and the
  variance-weighted standard errors.

- basis_order:

  Order of the polynomial basis used to approximate the outcome /
  propensity function `phi(p)` in Frandsen-Lefgren-Leslie (2023) step 1.
  Default `1L` reduces to the Sargan-Hansen overidentification form,
  which imposes constant treatment effects. Values above 1 relax this to
  `phi(p) = delta_0 + delta_1 p + delta_2 p^2 + ... + delta_m p^m` and
  test the joint-zero restriction on judge residuals under the richer
  fit. Only binary treatment is supported when `basis_order > 1`. The
  slope-bounded moment-inequality component of the FLL test is not
  implemented in v0.1.0 (deferred to v0.2.0).

- parallel:

  Logical. Run bootstrap replications in parallel on POSIX systems via
  [parallel::mclapply](https://rdrr.io/r/parallel/mclapply.html).
  Default `TRUE`.

## Value

An object of class `iv_test`; see
[iv_kitagawa](https://charlescoverdale.github.io/ivcheck/reference/iv_kitagawa.md)
for element descriptions. Additional elements:

- n_judges:

  Number of distinct judges / assignment groups.

- coef:

  Fitted weighted-LS slope and intercept of `mu_j` on `p_j`.

- pairwise_late:

  `K x K` matrix of pairwise Wald LATE estimates
  `(mu_j - mu_k) / (p_j - p_k)`. Under the null every entry estimates
  the common complier LATE.

- worst_pair:

  List identifying the judge pair with the largest deviation of its Wald
  LATE from the fitted slope; useful for diagnosing the source of a
  rejection.

## Details

Under the joint null, each pair of judges `(j, k)` identifies the same
complier LATE via the Wald estimator `(mu_j - mu_k) / (p_j - p_k)`. The
Frandsen-Lefgren-Leslie (2023) test is the overidentification test of
"all pairwise LATEs equal". Under binary treatment with WLS weighting,
that overidentification test is algebraically the weighted sum of
squared residuals from the linear fit `mu_j = alpha + beta * p_j`,
divided by a pooled variance estimator. `iv_testjfe` computes this
quadratic form and, by default, compares to a chi-squared distribution
with `K - 2` degrees of freedom (the FLL asymptotic form). The
multiplier bootstrap of the restricted residual process is available via
`method = "bootstrap"` for small-K robustness.

**Note on finite-sample size.** Per-judge propensities `p_j` enter the
test as estimated regressors. At modest per-judge sample sizes (`n_j`
below a few hundred), finite-sample binomial noise in `hat p_j`
compresses the distribution of the test statistic below the asymptotic
chi-squared reference, producing a test that is mildly conservative at
nominal 5 percent. Empirical size at `K = 20`, `N = 3000` is 1.5 percent
under the asymptotic method and 2.5 percent under the bootstrap. Both
methods sharpen toward nominal as `n_j` grows. The bootstrap is
recommended for publication-grade p-values at modest `n_j`.

The returned object includes `pairwise_late`, the `K x K` matrix of
pairwise Wald LATE estimates, and `worst_pair`, the judge pair with the
largest absolute deviation from the fitted slope. These are diagnostic
outputs in the sense of the paper's Figure 2: a pair whose Wald LATE
deviates far from the common slope is the first place to look when
investigating a rejection.

Multivalued treatment is supported: for `D` with `M + 1` distinct values
(`0, 1, ..., M`), the fit becomes a multiple WLS regression of `mu_j` on
the `M`-vector `(P(D = 1 | J), ..., P(D = M | J))` and the test
statistic is compared to `chi^2_{K - M - 1}` (FLL 2023 section 4).
`pairwise_late` and `worst_pair` are only defined for binary `D` and
return `NULL` otherwise.

## References

Frandsen, B. R., Lefgren, L. J., and Leslie, E. C. (2023). Judging Judge
Fixed Effects. *American Economic Review*, 113(1), 253-277.
[doi:10.1257/aer.20201860](https://doi.org/10.1257/aer.20201860)

Imbens, G. W. and Angrist, J. D. (1994). Identification and Estimation
of Local Average Treatment Effects. *Econometrica*, 62(2), 467-475.
[doi:10.2307/2951620](https://doi.org/10.2307/2951620)

## See also

[`iv_kitagawa()`](https://charlescoverdale.github.io/ivcheck/reference/iv_kitagawa.md)
for the unconditional binary-treatment test,
[`iv_mw()`](https://charlescoverdale.github.io/ivcheck/reference/iv_mw.md)
for the conditional version with covariates, and
[`iv_check()`](https://charlescoverdale.github.io/ivcheck/reference/iv_check.md)
for a one-shot wrapper that runs all applicable tests.

Other iv_tests:
[`iv_kitagawa()`](https://charlescoverdale.github.io/ivcheck/reference/iv_kitagawa.md),
[`iv_mw()`](https://charlescoverdale.github.io/ivcheck/reference/iv_mw.md)

## Examples

``` r
# \donttest{
set.seed(1)
n <- 2000
judge <- sample.int(20, n, replace = TRUE)
d <- rbinom(n, 1, 0.3 + 0.02 * judge)
y <- rnorm(n, mean = d)
iv_testjfe(y, d, judge, n_boot = 200, parallel = FALSE)
#> 
#> ── Frandsen-Lefgren-Leslie (2023) ──────────────────────────────────────────────
#> Sample size: 2000
#> Statistic: "27.2", p-value: "0.0751"
#> Verdict: cannot reject IV validity at 0.05
# }
```
