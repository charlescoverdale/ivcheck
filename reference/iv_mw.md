# Mourifie-Wan (2017) test for instrument validity

Reformulates the testable implications of Kitagawa (2015) as a set of
conditional moment inequalities and tests them in the intersection-
bounds framework of Chernozhukov, Lee, and Rosen (2013). Without
covariates `x`, `iv_mw` tests the same inequalities as `iv_kitagawa` and
reduces exactly to the variance-weighted Kitagawa test. With covariates,
`iv_mw` estimates the conditional CDFs `F(y, d | X = x, Z = z)`
nonparametrically via series regression, computes plug-in
heteroscedasticity-robust standard errors, and takes the sup over
`(y, x)` of the variance-weighted positive-part violation. Critical
values come from a multiplier bootstrap with adaptive moment selection
in the style of Andrews and Soares (2010).

## Usage

``` r
iv_mw(object, ...)

# Default S3 method
iv_mw(
  object,
  d,
  z,
  x = NULL,
  basis_order = 3L,
  x_grid_size = 20L,
  y_grid_size = 50L,
  adaptive = TRUE,
  grid = NULL,
  n_boot = 1000,
  alpha = 0.05,
  weights = NULL,
  parallel = TRUE,
  ...
)

# S3 method for class 'fixest'
iv_mw(
  object,
  x = NULL,
  basis_order = 3L,
  x_grid_size = 20L,
  y_grid_size = 50L,
  adaptive = TRUE,
  grid = NULL,
  n_boot = 1000,
  alpha = 0.05,
  weights = NULL,
  parallel = TRUE,
  ...
)

# S3 method for class 'ivreg'
iv_mw(
  object,
  x = NULL,
  basis_order = 3L,
  x_grid_size = 20L,
  y_grid_size = 50L,
  adaptive = TRUE,
  grid = NULL,
  n_boot = 1000,
  alpha = 0.05,
  weights = NULL,
  parallel = TRUE,
  ...
)
```

## Arguments

- object:

  For the default method: a numeric outcome vector. For the `fixest` and
  `ivreg` methods: a fitted instrumental variable model from
  [fixest::feols](https://lrberge.github.io/fixest/reference/feols.html)
  or
  [`ivreg::ivreg()`](https://zeileis.github.io/ivreg/reference/ivreg.html).

- ...:

  Further arguments passed to methods.

- d:

  Binary 0/1 treatment vector (default method only).

- z:

  Discrete instrument (numeric or factor, default method only).

- x:

  Optional numeric vector, matrix, or data frame of covariates. If
  supplied, the test is conditional on the first numeric column of `x`.
  If `NULL`, the test reduces to the unconditional Mourifie-Wan test.

- basis_order:

  Polynomial order of the series-regression basis used to estimate
  `F(y, d | X, Z)`. Default `3L` (cubic). Set to `"auto"` to select the
  basis order by 5-fold cross-validation over the candidates 2, 3, 4, 5
  with squared-error loss on the indicator regression. When `"auto"` is
  used, the bootstrap becomes post-selection-valid: the test statistic
  is compared to the maximum of the bootstrap statistics across the
  candidate orders, which controls size at the nominal level against any
  selection rule but is mildly conservative relative to a fixed-order
  test. Runtime with `"auto"` is approximately four times the
  fixed-order path.

- x_grid_size:

  Number of quantile points of `x` at which to evaluate the conditional
  CDFs. Default 20.

- y_grid_size:

  Number of quantile points of `y` at which to evaluate the
  inequalities. Default 50.

- adaptive:

  Logical. If `TRUE` (default), the bootstrap uses the adaptive moment
  selection of Andrews-Soares (2010) with tuning parameter
  `kappa_n = sqrt(log(log(n)))`. If `FALSE`, uses the plug-in
  least-favourable critical value (conservative).

- grid:

  Deprecated. Ignored; use `y_grid_size` and `x_grid_size` instead.

- n_boot:

  Number of multiplier-bootstrap replications. Default 1000.

- alpha:

  Significance level for the returned verdict. Default 0.05.

- weights:

  Optional survey weights. A non-negative numeric vector of length equal
  to the sample size. Scaled internally so the mean weight is 1.0
  (preserving effective sample-size interpretation). Applied to the
  empirical CDFs, the bootstrap multiplier process, and the
  variance-weighted standard errors.

- parallel:

  Logical. Run bootstrap replications in parallel on POSIX systems via
  [parallel::mclapply](https://rdrr.io/r/parallel/mclapply.html).
  Default `TRUE`.

## Value

An object of class `iv_test`; see
[iv_kitagawa](https://charlescoverdale.github.io/ivcheck/reference/iv_kitagawa.md)
for element descriptions. Additional elements:

- conditional:

  Logical, whether covariates were supplied.

- kappa_n:

  Andrews-Soares tuning parameter used (`NA` if not applicable).

## Details

The CLR framework targets conditional moment inequalities of the form
`E[m(W; theta) | X] <= 0` for all `X`. Applied to Kitagawa's (2015)
inequalities, the relevant moments are the positive-part differences of
the conditional joint CDFs `F(y, d | X, Z)` for each
`(d, z_low, z_high, y, x)` index. `iv_mw` estimates `F(y, d | X, Z)` by
series regression of the indicator `1{Y <= y, D = d}` on a polynomial
basis of `X` within each `Z` cell. Robust standard errors come from the
heteroscedasticity-consistent sandwich of the series regression.
Critical values are drawn by multiplier bootstrap: the bootstrap process
reuses the plug-in SE denominator and perturbs the residuals by
Rademacher weights, projected back through the basis. Adaptive moment
selection includes only moments whose observed studentised statistic is
within `kappa_n` of the inequality boundary, giving tighter critical
values when some inequalities are strictly slack.

## References

Mourifie, I. and Wan, Y. (2017). Testing Local Average Treatment Effect
Assumptions. *Review of Economics and Statistics*, 99(2), 305-313.
[doi:10.1162/REST_a_00622](https://doi.org/10.1162/REST_a_00622)

Chernozhukov, V., Lee, S., and Rosen, A. M. (2013). Intersection Bounds:
Estimation and Inference. *Econometrica*, 81(2), 667-737.
[doi:10.3982/ECTA8718](https://doi.org/10.3982/ECTA8718)

Imbens, G. W. and Angrist, J. D. (1994). Identification and Estimation
of Local Average Treatment Effects. *Econometrica*, 62(2), 467-475.
[doi:10.2307/2951620](https://doi.org/10.2307/2951620)

## See also

[`iv_kitagawa()`](https://charlescoverdale.github.io/ivcheck/reference/iv_kitagawa.md)
for the unconditional case,
[`iv_testjfe()`](https://charlescoverdale.github.io/ivcheck/reference/iv_testjfe.md)
for the judge-design test, and
[`iv_check()`](https://charlescoverdale.github.io/ivcheck/reference/iv_check.md)
for a one-shot wrapper that runs all applicable tests.

Other iv_tests:
[`iv_kitagawa()`](https://charlescoverdale.github.io/ivcheck/reference/iv_kitagawa.md),
[`iv_testjfe()`](https://charlescoverdale.github.io/ivcheck/reference/iv_testjfe.md)

## Examples

``` r
# \donttest{
set.seed(1)
n <- 500
z <- sample(0:1, n, replace = TRUE)
d <- rbinom(n, 1, 0.3 + 0.4 * z)
y <- rnorm(n, mean = d)
iv_mw(y, d, z, n_boot = 200, parallel = FALSE)
#> 
#> ── Mourifie-Wan (2017) ─────────────────────────────────────────────────────────
#> Sample size: 500
#> Statistic: "0.916", p-value: "1"
#> Verdict: cannot reject IV validity at 0.05
# }
```
