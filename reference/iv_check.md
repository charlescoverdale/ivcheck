# Run all applicable IV-validity tests on a fitted model

Detects which tests are applicable from the structure of the fitted
instrumental variable model and runs them. Returns a tidy summary with a
one-line verdict.

## Usage

``` r
iv_check(model, tests = "all", alpha = 0.05, n_boot = 1000, ...)
```

## Arguments

- model:

  A fitted IV model from
  [fixest::feols](https://lrberge.github.io/fixest/reference/feols.html)
  or
  [`ivreg::ivreg()`](https://zeileis.github.io/ivreg/reference/ivreg.html).

- tests:

  Character vector of test names to run, or `"all"` (the default) to run
  every applicable test.

- alpha:

  Significance level for the verdict. Default 0.05.

- n_boot:

  Number of bootstrap replications. Default 1000.

- ...:

  Further arguments passed to each underlying test.

## Value

An object of class `iv_check` containing a data frame with one row per
test (test name, statistic, p-value, verdict) plus an overall verdict
string.

## Details

Applicability is determined by:

- Kitagawa (2015) applies to any binary treatment with a discrete
  instrument.

- Mourifie-Wan (2017) applies to the same case, and additionally
  supports covariates.

- Frandsen-Lefgren-Leslie (2023) applies when the instrument is a set of
  mutually exclusive dummy variables (judge-IV / group design).

## Examples

``` r
# \donttest{
if (requireNamespace("fixest", quietly = TRUE)) {
  set.seed(1)
  n <- 500
  df <- data.frame(
    z = sample(0:1, n, replace = TRUE),
    x = rnorm(n)
  )
  df$d <- rbinom(n, 1, 0.3 + 0.4 * df$z)
  df$y <- rnorm(n, mean = df$d + 0.5 * df$x)
  m <- fixest::feols(y ~ x | d ~ z, data = df)
  iv_check(m, n_boot = 200)
}
#> 
#> ── IV validity diagnostic ──────────────────────────────────────────────────────
#> Kitagawa (2015): stat = "0.845", p = "1", pass
#> Mourifie-Wan (2017): stat = "42.4", p = "0.455", pass
#> Overall: cannot reject IV validity at 0.05.
# }
```
