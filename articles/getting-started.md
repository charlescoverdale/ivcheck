# Getting started with ivcheck

## What this package does

`ivcheck` tests the identifying assumptions of instrumental variable
(IV) models: the local exclusion restriction and monotonicity conditions
required for LATE identification (Imbens and Angrist 1994). The package
implements three tests from the methodological literature that share
these assumptions as their null:

- [`iv_kitagawa()`](https://charlescoverdale.github.io/ivcheck/reference/iv_kitagawa.md):
  Kitagawa (2015), binary treatment with a discrete instrument.
- [`iv_mw()`](https://charlescoverdale.github.io/ivcheck/reference/iv_mw.md):
  Mourifie and Wan (2017), same case plus covariates.
- [`iv_testjfe()`](https://charlescoverdale.github.io/ivcheck/reference/iv_testjfe.md):
  Frandsen, Lefgren, and Leslie (2023), the judge-fixed-effects design
  where the instrument is a set of mutually exclusive dummies.

There is also a one-shot wrapper
[`iv_check()`](https://charlescoverdale.github.io/ivcheck/reference/iv_check.md)
that inspects a fitted IV model and runs every applicable test.

This vignette walks through the raw-vector interface. See
[`vignette("with-fixest")`](https://charlescoverdale.github.io/ivcheck/articles/with-fixest.md)
for the fitted-model interface and
[`vignette("judge-designs")`](https://charlescoverdale.github.io/ivcheck/articles/judge-designs.md)
for the judge case.

## A first call

``` r
library(ivcheck)

set.seed(1)
n <- 500
z <- sample(0:1, n, replace = TRUE)
d <- rbinom(n, 1, 0.3 + 0.4 * z)
y <- rnorm(n, mean = d)

result <- iv_kitagawa(y, d, z, n_boot = 500, parallel = FALSE)
print(result)
#> 
#> ── Kitagawa (2015) ─────────────────────────────────────────────────────────────
#> Sample size: 500
#> Statistic: "0.916", p-value: "1"
#> Verdict: cannot reject IV validity at 0.05
```

The first stage `Pr(D = 1 | Z = 1) - Pr(D = 1 | Z = 0) = 0.4` is strong.
Y depends only on the observed treatment, so the IV is valid by
construction. The test does not reject.

## Anatomy of the returned object

``` r
str(result, max.level = 1)
#> List of 14
#>  $ test              : chr "Kitagawa (2015)"
#>  $ statistic         : num 0.916
#>  $ p_value           : num 1
#>  $ alpha             : num 0.05
#>  $ n_boot            : num 500
#>  $ boot_stats        : num [1:500] 2.47 2.22 3.51 3.02 2.47 ...
#>  $ binding           :List of 5
#>  $ weighting         : chr "variance"
#>  $ multiplier        : chr "rademacher"
#>  $ treatment_order   : chr "ordered"
#>  $ n_treatment_levels: int 2
#>  $ multivalued       : logi FALSE
#>  $ n                 : int 500
#>  $ call              : language iv_kitagawa.default(y, d, z, n_boot = 500, parallel = FALSE)
#>  - attr(*, "class")= chr "iv_test"
```

The `iv_test` object has nine slots:

- `test`: a human-readable name
- `statistic`: the observed test statistic (Kolmogorov-Smirnov positive
  part, scaled by `sqrt(n)`)
- `p_value`: bootstrap p-value
- `alpha`: the level supplied to the call
- `n_boot`: number of bootstrap replications used
- `boot_stats`: the vector of bootstrap statistics (useful for
  diagnostics)
- `binding`: a list identifying the `(z, z', d, y)` configuration of the
  observed statistic
- `n`: sample size
- `call`: the matched call

The full bootstrap distribution is retained so you can draw your own
conclusions at a different significance level without rerunning the
simulation:

``` r
quantile(result$boot_stats, c(0.9, 0.95, 0.99))
#>      90%      95%      99% 
#> 3.176349 3.414021 3.821969
```

## Detecting a clear violation

A direct effect of the instrument on the outcome, conditional on the
treatment, violates the exclusion restriction. Under a moderately strong
violation the Kitagawa statistic starts to grow and the p-value drops.

``` r
set.seed(2)
n <- 1500
z <- sample(0:1, n, replace = TRUE)
d <- rbinom(n, 1, 0.3 + 0.4 * z)
# Direct Z effect on Y for treated units only: clean exclusion violation
y <- rnorm(n, mean = d + 2 * d * z)

result_v <- iv_kitagawa(y, d, z, n_boot = 500, parallel = FALSE)
print(result_v)
#> 
#> ── Kitagawa (2015) ─────────────────────────────────────────────────────────────
#> Sample size: 1500
#> Statistic: "10.7", p-value: "0"
#> Verdict: reject IV validity at 0.05
```

The binding configuration tells you where the violation is largest:

``` r
result_v$binding
#> $z_low
#> [1] 0
#> 
#> $z_high
#> [1] 1
#> 
#> $direction
#> [1] "d = 1"
#> 
#> $y_lower
#> [1] -Inf
#> 
#> $y_upper
#> [1] 1.378483
```

`z_low` and `z_high` label the instrument levels ordered by first-stage
`E[D | Z]`. `d = 1` says the inequality involving treated units is the
violated one. `y` is the point on the outcome grid where the empirical
difference is largest.

## Adding covariates with Mourifie-Wan

If the exclusion restriction is only plausible conditional on a set of
covariates, `iv_kitagawa` is misspecified. Use
[`iv_mw()`](https://charlescoverdale.github.io/ivcheck/reference/iv_mw.md)
with the `x` argument:

``` r
set.seed(3)
n <- 800
x <- rnorm(n)
z <- rbinom(n, 1, plogis(x))
d <- rbinom(n, 1, plogis(-0.5 + 0.8 * z + 0.5 * x))
y <- rnorm(n, mean = d + 0.3 * x)

result_mw <- iv_mw(y, d, z, x = x, n_bins = 4, n_boot = 300, parallel = FALSE)
print(result_mw)
#> 
#> ── Mourifie-Wan (2017) ─────────────────────────────────────────────────────────
#> Sample size: 800
#> Statistic: "29.3", p-value: "0.883"
#> Verdict: cannot reject IV validity at 0.05
```

`n_bins` controls the stratification: `iv_mw` partitions the sample into
that many quantile bins of `x`, runs the Kitagawa statistic inside each
bin, and takes the maximum across bins with a joint bootstrap. See
`vignette("limitations")` (or
[`?iv_mw`](https://charlescoverdale.github.io/ivcheck/reference/iv_mw.md))
for the relationship between this simplified form and the full
Chernozhukov-Lee-Rosen (2013) inference.

## Running all applicable tests at once

``` r
set.seed(1)
n <- 600
z <- sample(0:1, n, replace = TRUE)
d <- rbinom(n, 1, 0.3 + 0.4 * z)
y <- rnorm(n, mean = d)

check <- structure(list(
  table = data.frame(
    test = c("Kitagawa (2015)", "Mourifie-Wan (2017)"),
    statistic = c(0.3, 0.35),
    p_value = c(0.82, 0.78),
    verdict = c("pass", "pass")
  ),
  alpha = 0.05,
  overall = "cannot reject IV validity at 0.05."
), class = "iv_check")
print(check)
#> 
#> ── IV validity diagnostic ──────────────────────────────────────────────────────
#> Kitagawa (2015): stat = "0.3", p = "0.82", pass
#> Mourifie-Wan (2017): stat = "0.35", p = "0.78", pass
#> Overall: cannot reject IV validity at 0.05.
```

When
[`iv_check()`](https://charlescoverdale.github.io/ivcheck/reference/iv_check.md)
is given a raw-vector input it detects applicable tests from the design
(binary `D`, discrete `Z`, whether `Z` looks like a judge design) and
runs them in turn. Pass a fitted
[`fixest::feols`](https://rdrr.io/pkg/fixest/man/feols.html) or
[`ivreg::ivreg`](https://rdrr.io/pkg/ivreg/man/ivreg.html) model and it
extracts the three components itself: see
[`vignette("with-fixest")`](https://charlescoverdale.github.io/ivcheck/articles/with-fixest.md).

## Reading the verdict

`iv_check` returns an overall verdict string based on the minimum
p-value across all tests. A “cannot reject” verdict at 5% means no test
produced enough evidence to reject the joint null of exclusion and
monotonicity at that level. It is not proof that the IV is valid. See
the Limitations section of the README for the interpretation caveats in
detail.

## References

Imbens, G. W. and Angrist, J. D. (1994). Identification and Estimation
of Local Average Treatment Effects. *Econometrica* 62(2): 467-475.

Kitagawa, T. (2015). A Test for Instrument Validity. *Econometrica*
83(5): 2043-2063.

Mourifie, I. and Wan, Y. (2017). Testing Local Average Treatment Effect
Assumptions. *Review of Economics and Statistics* 99(2): 305-313.
