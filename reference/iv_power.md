# Monte Carlo power curve for IV-validity tests

Simulates data under a user-specified deviation from validity and
estimates the rejection probability of the chosen test at each deviation
size. Useful for sample-size planning and for benchmarking different
tests on the same design.

## Usage

``` r
iv_power(
  y,
  d,
  z,
  method = c("kitagawa", "mw", "testjfe"),
  alpha = 0.05,
  n_sims = 500,
  delta_grid = NULL,
  n_boot = 200,
  parallel = TRUE,
  ...
)
```

## Arguments

- y, d, z:

  Observed data used to anchor the DGP (sample size, cell counts,
  empirical first-stage).

- method:

  Which test to benchmark. One of `"kitagawa"`, `"mw"`, or `"testjfe"`.

- alpha:

  Significance level.

- n_sims:

  Number of Monte Carlo simulations per deviation.

- delta_grid:

  Numeric vector of deviation sizes to evaluate. If `NULL`, defaults to
  `seq(0, 0.3, by = 0.05)`.

- n_boot:

  Number of bootstrap replications per simulation (for tests that use
  bootstrap). Default 200, which trades some Monte Carlo noise for
  tractable runtime.

- parallel:

  Logical. Run simulations in parallel on POSIX systems via
  [parallel::mclapply](https://rdrr.io/r/parallel/mclapply.html).
  Default `TRUE`.

- ...:

  Further arguments passed to the underlying test.

## Value

A data frame with columns `delta` (deviation size) and `power`
(estimated rejection probability at level `alpha`).

## Details

The deviation is parameterised as the size of a **D-specific direct
effect of the instrument on the outcome** (a clean exclusion violation
that the Kitagawa and Mourifie-Wan tests are designed to detect).
Specifically, the simulated outcome is
`Y = mu_hat[D + 1] + delta * sigma_hat * D * (Z - Z_low) + noise`, so
`delta = 0` corresponds to the null and larger values produce larger
violations of the testable inequality for the d = 1 cells. The simulator
preserves the observed sample size, first-stage propensities, and
outcome scale.

## Examples

``` r
# \donttest{
# Headline power curve for a small-N design
set.seed(1)
n <- 300
z <- sample(0:1, n, replace = TRUE)
d <- rbinom(n, 1, 0.3 + 0.4 * z)
y <- rnorm(n, mean = d)
iv_power(y, d, z, method = "kitagawa", n_sims = 50, n_boot = 100)
#>   delta power
#> 1  0.00     0
#> 2  0.05     0
#> 3  0.10     0
#> 4  0.15     0
#> 5  0.20     0
#> 6  0.25     0
#> 7  0.30     0
# }
```
