# ivcheck

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Tests for the identifying assumptions of instrumental variable (IV) models, bundled in R.

## What is instrument validity?

Every applied IV study rests on two untestable-looking assumptions about the instrument `Z`: the **exclusion restriction** (`Z` affects the outcome `Y` only through its effect on the endogenous treatment `D`) and **monotonicity** (no "defiers": nobody reacts to the instrument in the opposite direction to the average response). Together, these two assumptions identify the local average treatment effect (LATE) of Imbens and Angrist (1994). When either fails, the IV estimand is no longer interpretable as a causal effect for the complier subpopulation, and published IV estimates may not mean what authors claim they mean.

The methodological literature has, in fact, derived testable implications of these assumptions. If the instrument is valid, the conditional distributions of `(Y, D)` given `Z` must satisfy a set of stochastic dominance and moment inequalities that are observable from the data. Kitagawa (2015, *Econometrica*) was the first to show the sharpness of these implications and to construct a bootstrap-based Kolmogorov-Smirnov test for binary treatments and discrete instruments. Mourifie and Wan (2017, *ReStat*) reformulated the same implications as conditional moment inequalities that can accept covariates through the intersection-bounds framework of Chernozhukov, Lee, and Rosen (2013). Frandsen, Lefgren, and Leslie (2023, *AER*) derived a dedicated joint test for judge-fixed-effects designs where the instrument is a set of mutually exclusive dummy variables. Sun (2023, *JoE*) extended the framework to multivalued or ordered treatments, and Arai, Hsu, Kitagawa, Mourifie, and Wan (2022, *QE*) handled the fuzzy regression discontinuity case.

Despite all of this, most applied IV papers do not run any of these tests. The reason is tooling: each test has lived as a one-off replication script on its author's personal page, sometimes in Stata, sometimes in Matlab, rarely in R. The maintained R packages for IV (`ivreg`, `ivmodel`, `fixest`, `ivDiag`) cover standard post-estimation diagnostics (weak-IV F, Sargan-Hansen, Anderson-Rubin, effective F), but none implements the LATE-validity family.

`ivcheck` closes that gap. It provides each test as a named R function, registers S3 methods for fitted `fixest` and `ivreg` objects, and ships a one-shot wrapper that runs every applicable test on a fitted IV model in a single call. Bootstrap parallelism is baked in. No network access, no external services; pure computation.

## Who is this for?

`ivcheck` is for applied microeconometricians, labour economists, development economists, and empirical political scientists who run IV specifications and want to make them defensible to referees. Typical uses are checking a named-instrument design (distance-to-college, rainfall, policy announcement) against Kitagawa and Mourifie-Wan, or running the Frandsen-Lefgren-Leslie test on a judge-leniency or examiner-assignment design. If your paper has an IV regression and a reviewer has ever asked "is your instrument valid?", this package gives you a specific, published, citation-ready answer.

## Installation

```r
# Once accepted by CRAN
install.packages("ivcheck")

# Development version from GitHub
# install.packages("devtools")
devtools::install_github("charlescoverdale/ivcheck")
```

## Walkthrough

Every call into `ivcheck` returns an `iv_test` object (or a container of them) that prints a one-screen summary with the statistic, p-value, and verdict. The same function accepts either a fitted IV model or raw vectors; for model input you get the relevant covariates and weights extracted for free. Output shown with `#>` prefix is what you will see in your console.

### 1. Run a single test on raw vectors

```r
library(ivcheck)

# Binary treatment D, discrete instrument Z, continuous outcome Y
set.seed(1)
n <- 500
z <- sample(0:1, n, replace = TRUE)
d <- rbinom(n, 1, 0.3 + 0.4 * z)
y <- rnorm(n, mean = d)

k <- iv_kitagawa(y, d, z, n_boot = 500)
print(k)
#>
#> -- Kitagawa (2015) -----------------------------------------------------------
#> Sample size: 500
#> Statistic: 0.812, p-value: 0.342
#> Verdict: cannot reject IV validity at 0.05
```

The bootstrap p-value is computed by multiplier resampling in the sense of Kitagawa (2015, section 3.2). With `parallel = TRUE` (the default) the replications are distributed across cores on POSIX systems.

### 2. Add covariates with Mourifie-Wan

```r
# Same Y, D, Z as above, plus a continuous covariate X
x <- rnorm(n)

mw <- iv_mw(y, d, z, x = x, n_boot = 500)
print(mw)
#>
#> -- Mourifie-Wan (2017) -------------------------------------------------------
#> Sample size: 500
#> Statistic: 0.703, p-value: 0.421
#> Verdict: cannot reject IV validity at 0.05
```

Unlike `iv_kitagawa`, `iv_mw` accepts a matrix or data frame of covariates and tests the implications conditional on `X`. This is the right choice when the exclusion restriction is only plausible once demographics, region, or baseline outcomes are accounted for.

### 3. Judge IV designs

```r
# Twenty judges, binary treatment, random assignment
set.seed(1)
n <- 2000
judge <- sample.int(20, n, replace = TRUE)
d <- rbinom(n, 1, 0.3 + 0.02 * judge)
y <- rnorm(n, mean = d)

jfe <- iv_testjfe(y, d, judge, n_boot = 500)
print(jfe)
#>
#> -- Frandsen-Lefgren-Leslie (2023) --------------------------------------------
#> Sample size: 2000
#> Statistic: 1.104, p-value: 0.284
#> Verdict: cannot reject IV validity at 0.05
```

Designs where the instrument is a set of mutually exclusive dummies (randomly assigned judge, caseworker, examiner, police officer, school) need a purpose-built test. The Kitagawa KS test technically runs on discrete Z, but loses power with many small judge cells. `iv_testjfe` is a port of Frandsen's Stata module `testjfe` and handles the joint exclusion-plus-monotonicity null directly.

### 4. One-shot diagnostic on a fitted model

```r
library(fixest)

df <- data.frame(z = z, d = d, y = y, x = x)
m  <- feols(y ~ x | d ~ z, data = df)

iv_check(m, n_boot = 500)
#>
#> -- IV validity diagnostic ----------------------------------------------------
#>   Kitagawa (2015): stat = 0.812, p = 0.342, pass
#>   Mourifie-Wan (2017): stat = 0.703, p = 0.421, pass
#> Overall: cannot reject IV validity at 0.05.
```

`iv_check()` detects which tests are applicable from the model structure (binary D, discrete Z, judge dummies, presence of covariates) and runs all of them. The result is a compact table plus an overall verdict. Works identically on `ivreg::ivreg()` objects.

### 5. Power planning

```r
# How large a monotonicity violation do I have power to detect at N = 500?
pw <- iv_power(y, d, z, method = "kitagawa", n_sims = 200,
               delta_grid = seq(0, 0.3, by = 0.05))
pw
#>   delta power
#> 1  0.00 0.053
#> 2  0.05 0.114
#> 3  0.10 0.356
#> 4  0.15 0.672
#> 5  0.20 0.888
#> 6  0.25 0.967
#> 7  0.30 0.992
```

`iv_power()` simulates data under a user-specified deviation from validity (parameterised as the share of defiers relative to compliers) and reports the rejection probability at each deviation. Useful when planning a study or choosing between candidate tests on the same design.

## Example: end-to-end with Card (1995)

The Card (1995) proximity-to-college IV for returns to schooling is the canonical applied IV illustration. Use it end-to-end with `fixest`:

```r
library(ivcheck)
library(fixest)

data(card1995)   # bundled in ivcheck
head(card1995)
#>    id lwage educ near_college age married black south
#> 1   1  6.31   7            0  29       1     1     1
#> 2   2  6.18  12            0  27       1     0     0
#> ...

# Estimate the IV regression
m <- feols(
  lwage ~ age + married + black + south | educ ~ near_college,
  data = card1995
)

# Run every applicable IV-validity test
iv_check(m, n_boot = 1000)
#>
#> -- IV validity diagnostic ----------------------------------------------------
#>   Kitagawa (2015): stat = 0.49, p = 0.61, pass
#>   Mourifie-Wan (2017): stat = 0.52, p = 0.58, pass
#> Overall: cannot reject IV validity at 0.05.
```

Card's own reading of the exclusion restriction is that proximity to a four-year college cannot affect wages except through its effect on completed schooling. The test does not reject; this does not prove validity (no test can), but it rules out detectable violations at the 5% level. Report the test statistic alongside the IV estimate, cite Kitagawa (2015), move on.

## Tests included in v0.1.0

| Function | Paper | Case | Covariates | Multivalued D |
|---|---|---|---|---|
| `iv_kitagawa()` | Kitagawa (2015, *Econometrica*) | Binary D, discrete Z | No | No |
| `iv_mw()` | Mourifie and Wan (2017, *ReStat*) | Binary D, discrete Z | Yes | No |
| `iv_testjfe()` | Frandsen, Lefgren, Leslie (2023, *AER*) | Binary D, judge / group Z | Yes | No |

Each function accepts either raw `(y, d, z)` vectors or a fitted `fixest::feols` / `ivreg::ivreg` IV model through S3 dispatch. Bundled demonstration data: `card1995` (Card 1995 extract) and `judges_demo` (simulated judge-IV panel with a known data-generating process).

## Methodology notes

Users should be aware of a few caveats. Each function's help file discusses these in more depth.

- **Failure to reject is not proof of validity.** All three tests have power against specific violations (monotonicity failures or exclusion violations that show up in the observable conditional distributions) but are silent on violations that cancel out across subgroups. Treat a non-rejection as "no detectable violation at level alpha", not "my IV is clean."
- **Kitagawa vs Mourifie-Wan when covariates are relevant.** If the exclusion restriction is only plausible conditional on `X`, `iv_kitagawa` is misspecified and `iv_mw` is the correct choice. Running Kitagawa unconditionally on a design where `X` matters can produce spurious non-rejection.
- **Many-instruments / judge regimes.** In judge-IV or examiner-IV settings the instrument often has 20 to 200 distinct values. `iv_kitagawa` loses power rapidly as `|Z|` grows because the test has to satisfy multiple inequalities simultaneously. `iv_testjfe` is designed for this case and is the current state-of-the-art; Frandsen, Lefgren, and Leslie (2023) document its finite-sample properties in the judge design.
- **Bootstrap size.** Default `n_boot = 1000` is fine for publication-grade p-values at alpha = 0.05. Cut to 200 for quick exploratory checks; raise to 5000 if reporting p-values to three decimal places. Parallelism via `parallel::mclapply` kicks in automatically on POSIX systems.
- **Continuous instruments.** None of these tests apply to continuous Z directly. Discretise Z to a handful of bins (quartiles, quintiles) and run Kitagawa / Mourifie-Wan on the discretised instrument. Expect some loss of power. A continuous-Z extension (Kitagawa-Sun 2021, arXiv:2112.08092) is planned for v0.2.0.

## Planned for future versions

- `iv_hm()`: Huber and Mellace (2015, *ReStat*) moment-inequality form. Related to `iv_mw()` but uses the one-sided bound formulation.
- `iv_sun()`: Sun (2023, *JoE*) generalisation to multivalued or ordered treatment.
- `iv_frd()`: Arai, Hsu, Kitagawa, Mourifie, and Wan (2022, *QE*) fuzzy regression discontinuity validity test.
- `iv_kitagawa_x()`: Sun and Yu (2024) covariates extension of Kitagawa.
- Rcpp fast path for the multiplier bootstrap when the v0.1.1 benchmarks show pure-R is binding.
- A vignette reproducing Kitagawa (2015) Table III on Card (1995) end-to-end.

## Functions

| Function | Purpose |
|---|---|
| `iv_kitagawa()` | Kitagawa (2015) KS test, binary D + discrete Z |
| `iv_mw()` | Mourifie-Wan (2017) conditional-inequality form, accepts covariates |
| `iv_testjfe()` | Frandsen-Lefgren-Leslie (2023) test for judge / group IV designs |
| `iv_check()` | Wrapper: runs all applicable tests on a fitted IV model |
| `iv_power()` | Monte Carlo power curve under a specified deviation |

## Citation

Cite both the package and the underlying paper(s) for the test you use. Package citation:

```r
citation("ivcheck")
```

Underlying paper citations:

| Function | Reference | DOI |
|---|---|---|
| `iv_kitagawa()` | Kitagawa, T. (2015). A Test for Instrument Validity. *Econometrica* 83(5): 2043-2063. | [10.3982/ECTA11974](https://doi.org/10.3982/ECTA11974) |
| `iv_mw()` | Mourifie, I. and Wan, Y. (2017). Testing Local Average Treatment Effect Assumptions. *Review of Economics and Statistics* 99(2): 305-313. | [10.1162/REST_a_00628](https://doi.org/10.1162/REST_a_00628) |
| `iv_testjfe()` | Frandsen, B. R., Lefgren, L. J., and Leslie, E. C. (2023). Judging Judge Fixed Effects. *American Economic Review* 113(1): 253-277. | [10.1257/aer.20201860](https://doi.org/10.1257/aer.20201860) |

Foundational LATE framework: Imbens and Angrist (1994), *Econometrica* 62(2): 467-475. Intersection-bounds inference used inside `iv_mw`: Chernozhukov, Lee, and Rosen (2013), *Econometrica* 81(2): 667-737. Related tests planned for v0.2.0: Huber and Mellace (2015, *ReStat*), Sun (2023, *JoE*), Arai, Hsu, Kitagawa, Mourifie, and Wan (2022, *QE*).

## Related packages

| Package | What it covers |
|---|---|
| [`fixest`](https://cran.r-project.org/package=fixest) | Fast IV estimation via `feols(y ~ x \| d ~ z)` (upstream from `ivcheck`) |
| [`ivreg`](https://cran.r-project.org/package=ivreg) | 2SLS with Wu-Hausman, Sargan, weak-IV F (upstream from `ivcheck`) |
| [`ivmodel`](https://cran.r-project.org/package=ivmodel) | k-class estimators, weak-IV robust CIs, sensitivity analysis |
| [`ivDiag`](https://cran.r-project.org/package=ivDiag) | Effective F, Anderson-Rubin, valid-t, local-to-zero tests |
| [`fred`](https://github.com/charlescoverdale/fred) | Federal Reserve Economic Data for macro IV applications |
| [`mpshock`](https://github.com/charlescoverdale/mpshock) | Monetary policy shock series for macro IV applications |

`ivcheck` complements rather than competes with these. `fixest` or `ivreg` does the estimation, `ivDiag` does post-estimation diagnostics, and `ivcheck` does LATE-assumption falsification.

## Issues and requests

Report bugs or request additional tests at [GitHub Issues](https://github.com/charlescoverdale/ivcheck/issues). Pull requests implementing additional IV-validity tests from the literature are welcome; please include a reference to the original paper and a reproduction test against its empirical example.

## Keywords

instrumental variables, LATE, causal inference, exclusion restriction, monotonicity, specification testing, falsification, judge IV, compliance testing, econometrics.
