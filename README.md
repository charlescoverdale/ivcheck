# ivcheck

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Tests of the LATE identifying assumptions for instrumental variable (IV) models: three published falsification tests with S3 methods for fitted `fixest` and `ivreg` models, plus a one-shot wrapper. Pure computation, no network.

## What goes in your replication code

```r
library(fixest)
library(ivcheck)

m <- feols(lwage ~ controls | educ ~ near_college, data = card1995)
iv_check(m)
#> IV validity diagnostic
#>   Kitagawa (2015):     stat = 0.01, p = 1.00, pass
#>   Mourifie-Wan (2017): stat = 0.65, p = 0.99, pass
#> Overall: cannot reject IV validity at 0.05.
```

Two added lines, a falsification test the referee is almost guaranteed to ask about, citation-ready output.

## What's in the box

| Function | Paper | Case |
|---|---|---|
| `iv_kitagawa()` | Kitagawa (2015) *Econometrica* / Sun (2023) *JoE* | Binary or multivalued ordered D, discrete Z |
| `iv_mw()` | Mourifie-Wan (2017) *ReStat* | Binary D, discrete Z, with covariates |
| `iv_testjfe()` | Frandsen-Lefgren-Leslie (2023) *AER* | Judge / group IV design |
| `iv_check()` | (this package) | Runs all applicable tests on a fitted IV model |
| `iv_power()` | (this package) | Monte Carlo power curve under a specified violation |

None of these tests had a maintained R package before `ivcheck`. Replication scripts live as Stata modules or bespoke Matlab code on author websites. Most applied IV papers do not run any of them, for reasons of tooling rather than conviction.

## 101: what these tests check

Every IV estimate rests on two assumptions: the **exclusion restriction** (`Z` affects `Y` only through its effect on `D`) and **monotonicity** (no defiers). Together with independence, these are the Imbens-Angrist (1994) conditions that identify the LATE for compliers.

The tests here check necessary *implications* of those assumptions. If the null is true, the joint distribution of `(Y, D, Z)` must satisfy certain inequalities on the conditional CDFs (Kitagawa 2015, Mourifie-Wan 2017) or linearity conditions on per-group means (Frandsen-Lefgren-Leslie 2023). Empirical failure of any of these is evidence that exclusion or monotonicity (or both) has failed.

## What rejection and non-rejection actually mean

**Rejection**: the data show a detectable violation of the testable implication. That rules out at least some forms of exclusion or monotonicity failure. It does not, by itself, localise which assumption failed or in which direction.

**Non-rejection**: there is no detectable violation at the chosen level. It is not proof of validity. The tests have power against violations that show up in the observable distributions but are silent on violations that cancel out across subgroups. Report as "no detectable violation at level 0.05", not "my IV is clean."

## Why trust this implementation

- **Kitagawa statistic matches equation 2.1** of the paper. The sup is taken over the full class of intervals `[y, y']` with `y <= y'`, not just over single lower tails; normalised by the binomial-mixture plug-in standard error of Kitagawa's equation 2.1. Variance-weighted form is the default (`weighting = "variance"`); the unweighted form of equation 2.2 is available via `weighting = "unweighted"`.
- **iv_mw with covariates implements the full Chernozhukov-Lee-Rosen (2013) intersection-bounds framework**: series-regression conditional CDF estimation, heteroscedasticity-robust plug-in standard errors, multiplier bootstrap with adaptive moment selection (Andrews-Soares 2010, `kappa_n = sqrt(log(log(n)))`). Without covariates, `iv_mw` reduces exactly to the variance-weighted Kitagawa test (unit-tested).
- **iv_testjfe null distribution matches `chi^2_{K-2}` to Monte Carlo precision.** At `K=20`, `N=3000` (200 replications): empirical mean 18.01 vs target 18.0, variance 35.1 vs target 36.0, 95th percentile 29.4 vs target 28.9, empirical size 6.5% vs nominal 5%. A `method = "bootstrap"` option is available for small-K robustness.
- **All DOIs crossref-verified.** The audit caught a silent bug: Mourifie-Wan's DOI had been cited as `10.1162/REST_a_00628` (which resolves to a different paper entirely). The correct DOI is `00622`. Fixed before first submission.
- **R CMD check --as-cran**: 0 errors, 0 warnings, 0 notes. 87 unit tests covering structure, invariants, known-value cases, edge cases, and end-to-end S3 dispatch against `fixest` and `ivreg`.

## Installation

```r
# Once accepted by CRAN
install.packages("ivcheck")

# Development version from GitHub
# install.packages("devtools")
devtools::install_github("charlescoverdale/ivcheck")
```

## Walkthrough

Output lines prefixed with `#>` are what the console shows.

### A single test on raw vectors

```r
library(ivcheck)

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
#> Statistic: 0.04, p-value: 0.91
#> Verdict: cannot reject IV validity at 0.05
```

The bootstrap p-value comes from the multiplier resampling of Kitagawa (2015) section 3.2. With `parallel = TRUE` (default) replications run across cores on POSIX.

### With covariates (Mourifie-Wan)

```r
x <- rnorm(n)
mw <- iv_mw(y, d, z, x = x, n_boot = 500)
print(mw)
```

`iv_mw()` with covariates estimates `F(y, d | X = x, Z = z)` by cubic-polynomial series regression, computes heteroscedasticity-robust standard errors, and takes the sup of the studentised positive-part over a grid of `(y, x)` points. Critical values use adaptive moment selection. Without covariates it reduces to Kitagawa.

### Judge designs (Frandsen-Lefgren-Leslie)

```r
set.seed(1)
n <- 2000
judge <- sample.int(20, n, replace = TRUE)
d <- rbinom(n, 1, 0.3 + 0.02 * judge)
y <- rnorm(n, mean = d)

jfe <- iv_testjfe(y, d, judge, n_boot = 500)
```

Designs where the instrument is a set of mutually exclusive dummies (judge, caseworker, examiner) need a purpose-built test. `iv_testjfe()` fits a weighted-LS regression of per-judge `mu_j` on per-judge `p_j` and tests the implied linearity via chi-squared with `K - 2` degrees of freedom (default) or multiplier bootstrap (`method = "bootstrap"`).

### One-shot diagnostic on a fitted model

```r
library(fixest)

df <- data.frame(z = z, d = d, y = y, x = x)
m  <- feols(y ~ x | d ~ z, data = df)

iv_check(m, n_boot = 500)
```

`iv_check()` detects applicable tests from the model structure (binary D, discrete Z, judge dummies) and runs all of them. Works identically on `ivreg::ivreg()` objects.

### Power planning

```r
pw <- iv_power(y, d, z, method = "kitagawa", n_sims = 200)
```

Simulates data under a parametric exclusion violation and reports rejection probability at a grid of deviation sizes. Useful when choosing between candidate tests on the same design.

## Example: end-to-end with Card (1995)

```r
library(ivcheck)
library(fixest)

data(card1995)   # bundled
m <- feols(
  lwage ~ age + married + black + south | college ~ near_college,
  data = card1995
)

iv_check(m, n_boot = 1000)
#> IV validity diagnostic
#>   Kitagawa (2015):     stat = 0.01, p = 1.00, pass
#>   Mourifie-Wan (2017): stat = 32.4, p = 0.67, pass
#> Overall: cannot reject IV validity at 0.05.
```

Card's proximity-to-college instrument passes both tests. This does not prove validity (no test can) but it rules out detectable violations at the 5% level. Report the test statistic alongside the IV estimate, cite Kitagawa (2015) or Mourifie-Wan (2017), move on.

## Limitations

Read before using in published work.

### Scope (v0.1.0 does not cover)

- **Continuous instruments.** All three tests require a discrete `Z`. For continuous instruments, discretise into quantile bins (quartiles or quintiles) before passing to `iv_kitagawa` or `iv_mw`. This is standard practice in the applied literature and costs a modest amount of power. A formal nonparametric continuous-`Z` extension (e.g. the local-moment-inequality framework of Andrews and Shi 2013) is on the v0.2.0 roadmap.
- **Fuzzy regression discontinuity.** FRD has its own testable implications at the cutoff [@arai2022]. Handling them requires different infrastructure (running variable, bandwidth selection, bias correction) that does not fit the current `iv_test`-on-fitted-IV-model spine; a dedicated `iv_frd()` function is planned for v0.2.0.
- **`iv_mw` with covariates under weights.** The `weights` argument is fully implemented for the unconditional path and for the FLL judge test, but the CLR series-regression path for `iv_mw` with covariates does not yet propagate weights. Use `iv_kitagawa` directly with weights if you need weighted inference without covariates, or wait for v0.1.1.

### Faithfulness to the published tests

All three test functions implement the asymptotic form their authors publish. Under binary treatment with moderate `K`, `iv_testjfe()` is algebraically the FLL (2023) overidentification test: the weighted sum of squared residuals from the linear fit `mu_j = alpha + beta * p_j` is, via WLS, the quadratic form implied by pairwise Wald LATE equality. `iv_testjfe()` also returns the `pairwise_late` matrix and `worst_pair` as diagnostic output.

The v0.1.0 scope limitation for `iv_testjfe()` is **multivalued treatment** only: Frandsen, Lefgren, and Leslie (2023) Section 4 extends to multi-valued `D`, and this package does not yet support it. `iv_testjfe()` errors on non-binary `D`. Use the Stata `testjfe` module until the v0.2.0 port lands.

`iv_kitagawa()` and `iv_mw()` are faithful implementations of their published forms (variance-weighted KS per Kitagawa 2015 §4 and full CLR intersection-bounds inference with adaptive moment selection per Mourifie-Wan 2017 and Andrews-Soares 2010 respectively).

### Interpretation

- **Non-rejection is not proof of validity.** The tests have power against violations in the observable conditional distributions but are silent on violations that cancel out.
- **Kitagawa vs Mourifie-Wan with covariates.** If the exclusion restriction is only plausible conditional on `X`, run `iv_mw` with `x`. Running `iv_kitagawa` unconditionally on an X-dependent design can give spurious non-rejection.
- **Many-instrument / judge regimes.** For 20+ judge levels, prefer `iv_testjfe` over `iv_kitagawa`; the KS test loses power rapidly as `|Z|` grows.
- **Bootstrap size.** `n_boot = 1000` (default) is fine for publication-grade p-values. Cut to 200 for exploration; raise to 5000 if reporting p-values to three decimal places.

### Language implementations

As of 2026-04-19 there is no equivalent package on PyPI. Kitagawa (2015) has Matlab supplementary code; Mourifie-Wan (2017) has a Stata `clrtest` implementation; Frandsen-Lefgren-Leslie (2023) has the Stata `testjfe` SSC module. `ivcheck` is the first R-native implementation of this family.

## Planned for future versions

- `iv_hm()`: Huber-Mellace (2015, *ReStat*) mean-based moment-inequality form, complementary to the CDF form of Kitagawa
- `iv_frd()`: Arai, Hsu, Kitagawa, Mourifie, Wan (2022, *QE*) fuzzy regression discontinuity
- Continuous-instrument extension via Andrews and Shi (2013) conditional-moment-inequality inference
- Weighted inference in the `iv_mw` conditional (x) series-regression path
- Rcpp fast path for the interval-sup multiplier bootstrap
- Full flexible-basis FLL restricted-LS test with Andrews-Soares bounded-slope moment selection

## Functions

| Function | Purpose |
|---|---|
| `iv_kitagawa()` | Kitagawa (2015) variance-weighted KS test |
| `iv_mw()` | Mourifie-Wan (2017) conditional-inequality test with adaptive GMS |
| `iv_testjfe()` | Frandsen-Lefgren-Leslie (2023) judge / group IV test |
| `iv_check()` | Wrapper: runs all applicable tests on a fitted IV model |
| `iv_power()` | Monte Carlo power curve under a specified violation |

## Citation

Cite both the package and the underlying paper(s) for the test you use. Package citation:

```r
citation("ivcheck")
```

Underlying paper citations (DOIs verified via crossref.org):

| Function | Reference | DOI |
|---|---|---|
| `iv_kitagawa()` | Kitagawa, T. (2015). A Test for Instrument Validity. *Econometrica* 83(5): 2043-2063. | [10.3982/ECTA11974](https://doi.org/10.3982/ECTA11974) |
| `iv_mw()` | Mourifie, I. and Wan, Y. (2017). Testing Local Average Treatment Effect Assumptions. *Review of Economics and Statistics* 99(2): 305-313. | [10.1162/REST_a_00622](https://doi.org/10.1162/REST_a_00622) |
| `iv_testjfe()` | Frandsen, B. R., Lefgren, L. J., Leslie, E. C. (2023). Judging Judge Fixed Effects. *American Economic Review* 113(1): 253-277. | [10.1257/aer.20201860](https://doi.org/10.1257/aer.20201860) |

Foundational LATE framework: Imbens and Angrist (1994), *Econometrica* 62(2): 467-475. Intersection-bounds inference inside `iv_mw`: Chernozhukov, Lee, and Rosen (2013), *Econometrica* 81(2): 667-737. Adaptive moment selection: Andrews and Soares (2010), *Econometrica* 78(1): 119-157.

## Related packages

| Package | What it covers |
|---|---|
| [`fixest`](https://cran.r-project.org/package=fixest) | Fast IV estimation via `feols(y ~ x \| d ~ z)` (upstream from `ivcheck`) |
| [`ivreg`](https://cran.r-project.org/package=ivreg) | 2SLS with Wu-Hausman, Sargan, weak-IV F (upstream from `ivcheck`) |
| [`ivmodel`](https://cran.r-project.org/package=ivmodel) | k-class estimators, weak-IV robust CIs, sensitivity analysis |
| [`ivDiag`](https://cran.r-project.org/package=ivDiag) | Effective F, Anderson-Rubin, valid-t, local-to-zero tests |

`ivcheck` complements rather than competes: `fixest` or `ivreg` does the estimation, `ivDiag` does post-estimation diagnostics, and `ivcheck` does LATE-assumption falsification.

## Issues and requests

Report bugs or request additional tests at [GitHub Issues](https://github.com/charlescoverdale/ivcheck/issues). Pull requests implementing additional IV-validity tests from the literature are welcome; please include a reference to the original paper and a reproduction test against its empirical example.

## Keywords

instrumental variables, LATE, causal inference, exclusion restriction, monotonicity, specification testing, falsification, judge IV, Kitagawa test, Mourifie-Wan test, FLL test, econometrics.
