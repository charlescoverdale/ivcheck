# ivcheck

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/ivcheck)](https://CRAN.R-project.org/package=ivcheck)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Tests for the identifying assumptions of instrumental variable (IV) models: the local exclusion restriction and monotonicity conditions required for LATE identification.

## Installation

```r
install.packages("ivcheck")

# Or the development version from GitHub
# install.packages("devtools")
devtools::install_github("charlescoverdale/ivcheck")
```

## Why this package?

Every applied IV paper rests on two untestable-looking identifying assumptions: exclusion and monotonicity. In fact, both have testable implications. Kitagawa (2015), Mourifie and Wan (2017), and Frandsen, Lefgren, and Leslie (2023) each provide falsification tests. None of these has a maintained R package. Replication scripts live on author websites and almost no one runs them.

`ivcheck` closes that gap. It provides each test as a named R function, registers S3 methods for fitted `fixest` and `ivreg` objects, and ships a one-shot wrapper that runs every applicable test on a fitted IV model in a single call.

## Quick start

```r
library(ivcheck)
library(fixest)

# Card (1995) proximity-to-college IV
data(card1995)
m <- feols(lwage ~ 1 | educ ~ near_college, data = card1995)

iv_check(m)
#> IV validity diagnostic
#>   Kitagawa (2015):       stat = 0.82, p = 0.34
#>   Mourifie-Wan (2017):   stat = 0.71, p = 0.48
#> Verdict: cannot reject IV validity at 5%.
```

## Functions

| Function | Paper | Case |
|---|---|---|
| `iv_kitagawa()` | Kitagawa (2015, Econometrica) | Binary D, discrete Z |
| `iv_mw()` | Mourifie & Wan (2017, ReStat) | Binary D, discrete Z, with covariates |
| `iv_testjfe()` | Frandsen, Lefgren & Leslie (2023, AER) | Judge IV / mutually exclusive dummy instruments |
| `iv_check()` | This package | One-shot wrapper with auto-detection |
| `iv_power()` | This package | Monte Carlo power simulator |

## References

Kitagawa, T. (2015). A Test for Instrument Validity. *Econometrica*, 83(5), 2043-2063. <https://doi.org/10.3982/ECTA11974>

Mourifie, I. and Wan, Y. (2017). Testing Local Average Treatment Effect Assumptions. *Review of Economics and Statistics*, 99(2), 305-313. <https://doi.org/10.1162/REST_a_00628>

Frandsen, B. R., Lefgren, L. J., and Leslie, E. C. (2023). Judging Judge Fixed Effects. *American Economic Review*, 113(1), 253-277. <https://doi.org/10.1257/aer.20201860>

## Issues

Report bugs or request features at [GitHub Issues](https://github.com/charlescoverdale/ivcheck/issues).

## Keywords

instrumental variables, LATE, causal inference, econometrics, specification testing, exclusion restriction, monotonicity
