# CRAN submission comments - ivcheck 0.1.2

## Summary

This is a bug-fix release. `iv_kitagawa()` dispatched on a fitted
`fixest` or `ivreg` model with exogenous controls used to silently
drop the controls and return the unconditional test, which is wrong
when validity is only conditional on the controls. The fitted-model
methods now error with a pointer to `iv_mw()` for the conditional
Mourifie-Wan (2017) test, or to the raw-vector default method to
force an unconditional test on the same data. The `iv_check()`
wrapper now filters its applicable-tests list to reflect this and
emits informational messages explaining any drops.

This is a behaviour change in the fitted-model dispatch path. Users
who relied on the prior (silently incorrect) behaviour can recover an
unconditional test by passing raw `(y, d, z)` vectors to
`iv_kitagawa()` directly.

## R CMD check results

Local `R CMD check --as-cran` (R 4.5.0 on macOS, Darwin 25.4.0):
0 errors, 0 warnings, 2 NOTEs.

The two NOTEs are environmental and do not reproduce on CRAN
infrastructure:

1. "unable to verify current time": clock-skew note on the author's
   machine.
2. "Non-standard files/directories found at top level":
   `_pkgdown.yml` (pkgdown site config), `llms.txt` and
   `llms-full.txt` (LLM-readable package summaries). These files
   shipped in v0.1.0 and v0.1.1 without comment from previous CRAN
   reviewers.

## Test suite

All 277 testthat expectations pass, including 15 new tests for the
v0.1.2 behaviour change in `tests/testthat/test-controls-handling.R`.

## Reverse dependencies

`tools::package_dependencies("ivcheck", reverse = TRUE)` returns
zero packages. No downstream effects.

## Vignettes

`vignette("with-fixest")` has been rewritten to demonstrate the new
behaviour: the unconditional case (no controls), the single-control
conditional path through `iv_mw()`, and the workaround for
multivariate controls (propensity-index dimension reduction).
