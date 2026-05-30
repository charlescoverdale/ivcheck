# CRAN submission comments - ivcheck 0.1.2

## Resubmission

This is a resubmission. The 2026-05-30 pre-test flagged one NOTE:

    Non-standard files/directories found at top level:
      '_pkgdown.yml' 'llms-full.txt' 'llms.txt'

These three files are now listed in `.Rbuildignore` so they are no
longer shipped in the package tarball. No other changes.

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
0 errors, 0 warnings, 1 NOTE.

The single NOTE is environmental and does not reproduce on CRAN
infrastructure:

1. "unable to verify current time": clock-skew note on the author's
   machine.

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
