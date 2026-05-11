# ivcheck 0.1.2

## Bug fix (behaviour change)

* `iv_kitagawa()` dispatched on a fitted `fixest` or `ivreg` model with
  exogenous controls used to silently drop the controls and return the
  unconditional test. This was wrong when validity is only conditional
  on X. The fitted-model methods now error and direct the user to
  `iv_mw()` for the conditional Mourifie-Wan (2017) test, or to the
  raw-vector default method to force an unconditional test on the same
  data.
* `iv_check()` now filters its applicable-tests list by the presence
  and dimensionality of exogenous controls. Kitagawa is dropped when
  any control is present; Mourifie-Wan is dropped when more than one
  control is present (multivariate X is planned for v0.2.0 via a
  tensor-product series basis). When a test is dropped, `iv_check()`
  emits an informational message explaining why.
* `iv_check()` no longer falls back to running Mourifie-Wan
  unconditionally when handed a model with multivariate controls. That
  downgrade was the same silent-drop bug as the Kitagawa one.

## Documentation

* `vignette("with-fixest")` rewritten to demonstrate the new behaviour:
  the unconditional case, the single-control conditional path through
  `iv_mw()`, and the multivariate workaround (propensity-index
  dimension reduction).

# ivcheck 0.1.1

CRAN resubmission addressing reviewer feedback (Benjamin Altmann, 2026-04-21):

* `iv_power()` example: replaced `\dontrun{}` with `\donttest{}` so the
  example is executed by default on CRAN checks (it exceeds the 5-second
  limit for inline examples due to 50 simulations x 100 bootstrap
  iterations).
* `judge-designs` vignette: save and restore `par()` around the two-panel
  diagnostic plot so user `par` settings are not permanently altered.

# ivcheck 0.1.0

First CRAN submission. Three falsification tests for the identifying
assumptions of instrumental variable estimation of the local average
treatment effect, plus a one-shot wrapper.

For known deviations from the published tests (richer family of
ordered-multivalued inequalities than Sun 2023 Lemma 2.1; linearity-test
form of Frandsen, Lefgren, and Leslie 2023; empirical `se_floor = 0.15`
versus Kitagawa 2015's 0.05-0.10 recommendation), see the "Notes on
fidelity" sub-section of the README and the companion R Journal paper
at `paper/rj/paper.Rmd`.

## New functions

* `iv_kitagawa()`: Kitagawa (2015) variance-weighted Kolmogorov-Smirnov
  test of instrument validity under binary (or ordered multivalued,
  per Sun 2023) treatment and discrete instrument. Supports survey
  weights, user-tunable `se_floor` trimming, and both the unweighted
  and variance-weighted forms. Returns the binding `(z, z', d, y)`
  interval so users can see which outcome region carries the violation.
* `iv_mw()`: Mourifie and Wan (2017) conditional moment-inequality
  reformulation. Without covariates, reduces exactly to the variance-
  weighted Kitagawa test (unit-tested). With a covariate, uses the full
  Chernozhukov-Lee-Rosen (2013) intersection-bounds inference with
  series-regression conditional CDF estimation, heteroscedasticity-
  robust plug-in standard errors, and Andrews-Soares (2010) adaptive
  moment selection.
* `iv_testjfe()`: Frandsen, Lefgren, and Leslie (2023) joint test for
  exclusion and monotonicity in judge-fixed-effects designs. Asymptotic
  chi-squared and multiplier-bootstrap p-values. Supports binary and
  multivalued ordered treatment, survey weights, and covariate
  residualisation. Returns the pairwise-Wald-LATE matrix and the worst
  deviating pair as diagnostics.
* `iv_check()`: one-shot wrapper that runs all applicable tests on a
  fitted IV model and returns a tidy summary with overall verdict.
* `iv_power()`: Monte Carlo power simulator anchored to the observed
  first-stage propensities, for sample-size planning and benchmarking.

## S3 methods

* Dispatch on `fixest::feols` and `ivreg::ivreg` fitted IV models.
* `print()`, `plot()`, `summary()` methods for `iv_test` and `iv_check`.
* Optional `modelsummary` glue registered on load.

## Bundled data

* `card1995`: cleaned 3,003-row extract of the Card (1995) National
  Longitudinal Survey of Young Men data, distributed via the
  `wooldridge` CRAN package. Includes a binary `college` indicator
  (`educ >= 16`) suitable for binary-treatment tests.

## Documentation

* Three vignettes: `getting-started`, `judge-designs`, `with-fixest`.
* A companion R Journal paper in `paper/rj/paper.Rmd`.
* A replication script in `paper/replicate.R` reproducing every
  numerical claim in the paper.

## Scope (v0.1.0 does not cover)

* Continuous instruments (discretise into bins).
* Fuzzy regression discontinuity (planned v0.2.0).
* Fixed-effects IV models: pre-demean within each FE cell and pass
  vectors to the default method.
* Multivariate conditioning in `iv_mw`: supports a single covariate;
  multivariate planned v0.2.0.
* Sun (2023) unordered multivalued D: ordered extension is implemented.
