# ivcheck 0.1.0

First CRAN submission. Three falsification tests for the identifying
assumptions of instrumental variable estimation of the local average
treatment effect, plus a one-shot wrapper.

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
