#' Kitagawa (2015) test for instrument validity
#'
#' Tests the joint implication of the local exclusion restriction and the
#' local monotonicity condition in a binary-treatment, discrete-instrument
#' setting. The null hypothesis is that the instrument is valid. Under the
#' null, the conditional distributions of the outcome given treatment
#' status, evaluated at each level of the instrument, must satisfy a set
#' of stochastic dominance inequalities. Rejection is evidence that at
#' least one of exclusion or monotonicity fails.
#'
#' @param object For the default method: a numeric outcome vector.
#'   For the `fixest` and `ivreg` methods: a fitted instrumental variable
#'   model from [fixest::feols] or `ivreg::ivreg()`.
#' @param d Binary 0/1 treatment vector (default method only).
#' @param z Discrete instrument (numeric or factor, default method only).
#' @param n_boot Number of multiplier-bootstrap replications. Default 1000.
#' @param alpha Significance level for the returned verdict. Default 0.05.
#' @param weighting Test-statistic weighting. `"variance"` (default)
#'   divides each pointwise difference by its plug-in standard error
#'   estimator before taking the sup, as in Kitagawa (2015) section 4.
#'   `"unweighted"` uses the raw positive-part KS of section 3. The
#'   two are asymptotically equivalent at the boundary of the null;
#'   `"variance"` has better finite-sample power when instrument cells
#'   have unequal sizes.
#' @param weights Optional survey weights. Not yet implemented; reserved
#'   for v0.2.0.
#' @param parallel Logical. Run bootstrap replications in parallel on
#'   POSIX systems via [parallel::mclapply]. Default `TRUE`.
#' @param ... Further arguments passed to methods.
#'
#' @return An object of class `iv_test` with elements:
#'   \item{test}{Character, always `"Kitagawa (2015)"`.}
#'   \item{statistic}{Numeric test statistic (Kolmogorov-Smirnov
#'     positive-part, scaled by sqrt(n)).}
#'   \item{p_value}{Bootstrap p-value.}
#'   \item{alpha}{Supplied significance level.}
#'   \item{n_boot}{Number of bootstrap replications used.}
#'   \item{boot_stats}{Numeric vector of bootstrap test statistics.}
#'   \item{binding}{List identifying the binding `(z, z', d, y)`
#'     configuration of the observed statistic.}
#'   \item{n}{Sample size.}
#'   \item{call}{Matched call.}
#'
#' @details
#' Let `F(y, d | z) = Pr(Y <= y, D = d | Z = z)` denote the empirical
#' joint CDF of outcome and treatment conditional on the instrument.
#' Kitagawa (2015) shows that under IV validity there exists an ordering
#' of the instrument levels such that for every pair `(z, z')` with `z`
#' ordered above `z'`, and for every `y`, both
#' `F(y, 1 | z') <= F(y, 1 | z)` and `F(y, 0 | z) <= F(y, 0 | z')`. The
#' one-sided Kolmogorov-Smirnov statistic
#' `T_n = sqrt(n) * max over (z, z', d, y) of [F(y, d | z) - F(y, d | z')]^+`
#' tests this implication in its most-violated direction. The
#' distribution of `T_n` is obtained by multiplier bootstrap with
#' Rademacher weights as in Kitagawa (2015) section 3.2.
#'
#' @references
#' Kitagawa, T. (2015). A Test for Instrument Validity. *Econometrica*,
#' 83(5), 2043-2063. \doi{10.3982/ECTA11974}
#'
#' Imbens, G. W. and Angrist, J. D. (1994). Identification and Estimation
#' of Local Average Treatment Effects. *Econometrica*, 62(2), 467-475.
#' \doi{10.2307/2951620}
#'
#' @family iv_tests
#' @seealso [iv_mw()] for the conditional version with covariates,
#'   [iv_testjfe()] for the judge-design test, and [iv_check()] for a
#'   one-shot wrapper that runs all applicable tests.
#'
#' @examples
#' \donttest{
#' # Valid IV: compliers exist, no violations
#' set.seed(1)
#' n <- 500
#' z <- sample(0:1, n, replace = TRUE)
#' d <- rbinom(n, 1, 0.3 + 0.4 * z)
#' y <- rnorm(n, mean = d)
#' iv_kitagawa(y, d, z, n_boot = 200, parallel = FALSE)
#' }
#'
#' @export
iv_kitagawa <- function(object, ...) {
  UseMethod("iv_kitagawa")
}

#' @rdname iv_kitagawa
#' @export
iv_kitagawa.default <- function(object, d, z, n_boot = 1000, alpha = 0.05,
                                weighting = c("variance", "unweighted"),
                                weights = NULL, parallel = TRUE, ...) {
  weighting <- match.arg(weighting)
  y <- object
  validate_numeric(y, "y")
  d_num <- validate_binary(d, "d")
  z_num <- validate_discrete(z, "z")
  check_lengths(y, d_num, z_num)

  if (!is.null(weights)) {
    cli::cli_warn(
      "The {.arg weights} argument is not yet implemented and is ignored in v0.1.0."
    )
  }

  core <- kitagawa_core_test(y, d_num, z_num, n_boot, parallel,
                             weighting = weighting)

  structure(
    list(
      test = "Kitagawa (2015)",
      statistic = core$statistic,
      p_value = core$p_value,
      alpha = alpha,
      n_boot = n_boot,
      boot_stats = core$boot_stats,
      binding = core$binding,
      weighting = weighting,
      n = core$n,
      call = sys.call()
    ),
    class = "iv_test"
  )
}

#' @rdname iv_kitagawa
#' @export
iv_kitagawa.fixest <- function(object, n_boot = 1000, alpha = 0.05,
                               weighting = c("variance", "unweighted"),
                               weights = NULL, parallel = TRUE, ...) {
  weighting <- match.arg(weighting)
  yz <- extract_iv_data(object)
  iv_kitagawa.default(
    object = yz$y, d = yz$d, z = yz$z,
    n_boot = n_boot, alpha = alpha, weighting = weighting,
    weights = weights, parallel = parallel, ...
  )
}

#' @rdname iv_kitagawa
#' @export
iv_kitagawa.ivreg <- function(object, n_boot = 1000, alpha = 0.05,
                              weighting = c("variance", "unweighted"),
                              weights = NULL, parallel = TRUE, ...) {
  weighting <- match.arg(weighting)
  yz <- extract_iv_data(object)
  iv_kitagawa.default(
    object = yz$y, d = yz$d, z = yz$z,
    n_boot = n_boot, alpha = alpha, weighting = weighting,
    weights = weights, parallel = parallel, ...
  )
}
