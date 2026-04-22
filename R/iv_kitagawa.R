#' Kitagawa (2015) / Sun (2023) test for instrument validity
#'
#' Tests the joint implication of the local exclusion restriction and the
#' local monotonicity condition in a discrete-instrument setting.
#' Supports binary treatment (Kitagawa 2015), ordered multivalued
#' treatment (Sun 2023 section 3), and unordered multivalued treatment
#' (Sun 2023 section 3.3) under a user-supplied monotonicity set.
#' The null is that the instrument is valid. Under the null, the
#' conditional joint distribution of `(Y, D | Z)` must satisfy
#' stochastic dominance inequalities on cumulative-tail events.
#' Rejection is evidence that at least one of exclusion or monotonicity
#' fails.
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
#' @param se_floor Trimming constant `\xi` for the plug-in standard-
#'   error denominator in the variance-weighted form. Default 0.15.
#'   Kitagawa (2015) section 4 informally recommends `\xi \in [0.05,
#'   0.1]` for balanced-Z designs. Monte Carlo at skewed Z-cell
#'   distributions with weak first stages suggests a slightly larger
#'   floor (0.15) keeps empirical size near nominal 5% without
#'   measurable power loss in the designs tested. Users reproducing
#'   Kitagawa's published examples may set `se_floor = 0.1` to match.
#' @param weights Optional survey weights. A non-negative numeric vector
#'   of length equal to the sample size. Scaled internally so the mean
#'   weight is 1.0 (preserving effective sample-size interpretation).
#'   Applied to the empirical CDFs, the bootstrap multiplier process,
#'   and the variance-weighted standard errors.
#' @param parallel Logical. Run bootstrap replications in parallel on
#'   POSIX systems via [parallel::mclapply]. Default `TRUE`.
#' @param treatment_order Either `"ordered"` (default) or `"unordered"`.
#'   Binary `D` is handled identically under both. For multivalued `D`,
#'   `"ordered"` uses cumulative-tail inequalities
#'   `P(Y <= y, D <= ell | Z)` and `P(Y <= y, D >= ell | Z)` across all
#'   pairs of instrument values, a stronger family of implications than
#'   Sun (2023) equation 10's `d_min`-and-`d_max` subset (but still
#'   valid under Sun's Assumption 2.2). `"unordered"` requires a
#'   user-specified `monotonicity_set` naming the (level, z_from, z_to)
#'   triples for which `1{D_{z_to} = d} <= 1{D_{z_from} = d}` is
#'   assumed almost surely (Sun 2023 Assumption 2.4(iii)).
#' @param monotonicity_set A `data.frame` with columns `d`, `z_from`,
#'   `z_to` listing the triples that pin down the direction of the
#'   monotonicity restriction for `treatment_order = "unordered"`.
#'   Ignored when `treatment_order = "ordered"`.
#' @param multiplier Choice of bootstrap multiplier: `"rademacher"`
#'   (default; +/-1 two-point), `"gaussian"` (standard normal), or
#'   `"mammen"` (Mammen 1993 asymmetric two-point).
#' @param ... Further arguments passed to methods.
#'
#' @return An object of class `iv_test` with elements:
#'   \item{test}{`"Kitagawa (2015)"` for binary treatment;
#'     `"Sun (2023)"` for multivalued ordered treatment.}
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
#' Kitagawa (2015) equation 2.1 defines the statistic as the max over
#' instrument-level pairs `(z_low, z_high)`, treatment status
#' `d in {0, 1}`, and intervals `[y, y']` with `y <= y'`, of the
#' positive-part interval-probability difference normalised by the
#' binomial-mixture plug-in standard error:
#' `T_n = sqrt(n_low * n_high / (n_low + n_high))`
#' `* max [P([y, y'], d | z_low) - P([y, y'], d | z_high)]^+ / sigma_hat`.
#' (The denominator is the pair total, not the full sample size.) The sign flips for
#' `d = 0`. Instrument levels are pre-ordered by first-stage
#' `E_hat[D | Z]` so the inequalities are one-sided and `T_n >= 0`.
#' The implementation evaluates the sup on a quantile grid of observed
#' outcomes (default 50 points); this is equivalent to evaluation at
#' every sample-point pair under Kitagawa's Theorem 2.1. Critical
#' values come from a multiplier bootstrap (section 3.2) of the pooled
#' empirical distribution; bootstrap statistics reuse the data-derived
#' standard-error denominator.
#'
#' @references
#' Kitagawa, T. (2015). A Test for Instrument Validity. *Econometrica*,
#' 83(5), 2043-2063. \doi{10.3982/ECTA11974}
#'
#' Sun, Z. (2023). Instrument validity for heterogeneous causal effects.
#' *Journal of Econometrics*, 237(2), 105523. \doi{10.1016/j.jeconom.2023.105523}
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
                                weights = NULL, parallel = TRUE,
                                se_floor = 0.15,
                                treatment_order = c("ordered", "unordered"),
                                monotonicity_set = NULL,
                                multiplier = c("rademacher", "gaussian",
                                               "mammen"),
                                ...) {
  weighting <- match.arg(weighting)
  treatment_order <- match.arg(treatment_order)
  multiplier <- match.arg(multiplier)
  y <- object
  validate_numeric(y, "y")
  d_info <- validate_treatment_discrete(d, "d", max_levels = 20L)
  d_num <- d_info$d_num
  z_num <- validate_discrete(z, "z")
  check_lengths(y, d_num, z_num)

  # Re-map user-supplied monotonicity_set d values to the internal 0..k-1
  # coding used by kitagawa_core_test. Validate structure here so the
  # error arrives before the internal remap runs.
  ms_internal <- monotonicity_set
  if (!is.null(ms_internal)) {
    required <- c("d", "z_from", "z_to")
    if (!all(required %in% names(ms_internal))) {
      cli::cli_abort(
        "{.arg monotonicity_set} must have columns {.val {required}}."
      )
    }
    if (!is.null(d_info$original_levels)) {
      original_to_internal <- stats::setNames(
        seq_along(d_info$original_levels) - 1L,
        d_info$original_levels
      )
      ms_internal$d <- original_to_internal[as.character(ms_internal$d)]
    }
  }

  core <- kitagawa_core_test(y, d_num, z_num, n_boot, parallel,
                             weighting = weighting, weights = weights,
                             se_floor = se_floor,
                             d_labels = d_info$original_levels,
                             treatment_order = treatment_order,
                             monotonicity_set = ms_internal,
                             multiplier = multiplier)

  test_name <- if (core$multivalued) "Sun (2023)" else "Kitagawa (2015)"
  structure(
    list(
      test = test_name,
      statistic = core$statistic,
      p_value = core$p_value,
      alpha = alpha,
      n_boot = n_boot,
      boot_stats = core$boot_stats,
      binding = core$binding,
      weighting = weighting,
      multiplier = multiplier,
      treatment_order = treatment_order,
      n_treatment_levels = core$n_treatment_levels,
      multivalued = core$multivalued,
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
                               weights = NULL, parallel = TRUE,
                               treatment_order = c("ordered", "unordered"),
                               monotonicity_set = NULL,
                               multiplier = c("rademacher", "gaussian",
                                              "mammen"),
                               ...) {
  weighting <- match.arg(weighting)
  treatment_order <- match.arg(treatment_order)
  multiplier <- match.arg(multiplier)
  yz <- extract_iv_data(object)
  iv_kitagawa.default(
    object = yz$y, d = yz$d, z = yz$z,
    n_boot = n_boot, alpha = alpha, weighting = weighting,
    weights = weights, parallel = parallel,
    treatment_order = treatment_order,
    monotonicity_set = monotonicity_set,
    multiplier = multiplier, ...
  )
}

#' @rdname iv_kitagawa
#' @export
iv_kitagawa.ivreg <- function(object, n_boot = 1000, alpha = 0.05,
                              weighting = c("variance", "unweighted"),
                              weights = NULL, parallel = TRUE,
                              treatment_order = c("ordered", "unordered"),
                              monotonicity_set = NULL,
                              multiplier = c("rademacher", "gaussian",
                                             "mammen"),
                              ...) {
  weighting <- match.arg(weighting)
  treatment_order <- match.arg(treatment_order)
  multiplier <- match.arg(multiplier)
  yz <- extract_iv_data(object)
  iv_kitagawa.default(
    object = yz$y, d = yz$d, z = yz$z,
    n_boot = n_boot, alpha = alpha, weighting = weighting,
    weights = weights, parallel = parallel,
    treatment_order = treatment_order,
    monotonicity_set = monotonicity_set,
    multiplier = multiplier, ...
  )
}
