#' Kitagawa (2015) test for instrument validity
#'
#' Tests the joint implication of the local exclusion restriction and the
#' local monotonicity condition in a binary-treatment, discrete-instrument
#' setting. The null hypothesis is that the instrument is valid. Under the
#' null, the conditional densities of the outcome given treatment status,
#' evaluated at each level of the instrument, must satisfy a set of
#' stochastic dominance inequalities. Rejection is evidence that at least
#' one of exclusion or monotonicity fails.
#'
#' @param object For the default method: a numeric outcome vector.
#'   For the `fixest` and `ivreg` methods: a fitted instrumental variable
#'   model from [fixest::feols] or `ivreg::ivreg()`.
#' @param d Binary 0/1 treatment vector (default method only).
#' @param z Discrete instrument (numeric or factor, default method only).
#' @param n_boot Number of multiplier-bootstrap replications. Default 1000.
#' @param alpha Significance level for the returned verdict. Default 0.05.
#' @param weights Optional survey weights. If `NULL`, unweighted.
#' @param parallel Logical. Run bootstrap replications in parallel on
#'   POSIX systems via [parallel::mclapply]. Default `TRUE`.
#' @param ... Further arguments passed to methods.
#'
#' @return An object of class `iv_test` with elements:
#'   \item{test}{Character, always `"Kitagawa (2015)"`.}
#'   \item{statistic}{Numeric test statistic (variance-weighted
#'     Kolmogorov-Smirnov).}
#'   \item{p_value}{Bootstrap p-value.}
#'   \item{alpha}{Supplied significance level.}
#'   \item{n_boot}{Number of bootstrap replications used.}
#'   \item{boot_stats}{Numeric vector of bootstrap test statistics.}
#'   \item{binding}{The (z, z') pair giving the binding inequality.}
#'   \item{n}{Sample size.}
#'   \item{call}{Matched call.}
#'
#' @references
#' Kitagawa, T. (2015). A Test for Instrument Validity. *Econometrica*,
#' 83(5), 2043-2063. \doi{10.3982/ECTA11974}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 500
#' z <- sample(0:1, n, replace = TRUE)
#' d <- rbinom(n, 1, 0.3 + 0.4 * z)
#' y <- rnorm(n, mean = d)
#' iv_kitagawa(y, d, z, n_boot = 200)
#' }
#'
#' @export
iv_kitagawa <- function(object, ...) {
  UseMethod("iv_kitagawa")
}

#' @rdname iv_kitagawa
#' @export
iv_kitagawa.default <- function(object, d, z, n_boot = 1000, alpha = 0.05,
                                weights = NULL, parallel = TRUE, ...) {
  y <- object
  validate_numeric(y, "y")
  d_num <- validate_binary(d, "d")
  z_num <- validate_discrete(z, "z")
  n <- check_lengths(y, d_num, z_num)

  # TODO v0.1.0: implement Kitagawa (2015) eq. 3.1 KS statistic and
  # multiplier bootstrap per section 3.2.
  boot_stats <- rep(NA_real_, n_boot)

  structure(
    list(
      test = "Kitagawa (2015)",
      statistic = NA_real_,
      p_value = NA_real_,
      alpha = alpha,
      n_boot = n_boot,
      boot_stats = boot_stats,
      binding = NULL,
      n = n,
      call = sys.call()
    ),
    class = "iv_test"
  )
}

#' @rdname iv_kitagawa
#' @export
iv_kitagawa.fixest <- function(object, n_boot = 1000, alpha = 0.05,
                               weights = NULL, parallel = TRUE, ...) {
  yz <- extract_iv_data(object)
  iv_kitagawa.default(
    object = yz$y, d = yz$d, z = yz$z,
    n_boot = n_boot, alpha = alpha,
    weights = weights, parallel = parallel, ...
  )
}

#' @rdname iv_kitagawa
#' @export
iv_kitagawa.ivreg <- function(object, n_boot = 1000, alpha = 0.05,
                              weights = NULL, parallel = TRUE, ...) {
  yz <- extract_iv_data(object)
  iv_kitagawa.default(
    object = yz$y, d = yz$d, z = yz$z,
    n_boot = n_boot, alpha = alpha,
    weights = weights, parallel = parallel, ...
  )
}
