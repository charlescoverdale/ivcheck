#' Frandsen-Lefgren-Leslie (2023) test for instrument validity in
#' judge-fixed-effects designs
#'
#' Jointly tests the local exclusion and monotonicity assumptions when
#' treatment is binary and the instruments are a set of mutually
#' exclusive dummy variables (the leniency-of-assigned-judge design).
#' Rejection is evidence that at least one of exclusion or monotonicity
#' fails.
#'
#' @inheritParams iv_kitagawa
#' @param z Factor, integer, or matrix of mutually exclusive dummy
#'   variables identifying the judge (or other random-assignment unit).
#' @param x Optional covariates to condition on.
#'
#' @return An object of class `iv_test`; see [iv_kitagawa] for element
#'   descriptions.
#'
#' @references
#' Frandsen, B. R., Lefgren, L. J., and Leslie, E. C. (2023). Judging
#' Judge Fixed Effects. *American Economic Review*, 113(1), 253-277.
#' \doi{10.1257/aer.20201860}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 1000
#' judge <- sample.int(20, n, replace = TRUE)
#' d <- rbinom(n, 1, 0.3 + 0.02 * judge)
#' y <- rnorm(n, mean = d)
#' iv_testjfe(y, d, judge, n_boot = 200)
#' }
#'
#' @export
iv_testjfe <- function(object, ...) {
  UseMethod("iv_testjfe")
}

#' @rdname iv_testjfe
#' @export
iv_testjfe.default <- function(object, d, z, x = NULL, n_boot = 1000,
                               alpha = 0.05, weights = NULL,
                               parallel = TRUE, ...) {
  y <- object
  validate_numeric(y, "y")
  d_num <- validate_binary(d, "d")
  z_num <- validate_discrete(z, "z", max_levels = 500L)
  n <- check_lengths(y, d_num, z_num)

  # TODO v0.1.0: port Frandsen-Lefgren-Leslie restricted-LS test from
  # the Stata `testjfe` module.
  boot_stats <- rep(NA_real_, n_boot)

  structure(
    list(
      test = "Frandsen-Lefgren-Leslie (2023)",
      statistic = NA_real_,
      p_value = NA_real_,
      alpha = alpha,
      n_boot = n_boot,
      boot_stats = boot_stats,
      binding = NULL,
      n_judges = length(unique(z_num)),
      n = n,
      call = sys.call()
    ),
    class = "iv_test"
  )
}

#' @rdname iv_testjfe
#' @export
iv_testjfe.fixest <- function(object, x = NULL, n_boot = 1000,
                              alpha = 0.05, weights = NULL,
                              parallel = TRUE, ...) {
  yz <- extract_iv_data(object)
  iv_testjfe.default(
    object = yz$y, d = yz$d, z = yz$z, x = x %||% yz$x,
    n_boot = n_boot, alpha = alpha,
    weights = weights, parallel = parallel, ...
  )
}

#' @rdname iv_testjfe
#' @export
iv_testjfe.ivreg <- function(object, x = NULL, n_boot = 1000,
                             alpha = 0.05, weights = NULL,
                             parallel = TRUE, ...) {
  yz <- extract_iv_data(object)
  iv_testjfe.default(
    object = yz$y, d = yz$d, z = yz$z, x = x %||% yz$x,
    n_boot = n_boot, alpha = alpha,
    weights = weights, parallel = parallel, ...
  )
}
