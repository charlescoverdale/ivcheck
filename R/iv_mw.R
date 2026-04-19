#' Mourifie-Wan (2017) test for instrument validity
#'
#' Reformulates the testable implications of Kitagawa (2015) as a set of
#' conditional moment inequalities and tests them in the intersection-
#' bounds framework of Chernozhukov, Lee, and Rosen (2013). Unlike
#' [iv_kitagawa], `iv_mw` accepts covariates `x` and tests conditional
#' on them.
#'
#' @inheritParams iv_kitagawa
#' @param x Optional numeric matrix or data frame of covariates. If
#'   supplied, the test is conditional on `x`.
#' @param grid Optional numeric grid of outcome values at which to
#'   evaluate the inequalities. If `NULL`, a data-adaptive grid is used.
#'
#' @return An object of class `iv_test`; see [iv_kitagawa] for element
#'   descriptions. Additionally:
#'   \item{conditional}{Logical, whether covariates were supplied.}
#'
#' @references
#' Mourifie, I. and Wan, Y. (2017). Testing Local Average Treatment
#' Effect Assumptions. *Review of Economics and Statistics*, 99(2),
#' 305-313. \doi{10.1162/REST_a_00628}
#'
#' Chernozhukov, V., Lee, S., and Rosen, A. M. (2013). Intersection
#' Bounds: Estimation and Inference. *Econometrica*, 81(2), 667-737.
#' \doi{10.3982/ECTA8718}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 500
#' z <- sample(0:1, n, replace = TRUE)
#' d <- rbinom(n, 1, 0.3 + 0.4 * z)
#' y <- rnorm(n, mean = d)
#' iv_mw(y, d, z, n_boot = 200)
#' }
#'
#' @export
iv_mw <- function(object, ...) {
  UseMethod("iv_mw")
}

#' @rdname iv_mw
#' @export
iv_mw.default <- function(object, d, z, x = NULL, grid = NULL,
                          n_boot = 1000, alpha = 0.05,
                          weights = NULL, parallel = TRUE, ...) {
  y <- object
  validate_numeric(y, "y")
  d_num <- validate_binary(d, "d")
  z_num <- validate_discrete(z, "z")
  n <- check_lengths(y, d_num, z_num)

  # TODO v0.1.0: implement MW 2017 conditional inequalities + CLR-style
  # inference.
  boot_stats <- rep(NA_real_, n_boot)

  structure(
    list(
      test = "Mourifie-Wan (2017)",
      statistic = NA_real_,
      p_value = NA_real_,
      alpha = alpha,
      n_boot = n_boot,
      boot_stats = boot_stats,
      binding = NULL,
      conditional = !is.null(x),
      n = n,
      call = sys.call()
    ),
    class = "iv_test"
  )
}

#' @rdname iv_mw
#' @export
iv_mw.fixest <- function(object, x = NULL, grid = NULL, n_boot = 1000,
                         alpha = 0.05, weights = NULL,
                         parallel = TRUE, ...) {
  yz <- extract_iv_data(object)
  iv_mw.default(
    object = yz$y, d = yz$d, z = yz$z, x = x %||% yz$x, grid = grid,
    n_boot = n_boot, alpha = alpha,
    weights = weights, parallel = parallel, ...
  )
}

#' @rdname iv_mw
#' @export
iv_mw.ivreg <- function(object, x = NULL, grid = NULL, n_boot = 1000,
                        alpha = 0.05, weights = NULL,
                        parallel = TRUE, ...) {
  yz <- extract_iv_data(object)
  iv_mw.default(
    object = yz$y, d = yz$d, z = yz$z, x = x %||% yz$x, grid = grid,
    n_boot = n_boot, alpha = alpha,
    weights = weights, parallel = parallel, ...
  )
}
