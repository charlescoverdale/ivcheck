#' Mourifie-Wan (2017) test for instrument validity
#'
#' Reformulates the testable implications of Kitagawa (2015) as a set of
#' conditional moment inequalities and tests them in the intersection-
#' bounds framework of Chernozhukov, Lee, and Rosen (2013). Without
#' covariates `x`, `iv_mw` tests the same inequalities as `iv_kitagawa`
#' and reduces exactly to the variance-weighted Kitagawa test. With
#' covariates, `iv_mw` estimates the conditional CDFs
#' `F(y, d | X = x, Z = z)` nonparametrically via series regression,
#' computes plug-in heteroscedasticity-robust standard errors, and
#' takes the sup over `(y, x)` of the variance-weighted positive-part
#' violation. Critical values come from a multiplier bootstrap with
#' adaptive moment selection in the style of Andrews and Soares (2010).
#'
#' @inheritParams iv_kitagawa
#' @param x Optional numeric vector, matrix, or data frame of covariates.
#'   If supplied, the test is conditional on the first numeric column
#'   of `x`. If `NULL`, the test reduces to the unconditional
#'   Mourifie-Wan test.
#' @param basis_order Polynomial order of the series-regression basis
#'   used to estimate `F(y, d | X, Z)`. Default `3L` (cubic). Set to
#'   `"auto"` to select the basis order by 5-fold cross-validation over
#'   the candidates 2, 3, 4, 5 with squared-error loss on the indicator
#'   regression. When `"auto"` is used, the bootstrap becomes
#'   post-selection-valid: the test statistic is compared to the
#'   maximum of the bootstrap statistics across the candidate orders,
#'   which controls size at the nominal level against any selection
#'   rule but is mildly conservative relative to a fixed-order test.
#'   Runtime with `"auto"` is approximately four times the fixed-order
#'   path.
#' @param x_grid_size Number of quantile points of `x` at which to
#'   evaluate the conditional CDFs. Default 20.
#' @param y_grid_size Number of quantile points of `y` at which to
#'   evaluate the inequalities. Default 50.
#' @param adaptive Logical. If `TRUE` (default), the bootstrap uses the
#'   adaptive moment selection of Andrews-Soares (2010) with tuning
#'   parameter `kappa_n = sqrt(log(log(n)))`. If `FALSE`, uses the
#'   plug-in least-favourable critical value (conservative).
#' @param grid Deprecated. Ignored; use `y_grid_size` and
#'   `x_grid_size` instead.
#'
#' @return An object of class `iv_test`; see [iv_kitagawa] for element
#'   descriptions. Additional elements:
#'   \item{conditional}{Logical, whether covariates were supplied.}
#'   \item{kappa_n}{Andrews-Soares tuning parameter used
#'     (`NA` if not applicable).}
#'
#' @details
#' The CLR framework targets conditional moment inequalities of the
#' form `E[m(W; theta) | X] <= 0` for all `X`. Applied to Kitagawa's
#' (2015) inequalities, the relevant moments are the positive-part
#' differences of the conditional joint CDFs `F(y, d | X, Z)` for each
#' `(d, z_low, z_high, y, x)` index. `iv_mw` estimates `F(y, d | X, Z)`
#' by series regression of the indicator `1{Y <= y, D = d}` on a
#' polynomial basis of `X` within each `Z` cell. Robust standard errors
#' come from the heteroscedasticity-consistent sandwich of the series
#' regression. Critical values are drawn by multiplier bootstrap: the
#' bootstrap process reuses the plug-in SE denominator and perturbs
#' the residuals by Rademacher weights, projected back through the
#' basis. Adaptive moment selection includes only moments whose
#' observed studentised statistic is within `kappa_n` of the
#' inequality boundary, giving tighter critical values when some
#' inequalities are strictly slack.
#'
#' @references
#' Mourifie, I. and Wan, Y. (2017). Testing Local Average Treatment
#' Effect Assumptions. *Review of Economics and Statistics*, 99(2),
#' 305-313. \doi{10.1162/REST_a_00622}
#'
#' Chernozhukov, V., Lee, S., and Rosen, A. M. (2013). Intersection
#' Bounds: Estimation and Inference. *Econometrica*, 81(2), 667-737.
#' \doi{10.3982/ECTA8718}
#'
#' Imbens, G. W. and Angrist, J. D. (1994). Identification and Estimation
#' of Local Average Treatment Effects. *Econometrica*, 62(2), 467-475.
#' \doi{10.2307/2951620}
#'
#' @family iv_tests
#' @seealso [iv_kitagawa()] for the unconditional case,
#'   [iv_testjfe()] for the judge-design test, and [iv_check()] for a
#'   one-shot wrapper that runs all applicable tests.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 500
#' z <- sample(0:1, n, replace = TRUE)
#' d <- rbinom(n, 1, 0.3 + 0.4 * z)
#' y <- rnorm(n, mean = d)
#' iv_mw(y, d, z, n_boot = 200, parallel = FALSE)
#' }
#'
#' @export
iv_mw <- function(object, ...) {
  UseMethod("iv_mw")
}

#' @rdname iv_mw
#' @export
iv_mw.default <- function(object, d, z, x = NULL,
                          basis_order = 3L,
                          x_grid_size = 20L,
                          y_grid_size = 50L,
                          adaptive = TRUE,
                          grid = NULL,
                          n_boot = 1000, alpha = 0.05,
                          weights = NULL, parallel = TRUE, ...) {
  y <- object
  validate_numeric(y, "y")
  d_num <- validate_binary(d, "d")
  z_num <- validate_discrete(z, "z")
  n <- as.integer(check_lengths(y, d_num, z_num))

  if (!is.null(grid)) {
    cli::cli_warn(
      "The {.arg grid} argument is deprecated; use {.arg x_grid_size} and {.arg y_grid_size}."
    )
  }
  if (!is.null(weights) && !is.null(x)) {
    cli::cli_warn(
      "The {.arg weights} argument is not yet implemented for the conditional (x) path; it will be applied only to the no-x delegation or ignored under conditional inference."
    )
  }

  if (is.null(x)) {
    core <- kitagawa_core_test(y, d_num, z_num, n_boot, parallel,
                               weighting = "variance", weights = weights)
    return(structure(
      list(
        test = "Mourifie-Wan (2017)",
        statistic = core$statistic,
        p_value = core$p_value,
        alpha = alpha,
        n_boot = n_boot,
        boot_stats = core$boot_stats,
        binding = core$binding,
        conditional = FALSE,
        kappa_n = NA_real_,
        n = core$n,
        call = sys.call()
      ),
      class = "iv_test"
    ))
  }

  x_mat <- if (is.matrix(x) || is.data.frame(x)) as.matrix(x) else cbind(x)
  validate_numeric(as.vector(x_mat), "x")
  if (nrow(x_mat) != n) {
    cli::cli_abort("{.arg x} must have the same number of rows as {.arg y}.")
  }
  if (ncol(x_mat) > 1L) {
    cli::cli_abort(c(
      "{.fn iv_mw} v0.1.0 supports a single conditioning covariate.",
      i = "Multivariate conditioning requires a tensor-product basis (planned for v0.2.0).",
      i = "Workaround: condition on the covariate most plausibly driving heterogeneity in compliance, or reduce {.arg x} to a single pre-computed index."
    ))
  }

  use_auto <- identical(basis_order, "auto")
  if (use_auto) {
    core <- mw_clr_test_posi(y, d_num, z_num, x_mat, n_boot, parallel,
                             candidates = 2:5,
                             x_grid_size = x_grid_size,
                             y_grid_size = y_grid_size,
                             adaptive = adaptive)
  } else {
    if (!is.numeric(basis_order) || length(basis_order) != 1L ||
        basis_order < 1L || basis_order != as.integer(basis_order)) {
      cli::cli_abort(
        "{.arg basis_order} must be a positive integer or the string {.val auto}."
      )
    }
    core <- mw_clr_test(y, d_num, z_num, x_mat, n_boot, parallel,
                        basis_order = as.integer(basis_order),
                        x_grid_size = x_grid_size,
                        y_grid_size = y_grid_size,
                        adaptive = adaptive)
  }

  structure(
    list(
      test = "Mourifie-Wan (2017)",
      statistic = core$statistic,
      p_value = core$p_value,
      alpha = alpha,
      n_boot = n_boot,
      boot_stats = core$boot_stats,
      binding = core$binding,
      conditional = TRUE,
      kappa_n = core$kappa_n,
      basis_order = if (use_auto) core$basis_order else as.integer(basis_order),
      basis_order_auto = use_auto,
      basis_order_cv = if (use_auto) core$basis_order_cv else NULL,
      post_selection_valid = use_auto,
      n = core$n,
      call = sys.call()
    ),
    class = "iv_test"
  )
}

#' @rdname iv_mw
#' @export
iv_mw.fixest <- function(object, x = NULL,
                         basis_order = 3L, x_grid_size = 20L,
                         y_grid_size = 50L, adaptive = TRUE,
                         grid = NULL, n_boot = 1000, alpha = 0.05,
                         weights = NULL, parallel = TRUE, ...) {
  yz <- extract_iv_data(object)
  iv_mw.default(
    object = yz$y, d = yz$d, z = yz$z,
    x = x %||% yz$x,
    basis_order = basis_order, x_grid_size = x_grid_size,
    y_grid_size = y_grid_size, adaptive = adaptive, grid = grid,
    n_boot = n_boot, alpha = alpha,
    weights = weights, parallel = parallel, ...
  )
}

#' @rdname iv_mw
#' @export
iv_mw.ivreg <- function(object, x = NULL,
                        basis_order = 3L, x_grid_size = 20L,
                        y_grid_size = 50L, adaptive = TRUE,
                        grid = NULL, n_boot = 1000, alpha = 0.05,
                        weights = NULL, parallel = TRUE, ...) {
  yz <- extract_iv_data(object)
  iv_mw.default(
    object = yz$y, d = yz$d, z = yz$z,
    x = x %||% yz$x,
    basis_order = basis_order, x_grid_size = x_grid_size,
    y_grid_size = y_grid_size, adaptive = adaptive, grid = grid,
    n_boot = n_boot, alpha = alpha,
    weights = weights, parallel = parallel, ...
  )
}
