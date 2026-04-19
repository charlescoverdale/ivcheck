#' Mourifie-Wan (2017) test for instrument validity
#'
#' Reformulates the testable implications of Kitagawa (2015) as a set of
#' conditional moment inequalities. Without covariates `x`, `iv_mw` tests
#' the same inequalities as `iv_kitagawa` and uses the same multiplier-
#' bootstrap critical values (the two tests are asymptotically equivalent
#' at the boundary of the null). With covariates, `iv_mw` stratifies the
#' sample by quantile bins of `x`, runs the inequality test within each
#' bin, and takes the max statistic across bins with a jointly-drawn
#' multiplier bootstrap.
#'
#' @inheritParams iv_kitagawa
#' @param x Optional numeric vector, matrix, or data frame of covariates.
#'   If supplied, the test is stratified by quantile bins of the first
#'   numeric column of `x`. If `NULL`, the test reduces to the
#'   unconditional Mourifie-Wan test.
#' @param n_bins Number of quantile bins of `x`. Default 5. Ignored when
#'   `x` is `NULL`.
#' @param grid Reserved. Custom outcome-evaluation grid for the CLR-style
#'   inference in v0.2.0.
#'
#' @return An object of class `iv_test`; see [iv_kitagawa] for element
#'   descriptions. Additional elements:
#'   \item{conditional}{Logical, whether covariates were supplied.}
#'   \item{n_bins}{Number of x-bins used when `conditional` is `TRUE`.}
#'
#' @details
#' This is a first-pass implementation of Mourifie-Wan (2017). The
#' unconditional path is exact. The conditional path uses a simple bin-
#' stratified version that preserves the spirit of the CLR (2013)
#' intersection-bounds framework without yet implementing the full
#' Chernozhukov-Lee-Rosen inference. A v0.2.0 release will add the full
#' CLR machinery; users needing that today should consult the Stata
#' implementation accompanying the CLR (2015) paper.
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
iv_mw.default <- function(object, d, z, x = NULL, n_bins = 5L, grid = NULL,
                          n_boot = 1000, alpha = 0.05,
                          weights = NULL, parallel = TRUE, ...) {
  y <- object
  validate_numeric(y, "y")
  d_num <- validate_binary(d, "d")
  z_num <- validate_discrete(z, "z")
  n <- check_lengths(y, d_num, z_num)

  if (!is.null(weights)) {
    cli::cli_warn(
      "The {.arg weights} argument is not yet implemented and is ignored in v0.1.0."
    )
  }
  if (!is.null(grid)) {
    cli::cli_warn(
      "The {.arg grid} argument is reserved for the v0.2.0 CLR-style path and is ignored."
    )
  }

  if (is.null(x)) {
    core <- kitagawa_core_test(y, d_num, z_num, n_boot, alpha, parallel)
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
        n_bins = NA_integer_,
        n = core$n,
        call = sys.call()
      ),
      class = "iv_test"
    ))
  }

  # Conditional path: stratify by quantile bins of the first x-column.
  x_mat <- if (is.matrix(x) || is.data.frame(x)) as.matrix(x) else cbind(x)
  x_stratify <- x_mat[, 1]
  validate_numeric(x_stratify, "x")
  if (length(x_stratify) != n) {
    cli::cli_abort("{.arg x} must have the same number of rows as {.arg y}.")
  }

  bin_edges <- stats::quantile(
    x_stratify, probs = seq(0, 1, length.out = n_bins + 1L), na.rm = TRUE
  )
  bin_edges[1] <- bin_edges[1] - .Machine$double.eps * abs(bin_edges[1] + 1)
  bin <- as.integer(cut(x_stratify, breaks = bin_edges, include.lowest = TRUE))

  per_bin <- vector("list", n_bins)
  for (b in seq_len(n_bins)) {
    idx <- which(bin == b)
    if (length(idx) < 20L) {
      per_bin[[b]] <- NULL
      next
    }
    # Respect the user's parallel choice but serialise inside the bin
    # loop to keep total cores bounded.
    per_bin[[b]] <- kitagawa_core_test(
      y[idx], d_num[idx], z_num[idx],
      n_boot = n_boot, alpha = alpha, parallel = FALSE
    )
  }
  valid <- which(!vapply(per_bin, is.null, logical(1)))
  if (length(valid) == 0L) {
    cli::cli_abort("No x-bin has sufficient observations (>= 20) to run the test.")
  }

  stats_by_bin <- vapply(per_bin[valid],
                         function(r) r$statistic, numeric(1))
  boot_mat <- do.call(cbind, lapply(per_bin[valid], function(r) r$boot_stats))
  joint_obs <- max(stats_by_bin)
  joint_boot <- apply(boot_mat, 1, max)
  p_value <- mean(joint_boot >= joint_obs)

  binding_bin <- valid[which.max(stats_by_bin)]

  structure(
    list(
      test = "Mourifie-Wan (2017)",
      statistic = joint_obs,
      p_value = p_value,
      alpha = alpha,
      n_boot = n_boot,
      boot_stats = joint_boot,
      binding = c(per_bin[[binding_bin]]$binding, list(x_bin = binding_bin)),
      conditional = TRUE,
      n_bins = length(valid),
      n = n,
      call = sys.call()
    ),
    class = "iv_test"
  )
}

#' @rdname iv_mw
#' @export
iv_mw.fixest <- function(object, x = NULL, n_bins = 5L, grid = NULL,
                         n_boot = 1000, alpha = 0.05, weights = NULL,
                         parallel = TRUE, ...) {
  yz <- extract_iv_data(object)
  iv_mw.default(
    object = yz$y, d = yz$d, z = yz$z,
    x = x %||% yz$x, n_bins = n_bins, grid = grid,
    n_boot = n_boot, alpha = alpha,
    weights = weights, parallel = parallel, ...
  )
}

#' @rdname iv_mw
#' @export
iv_mw.ivreg <- function(object, x = NULL, n_bins = 5L, grid = NULL,
                        n_boot = 1000, alpha = 0.05, weights = NULL,
                        parallel = TRUE, ...) {
  yz <- extract_iv_data(object)
  iv_mw.default(
    object = yz$y, d = yz$d, z = yz$z,
    x = x %||% yz$x, n_bins = n_bins, grid = grid,
    n_boot = n_boot, alpha = alpha,
    weights = weights, parallel = parallel, ...
  )
}
