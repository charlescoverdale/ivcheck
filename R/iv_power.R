#' Monte Carlo power curve for IV-validity tests
#'
#' Simulates data under a user-specified deviation from validity and
#' estimates the rejection probability of the chosen test at each
#' deviation size. Useful for sample-size planning and for benchmarking
#' different tests on the same design.
#'
#' The deviation is parameterised as the size of a **D-specific direct
#' effect of the instrument on the outcome** (a clean exclusion
#' violation that the Kitagawa and Mourifie-Wan tests are designed to
#' detect). Specifically, the simulated outcome is
#' `Y = mu_hat[D + 1] + delta * sigma_hat * D * (Z - Z_low) + noise`,
#' so `delta = 0` corresponds to the null and larger values produce
#' larger violations of the testable inequality for the d = 1 cells.
#' The simulator preserves the observed sample size, first-stage
#' propensities, and outcome scale.
#'
#' @param y,d,z Observed data used to anchor the DGP (sample size, cell
#'   counts, empirical first-stage).
#' @param method Which test to benchmark. One of `"kitagawa"`, `"mw"`,
#'   or `"testjfe"`.
#' @param alpha Significance level.
#' @param n_sims Number of Monte Carlo simulations per deviation.
#' @param delta_grid Numeric vector of deviation sizes to evaluate.
#'   If `NULL`, defaults to `seq(0, 0.3, by = 0.05)`.
#' @param n_boot Number of bootstrap replications per simulation (for
#'   tests that use bootstrap). Default 200, which trades some Monte
#'   Carlo noise for tractable runtime.
#' @param parallel Logical. Run simulations in parallel on POSIX systems
#'   via [parallel::mclapply]. Default `TRUE`.
#' @param ... Further arguments passed to the underlying test.
#'
#' @return A data frame with columns `delta` (deviation size) and
#'   `power` (estimated rejection probability at level `alpha`).
#'
#' @examples
#' \dontrun{
#' # Headline power curve for a small-N design
#' set.seed(1)
#' n <- 300
#' z <- sample(0:1, n, replace = TRUE)
#' d <- rbinom(n, 1, 0.3 + 0.4 * z)
#' y <- rnorm(n, mean = d)
#' iv_power(y, d, z, method = "kitagawa", n_sims = 50, n_boot = 100)
#' }
#'
#' @export
iv_power <- function(y, d, z, method = c("kitagawa", "mw", "testjfe"),
                     alpha = 0.05, n_sims = 500, delta_grid = NULL,
                     n_boot = 200, parallel = TRUE, ...) {
  method <- match.arg(method)
  validate_numeric(y, "y")
  d_num <- validate_binary(d, "d")
  z_num <- validate_discrete(z, "z", max_levels = 500L)
  n <- check_lengths(y, d_num, z_num)
  if (is.null(delta_grid)) {
    delta_grid <- seq(0, 0.3, by = 0.05)
  }

  # DGP anchored to observed data:
  #   - Sample size n kept fixed.
  #   - Principal-strata shares (AT, NT, C) estimated from empirical
  #     first-stage moments.
  #   - Outcome scale sigma estimated from within-(D, Z) variance of y.
  #   - Violation: monotonicity failure via a defier share `delta`,
  #     siphoned from compliers.
  z_vals <- sort(unique(z_num))
  K_z <- length(z_vals)
  p_by_z <- vapply(z_vals, function(zv) mean(d_num[z_num == zv]), numeric(1))
  mu_by_d <- vapply(c(0, 1), function(dv) mean(y[d_num == dv], na.rm = TRUE),
                    numeric(1))
  cells <- expand.grid(d = c(0, 1), z = z_vals)
  sds <- vapply(seq_len(nrow(cells)), function(i) {
    idx <- which(d_num == cells$d[i] & z_num == cells$z[i])
    if (length(idx) < 2L) return(NA_real_)
    stats::sd(y[idx])
  }, numeric(1))
  sigma_hat <- mean(sds, na.rm = TRUE)
  if (!is.finite(sigma_hat) || sigma_hat <= 0) sigma_hat <- stats::sd(y)

  test_fn <- switch(method,
    kitagawa = function(ys, ds, zs) iv_kitagawa(ys, ds, zs, n_boot = n_boot,
                                                alpha = alpha,
                                                parallel = FALSE, ...)$p_value,
    mw       = function(ys, ds, zs) iv_mw(ys, ds, zs, n_boot = n_boot,
                                          alpha = alpha,
                                          parallel = FALSE, ...)$p_value,
    testjfe  = function(ys, ds, zs) iv_testjfe(ys, ds, zs, n_boot = n_boot,
                                               alpha = alpha,
                                               parallel = FALSE, ...)$p_value
  )

  z_lo <- z_vals[1]

  simulate_one <- function(delta) {
    zs <- sample(z_vals, n, replace = TRUE)
    p_vec <- p_by_z[match(zs, z_vals)]
    ds <- stats::rbinom(n, 1, p_vec)
    # D-specific direct-effect exclusion violation: shifts Y upward for
    # the (D = 1, Z != Z_low) cells, which is exactly the direction the
    # Kitagawa inequality for d = 1 tests.
    ys <- mu_by_d[ds + 1L] +
      delta * sigma_hat * ds * (zs - z_lo) +
      stats::rnorm(n, 0, sigma_hat)
    test_fn(ys, ds, zs)
  }

  run_delta <- function(delta) {
    if (isTRUE(parallel) && .Platform$OS.type == "unix") {
      cran_env <- nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
      cores <- if (cran_env) 2L else min(4L, max(1L, parallel::detectCores() - 1L))
      ps <- unlist(parallel::mclapply(seq_len(n_sims),
                                      function(s) simulate_one(delta),
                                      mc.cores = cores),
                   use.names = FALSE)
    } else {
      ps <- vapply(seq_len(n_sims),
                   function(s) simulate_one(delta),
                   numeric(1))
    }
    mean(ps < alpha, na.rm = TRUE)
  }

  power <- vapply(delta_grid, run_delta, numeric(1))
  data.frame(delta = delta_grid, power = power)
}
