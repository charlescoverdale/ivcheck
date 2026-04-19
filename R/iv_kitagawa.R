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

  # Order instrument levels ascending by first-stage E[D | Z]. Under
  # monotonicity, the testable inequalities hold in this ordering:
  #   F(y, 1 | z_low) <= F(y, 1 | z_high)
  #   F(y, 0 | z_high) <= F(y, 0 | z_low)
  # Pairwise violations of these inequalities are the statistic.
  raw_levels <- sort(unique(z_num))
  e_d_by_z <- vapply(raw_levels,
                     function(zk) mean(d_num[z_num == zk]),
                     numeric(1))
  ord <- order(e_d_by_z)
  z_levels <- raw_levels[ord]
  K <- length(z_levels)
  y_grid <- sort(unique(y))
  G <- length(y_grid)

  # Indicator matrix leY[i, g] = 1{Y_i <= y_grid[g]}.
  leY <- outer(y, y_grid, "<=")

  # F_arr[g, k, d_idx]: empirical joint CDF F_hat(y_grid[g], d | z_levels[k])
  # where z_levels is sorted by first-stage E[D | Z] ascending.
  F_arr <- array(0, dim = c(G, K, 2))
  n_z <- integer(K)
  idx_by_z <- vector("list", K)

  for (k in seq_len(K)) {
    zk <- z_levels[k]
    idx <- which(z_num == zk)
    idx_by_z[[k]] <- idx
    n_z[k] <- length(idx)
    if (n_z[k] == 0L) next
    d_sub <- d_num[idx]
    leY_sub <- leY[idx, , drop = FALSE]
    F_arr[, k, 2] <- colSums(leY_sub * d_sub) / n_z[k]
    F_arr[, k, 1] <- colSums(leY_sub * (1 - d_sub)) / n_z[k]
  }

  # Observed statistic: positive-part KS across all ordered pairs
  # (k_low, k_high) with k_low < k_high (ascending E[D | Z]).
  # For d = 1: violation is F(y, 1 | low) - F(y, 1 | high) > 0.
  # For d = 0: violation is F(y, 0 | high) - F(y, 0 | low) > 0.
  compute_stat <- function(F_arr_local) {
    best <- 0
    binding_local <- NULL
    if (K < 2L) return(list(stat = 0, binding = NULL))
    for (k_low in seq_len(K - 1L)) {
      for (k_high in (k_low + 1L):K) {
        # d = 1 inequality
        diffs1 <- F_arr_local[, k_low, 2] - F_arr_local[, k_high, 2]
        g1 <- which.max(diffs1)
        if (diffs1[g1] > best) {
          best <- diffs1[g1]
          binding_local <- list(
            z_low = z_levels[k_low],
            z_high = z_levels[k_high],
            d = 1L,
            y = y_grid[g1]
          )
        }
        # d = 0 inequality
        diffs0 <- F_arr_local[, k_high, 1] - F_arr_local[, k_low, 1]
        g0 <- which.max(diffs0)
        if (diffs0[g0] > best) {
          best <- diffs0[g0]
          binding_local <- list(
            z_low = z_levels[k_low],
            z_high = z_levels[k_high],
            d = 0L,
            y = y_grid[g0]
          )
        }
      }
    }
    list(stat = best, binding = binding_local)
  }

  obs <- compute_stat(F_arr)
  T_n <- sqrt(n) * obs$stat

  # Precompute centred indicator matrices for the multiplier bootstrap.
  # For each d in {0, 1} and k in 1..K, define
  #   C[i, g] = 1{Y_i <= y_grid[g], D_i = d} - F_hat(y_grid[g], d | Z_i)
  # only for i in idx_by_z[[k]]. The bootstrap process at (y_grid[g], d, z_k)
  # is then (1 / n_z[k]) * sum over i in idx_by_z[[k]] of W_i * C[i, g].
  # We keep a (n x G) centred matrix per d for fast matrix-vector ops.
  centred <- list()
  for (d_idx in 1:2) {
    d_val <- d_idx - 1L
    mat <- matrix(0, nrow = n, ncol = G)
    raw <- leY * (d_num == d_val) # n x G, 0/1 indicators
    for (k in seq_len(K)) {
      idx <- idx_by_z[[k]]
      if (length(idx) == 0L) next
      # Subtract F_hat row for that z: broadcast F_arr[, k, d_idx]
      mat[idx, ] <- raw[idx, , drop = FALSE] -
        matrix(F_arr[, k, d_idx], nrow = length(idx), ncol = G, byrow = TRUE)
    }
    centred[[d_idx]] <- mat
  }

  # One bootstrap replication.
  one_boot <- function(seed) {
    if (!is.null(seed)) set.seed(seed)
    W <- sample(c(-1, 1), n, replace = TRUE)
    F_star <- array(0, dim = c(G, K, 2))
    for (d_idx in 1:2) {
      C <- centred[[d_idx]]
      # Per-z bootstrap process at each y_grid point.
      for (k in seq_len(K)) {
        idx <- idx_by_z[[k]]
        if (length(idx) == 0L) next
        F_star[, k, d_idx] <-
          crossprod(W[idx], C[idx, , drop = FALSE])[1, ] / n_z[k]
      }
    }
    compute_stat(F_star)$stat
  }

  # Parallelise on POSIX. CRAN's R CMD check sets _R_CHECK_LIMIT_CORES_
  # and caps at 2 cores; respect that. For interactive use we default
  # to detectCores() - 1, capped at 4 to avoid contention on shared
  # machines.
  use_parallel <- isTRUE(parallel) && .Platform$OS.type == "unix" && n_boot >= 100L
  boot_stats <- if (use_parallel) {
    cran_env <- nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
    cores <- if (cran_env) 2L else min(4L, max(1L, parallel::detectCores() - 1L))
    unlist(parallel::mclapply(seq_len(n_boot),
                              function(b) one_boot(NULL),
                              mc.cores = cores),
           use.names = FALSE)
  } else {
    vapply(seq_len(n_boot), function(b) one_boot(NULL), numeric(1))
  }
  boot_stats <- sqrt(n) * boot_stats

  p_value <- mean(boot_stats >= T_n)

  structure(
    list(
      test = "Kitagawa (2015)",
      statistic = T_n,
      p_value = p_value,
      alpha = alpha,
      n_boot = n_boot,
      boot_stats = boot_stats,
      binding = obs$binding,
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
