# Shared computational core for pairwise one-sided KS-type tests of
# the Imbens-Angrist (1994) IV-validity inequalities.
#
# Used internally by:
#   * iv_kitagawa() (unconditional case)
#   * iv_mw()       (both unconditional and X-stratified cases)
#
# Supports both the unweighted form (Kitagawa 2015 section 3) and the
# variance-weighted form (Kitagawa 2015 section 4). The variance-
# weighted form divides each pointwise difference by its plug-in
# standard error estimator before the sup-over-(y, pairs, d), which
# gives better finite-sample power when cell sizes are unequal.
#
# Not exported.

#' @noRd
kitagawa_core_test <- function(y, d_num, z_num, n_boot, parallel,
                               weighting = c("variance", "unweighted"),
                               se_floor = 0.001) {
  weighting <- match.arg(weighting)
  n <- length(y)

  # Order Z levels ascending by first-stage E[D | Z].
  raw_levels <- sort(unique(z_num))
  e_d_by_z <- vapply(raw_levels,
                     function(zk) mean(d_num[z_num == zk]),
                     numeric(1))
  ord <- order(e_d_by_z)
  z_levels <- raw_levels[ord]
  K <- length(z_levels)

  y_grid <- sort(unique(y))
  G <- length(y_grid)

  # Indicator matrix leY[i, g] = 1{Y_i <= y_grid[g]}
  leY <- outer(y, y_grid, "<=")

  # F_arr[g, k, d_idx] = F_hat(y_grid[g], d | z_levels[k])
  F_arr <- array(0, dim = c(G, K, 2))
  n_z <- integer(K)
  idx_by_z <- vector("list", K)

  for (k in seq_len(K)) {
    idx <- which(z_num == z_levels[k])
    idx_by_z[[k]] <- idx
    n_z[k] <- length(idx)
    if (n_z[k] == 0L) next
    d_sub <- d_num[idx]
    leY_sub <- leY[idx, , drop = FALSE]
    F_arr[, k, 2] <- colSums(leY_sub * d_sub) / n_z[k]
    F_arr[, k, 1] <- colSums(leY_sub * (1 - d_sub)) / n_z[k]
  }

  # SE matrix: SE[g, k, d_idx] = sqrt(n) * sd(F_hat(y, d | z))
  # under binomial sampling, = sqrt(n * F * (1 - F) / n_z)
  # = sqrt((n / n_z) * F * (1 - F)).
  # The variance of sqrt(n) * (F_low - F_high) is SE_low^2 + SE_high^2
  # (independent strata). Floor ensures we don't divide by zero where
  # both F's are 0 or 1.
  SE_arr <- if (weighting == "variance") {
    out <- array(0, dim = c(G, K, 2))
    for (k in seq_len(K)) {
      if (n_z[k] == 0L) next
      rf <- n / n_z[k]
      for (d_idx in 1:2) {
        Fv <- F_arr[, k, d_idx]
        out[, k, d_idx] <- sqrt(pmax(rf * Fv * (1 - Fv), 0))
      }
    }
    out
  } else {
    NULL
  }

  pair_se <- function(k_low, k_high, d_idx) {
    if (is.null(SE_arr)) return(1)
    s <- sqrt(SE_arr[, k_low, d_idx]^2 + SE_arr[, k_high, d_idx]^2)
    pmax(s, se_floor)
  }

  compute_stat <- function(F_arr_local) {
    best <- 0
    binding_local <- NULL
    if (K < 2L) return(list(stat = 0, binding = NULL))
    for (k_low in seq_len(K - 1L)) {
      for (k_high in (k_low + 1L):K) {
        # d = 1 inequality: F(y, 1 | low) - F(y, 1 | high) <= 0
        scale1 <- pair_se(k_low, k_high, 2)
        diffs1 <- (F_arr_local[, k_low, 2] - F_arr_local[, k_high, 2]) / scale1
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
        # d = 0 inequality: F(y, 0 | high) - F(y, 0 | low) <= 0
        scale0 <- pair_se(k_low, k_high, 1)
        diffs0 <- (F_arr_local[, k_high, 1] - F_arr_local[, k_low, 1]) / scale0
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
  # For the unweighted form we still scale by sqrt(n). For the
  # variance-weighted form the scaling is already absorbed into the SE
  # denominator (SE contains the sqrt(n) rescale).
  scale_n <- if (weighting == "variance") 1 else sqrt(n)
  T_n <- scale_n * obs$stat

  # Precompute centred indicator matrices for the multiplier bootstrap.
  centred <- list()
  for (d_idx in 1:2) {
    d_val <- d_idx - 1L
    mat <- matrix(0, nrow = n, ncol = G)
    raw <- leY * (d_num == d_val)
    for (k in seq_len(K)) {
      idx <- idx_by_z[[k]]
      if (length(idx) == 0L) next
      mat[idx, ] <- raw[idx, , drop = FALSE] -
        matrix(F_arr[, k, d_idx], nrow = length(idx), ncol = G, byrow = TRUE)
    }
    centred[[d_idx]] <- mat
  }

  one_boot <- function() {
    W <- sample(c(-1, 1), n, replace = TRUE)
    F_star <- array(0, dim = c(G, K, 2))
    for (d_idx in 1:2) {
      C <- centred[[d_idx]]
      for (k in seq_len(K)) {
        idx <- idx_by_z[[k]]
        if (length(idx) == 0L) next
        F_star[, k, d_idx] <-
          crossprod(W[idx], C[idx, , drop = FALSE])[1, ] / n_z[k]
      }
    }
    compute_stat(F_star)$stat
  }

  use_parallel <- isTRUE(parallel) && .Platform$OS.type == "unix" && n_boot >= 100L
  boot_stats <- if (use_parallel) {
    cran_env <- nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
    cores <- if (cran_env) 2L else min(4L, max(1L, parallel::detectCores() - 1L))
    unlist(parallel::mclapply(seq_len(n_boot),
                              function(b) one_boot(),
                              mc.cores = cores),
           use.names = FALSE)
  } else {
    vapply(seq_len(n_boot), function(b) one_boot(), numeric(1))
  }
  boot_stats <- scale_n * boot_stats

  p_value <- mean(boot_stats >= T_n)

  list(
    statistic = T_n,
    p_value = p_value,
    boot_stats = boot_stats,
    binding = obs$binding,
    n = n,
    weighting = weighting
  )
}
