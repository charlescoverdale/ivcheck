# Shared computational core for pairwise one-sided KS-type tests of
# the Imbens-Angrist (1994) IV-validity inequalities.
#
# Used internally by:
#   * iv_kitagawa() (unconditional case)
#   * iv_mw()       (both unconditional and X-stratified cases)
#
# Not exported.

#' @noRd
kitagawa_core_test <- function(y, d_num, z_num, n_boot, parallel) {
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

  compute_stat <- function(F_arr_local) {
    best <- 0
    binding_local <- NULL
    if (K < 2L) return(list(stat = 0, binding = NULL))
    for (k_low in seq_len(K - 1L)) {
      for (k_high in (k_low + 1L):K) {
        # d = 1 inequality: F(y, 1 | low) - F(y, 1 | high) <= 0
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
        # d = 0 inequality: F(y, 0 | high) - F(y, 0 | low) <= 0
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
  boot_stats <- sqrt(n) * boot_stats

  p_value <- mean(boot_stats >= T_n)

  list(
    statistic = T_n,
    p_value = p_value,
    boot_stats = boot_stats,
    binding = obs$binding,
    n = n
  )
}
