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
                               weights = NULL,
                               se_floor = 0.001,
                               y_grid_size = 50L) {
  weighting <- match.arg(weighting)
  n <- length(y)

  if (is.null(weights)) {
    w <- rep(1, n)
  } else {
    if (length(weights) != n) {
      cli::cli_abort("{.arg weights} must have length equal to the sample size.")
    }
    if (any(weights < 0) || any(!is.finite(weights))) {
      cli::cli_abort("{.arg weights} must be finite and non-negative.")
    }
    # Rescale so mean = 1 (preserves effective sample size interpretation)
    w <- weights * n / sum(weights)
  }

  # Order Z levels ascending by first-stage E[D | Z].
  raw_levels <- sort(unique(z_num))
  e_d_by_z <- vapply(raw_levels,
                     function(zk) mean(d_num[z_num == zk]),
                     numeric(1))
  ord <- order(e_d_by_z)
  z_levels <- raw_levels[ord]
  K <- length(z_levels)

  # Evaluation grid. Per Kitagawa (2015) eq. 2.1 the statistic is a sup
  # over intervals [y, y'] with y <= y'. For efficiency we evaluate on a
  # quantile grid; equivalent under Kitagawa's theorem 2.1 because the
  # sup is attained at sample Y values and the quantile grid is dense in
  # the observed support.
  y_unique <- sort(unique(y))
  if (length(y_unique) <= y_grid_size) {
    y_grid <- y_unique
  } else {
    probs <- seq(0, 1, length.out = y_grid_size)
    y_grid <- sort(unique(stats::quantile(y, probs = probs, names = FALSE)))
  }
  G <- length(y_grid)

  leY <- outer(y, y_grid, "<=")

  # Joint CDF F_arr[g, k, d_idx] = F_hat(y_grid[g], d | z_levels[k])
  F_arr <- array(0, dim = c(G, K, 2))
  n_z <- integer(K)
  idx_by_z <- vector("list", K)

  # Effective sample sizes (used for scaling under weighting).
  n_z_eff <- numeric(K)
  for (k in seq_len(K)) {
    idx <- which(z_num == z_levels[k])
    idx_by_z[[k]] <- idx
    n_z[k] <- length(idx)
    if (n_z[k] == 0L) next
    d_sub <- d_num[idx]
    leY_sub <- leY[idx, , drop = FALSE]
    w_sub <- w[idx]
    sw_k <- sum(w_sub)
    n_z_eff[k] <- sw_k
    # Weighted empirical joint CDF: (1/sum w_i) sum w_i * 1{...}
    F_arr[, k, 2] <- colSums(w_sub * leY_sub * d_sub) / sw_k
    F_arr[, k, 1] <- colSums(w_sub * leY_sub * (1 - d_sub)) / sw_k
  }

  # Interval probabilities P_k_d[g1, g2] = F_hat(y_grid[g2], d | z_k) -
  # F_hat(y_grid[g1], d | z_k) for g1 <= g2. Corresponds to
  # P(Y in (y_grid[g1], y_grid[g2]], D = d | Z = z_k) plus the point
  # probabilities at the grid endpoints (by convention for empirical
  # CDFs). We also cover the lower-tail case by prepending a phantom
  # -infinity grid point with F = 0.
  interval_diff <- function(F_arr_local, k, d_idx) {
    f <- F_arr_local[, k, d_idx]
    # Prepend 0 at the start to include intervals (-inf, y']
    f_ext <- c(0, f)
    # G+1 x G+1 matrix where entry (g1+1, g2+1) = f_ext[g2+1] - f_ext[g1+1]
    # for g1 < g2 (strict lower triangular). For efficiency we compute
    # only the upper triangle via outer.
    outer(f_ext, f_ext, "-") * -1  # P[g1, g2] = F(g2) - F(g1)
  }

  # Pre-compute the data-derived SE matrix once per (k_low, k_high, d).
  # Kitagawa (2015) eq. 2.1: with lambda_k = n_k / n for k in
  # {low, high}, sigma^2([y, y'], d) =
  #   (n_high / n) * P_low  * (1 - P_low)
  # + (n_low  / n) * P_high * (1 - P_high).
  # That is, the variance of sqrt(n_low * n_high / n) * (P_low - P_high)
  # asymptotic to a binomial mixture. Both the observed statistic and
  # every bootstrap statistic are normalised by this plug-in SE.
  SE_cache <- array(0, dim = c(G + 1L, G + 1L, K, K, 2))
  if (weighting == "variance") {
    for (k_low in seq_len(K - 1L)) {
      for (k_high in (k_low + 1L):K) {
        for (d_idx in 1:2) {
          P_low  <- interval_diff(F_arr, k_low,  d_idx)
          P_high <- interval_diff(F_arr, k_high, d_idx)
          # Weighted mixture variance using effective sample sizes.
          total_eff <- n_z_eff[k_low] + n_z_eff[k_high]
          w_low  <- n_z_eff[k_high] / total_eff
          w_high <- n_z_eff[k_low]  / total_eff
          var_mat <- w_low  * P_low  * (1 - P_low)  +
                     w_high * P_high * (1 - P_high)
          SE_cache[, , k_low, k_high, d_idx] <- pmax(sqrt(pmax(var_mat, 0)),
                                                     se_floor)
        }
      }
    }
  } else {
    SE_cache[] <- 1
  }

  pair_interval_se <- function(k_low, k_high, d_idx) {
    SE_cache[, , k_low, k_high, d_idx]
  }

  compute_stat <- function(F_arr_local) {
    best <- 0
    binding_local <- NULL
    if (K < 2L) return(list(stat = 0, binding = NULL))
    for (k_low in seq_len(K - 1L)) {
      for (k_high in (k_low + 1L):K) {
        # d = 1: P(Y in B, D=1 | z_low) - P(Y in B, D=1 | z_high) <= 0
        Pdiff_1_low  <- interval_diff(F_arr_local, k_low,  2)
        Pdiff_1_high <- interval_diff(F_arr_local, k_high, 2)
        viol_1 <- Pdiff_1_low - Pdiff_1_high
        SE_1   <- pair_interval_se(k_low, k_high, 2)
        # Only consider upper triangle (g1 < g2).
        UT <- upper.tri(viol_1, diag = TRUE)
        V1 <- viol_1 / SE_1
        V1[!UT] <- -Inf
        m1 <- max(V1)
        if (m1 > best) {
          best <- m1
          pos <- which(V1 == m1, arr.ind = TRUE)[1, ]
          binding_local <- list(
            z_low = z_levels[k_low],
            z_high = z_levels[k_high],
            d = 1L,
            y_lower = if (pos[1] == 1L) -Inf else y_grid[pos[1] - 1L],
            y_upper = y_grid[pos[2] - 1L]
          )
        }
        # d = 0: P(Y in B, D=0 | z_high) - P(Y in B, D=0 | z_low) <= 0
        Pdiff_0_low  <- interval_diff(F_arr_local, k_low,  1)
        Pdiff_0_high <- interval_diff(F_arr_local, k_high, 1)
        viol_0 <- Pdiff_0_high - Pdiff_0_low
        SE_0   <- pair_interval_se(k_low, k_high, 1)
        V0 <- viol_0 / SE_0
        V0[!UT] <- -Inf
        m0 <- max(V0)
        if (m0 > best) {
          best <- m0
          pos <- which(V0 == m0, arr.ind = TRUE)[1, ]
          binding_local <- list(
            z_low = z_levels[k_low],
            z_high = z_levels[k_high],
            d = 0L,
            y_lower = if (pos[1] == 1L) -Inf else y_grid[pos[1] - 1L],
            y_upper = y_grid[pos[2] - 1L]
          )
        }
      }
    }
    list(stat = best, binding = binding_local)
  }

  obs <- compute_stat(F_arr)
  # Kitagawa (eq. 2.1) scales by sqrt(n_low * n_high / n_total). Under
  # survey weights the effective sample sizes replace raw counts.
  n_eff_total <- sum(n_z_eff)
  scale_n <- if (weighting == "variance") {
    sqrt(n_eff_total / K)
  } else {
    sqrt(n_eff_total)
  }
  T_n <- scale_n * obs$stat

  # Precompute weight-scaled centred indicator matrices for the
  # multiplier bootstrap. Each observation's contribution to the bootstrap
  # F^* process is w_i * (indicator - F_hat) multiplied by a Rademacher
  # weight. Normalisation uses the per-judge sum of weights (sw_k).
  centred <- list()
  for (d_idx in 1:2) {
    d_val <- d_idx - 1L
    mat <- matrix(0, nrow = n, ncol = G)
    raw <- leY * (d_num == d_val)
    for (k in seq_len(K)) {
      idx <- idx_by_z[[k]]
      if (length(idx) == 0L) next
      # w * (indicator - F_hat)
      mat[idx, ] <- w[idx] * (raw[idx, , drop = FALSE] -
        matrix(F_arr[, k, d_idx], nrow = length(idx), ncol = G, byrow = TRUE))
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
          crossprod(W[idx], C[idx, , drop = FALSE])[1, ] / n_z_eff[k]
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
