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
                               se_floor = 0.15,
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
    w <- weights * n / sum(weights)
  }

  # Identify treatment levels. For binary D (Kitagawa 2015) thresholds
  # reduce to L = 1. For ordered multivalued D (Sun 2023) there are L+1
  # levels {0, 1, ..., L}; we test stochastic-dominance inequalities on
  # cumulative tails P(Y <= y, D >= ell | Z) for ell = 1, ..., L and
  # P(Y <= y, D <= ell | Z) for ell = 0, ..., L-1.
  d_vals <- sort(unique(d_num))
  L <- length(d_vals) - 1L
  multivalued <- L > 1L

  # Order Z levels ascending by first-stage E[D | Z].
  raw_levels <- sort(unique(z_num))
  e_d_by_z <- vapply(raw_levels,
                     function(zk) sum(w[z_num == zk] * d_num[z_num == zk]) /
                                  sum(w[z_num == zk]),
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

  # Joint CDF F_arr[g, k, l_idx] encodes the necessary conditional CDFs.
  # For binary D (L = 1), we use two indicators:
  #   F_arr[, , 1] = F_hat(y, D = 0 | z) (needed for d = 0 inequality)
  #   F_arr[, , 2] = F_hat(y, D = 1 | z) (needed for d = 1 inequality)
  # For multivalued ordered D (L > 1), we instead store the cumulative-tail
  # families:
  #   F_arr[, , ell]         = F_hat(y, D <= d_vals[ell]  | z)  for ell = 1..L
  #   F_arr[, , L + ell]     = F_hat(y, D >= d_vals[ell+1] | z) for ell = 1..L
  # Each "direction" gives L inequalities per instrument pair.
  n_L_planes <- if (multivalued) 2L * L else 2L
  F_arr <- array(0, dim = c(G, K, n_L_planes))
  n_z <- integer(K)
  idx_by_z <- vector("list", K)
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
    if (!multivalued) {
      F_arr[, k, 2] <- colSums(w_sub * leY_sub * d_sub) / sw_k
      F_arr[, k, 1] <- colSums(w_sub * leY_sub * (1 - d_sub)) / sw_k
    } else {
      # Ordered multivalued D (Sun 2023). For each threshold ell = 1..L:
      #   lower plane: P(Y <= y, D <= d_vals[ell] | z)
      #   upper plane: P(Y <= y, D >= d_vals[ell + 1] | z)
      for (ell in seq_len(L)) {
        d_thr <- d_vals[ell]
        ind_le <- (d_sub <= d_thr)
        F_arr[, k, ell] <- colSums(w_sub * leY_sub * ind_le) / sw_k
        d_thr_up <- d_vals[ell + 1L]
        ind_ge <- (d_sub >= d_thr_up)
        F_arr[, k, L + ell] <- colSums(w_sub * leY_sub * ind_ge) / sw_k
      }
    }
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

  # Pre-compute the data-derived SE matrix. Same binomial-mixture logic
  # per plane; the n_L_planes dimension replaces the old "2" (d) axis.
  SE_cache <- array(0, dim = c(G + 1L, G + 1L, K, K, n_L_planes))
  if (weighting == "variance") {
    for (k_low in seq_len(K - 1L)) {
      for (k_high in (k_low + 1L):K) {
        for (d_idx in seq_len(n_L_planes)) {
          P_low  <- interval_diff(F_arr, k_low,  d_idx)
          P_high <- interval_diff(F_arr, k_high, d_idx)
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

  # Violation-direction descriptor: for binary D these are literally
  # (d = 0, d = 1); for multivalued D they are cumulative-tail thresholds
  # (D <= ell vs D >= ell+1). Violations have opposite sign between the
  # two directions, encoded by `sign_low` (= +1 for the "low tail"
  # inequality and -1 for the "upper tail").
  direction_info <- function() {
    if (!multivalued) {
      list(
        list(plane = 2L, label = "d = 1", direction = +1L),  # low - high >= 0
        list(plane = 1L, label = "d = 0", direction = -1L)   # high - low >= 0
      )
    } else {
      out <- vector("list", 2L * L)
      for (ell in seq_len(L)) {
        out[[ell]] <- list(plane = ell, label = paste0("D <= ", d_vals[ell]),
                           direction = -1L)
        out[[L + ell]] <- list(plane = L + ell,
                               label = paste0("D >= ", d_vals[ell + 1L]),
                               direction = +1L)
      }
      out
    }
  }
  dirs <- direction_info()

  compute_stat <- function(F_arr_local) {
    best <- 0
    binding_local <- NULL
    if (K < 2L) return(list(stat = 0, binding = NULL))
    for (k_low in seq_len(K - 1L)) {
      for (k_high in (k_low + 1L):K) {
        for (dd in dirs) {
          P_low  <- interval_diff(F_arr_local, k_low,  dd$plane)
          P_high <- interval_diff(F_arr_local, k_high, dd$plane)
          diff_mat <- if (dd$direction > 0) {
            P_low - P_high
          } else {
            P_high - P_low
          }
          SE_m <- pair_interval_se(k_low, k_high, dd$plane)
          UT <- upper.tri(diff_mat, diag = TRUE)
          V <- diff_mat / SE_m
          V[!UT] <- -Inf
          m <- max(V)
          if (m > best) {
            best <- m
            pos <- which(V == m, arr.ind = TRUE)[1, ]
            binding_local <- list(
              z_low = z_levels[k_low],
              z_high = z_levels[k_high],
              direction = dd$label,
              y_lower = if (pos[1] == 1L) -Inf else y_grid[pos[1] - 1L],
              y_upper = y_grid[pos[2] - 1L]
            )
          }
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

  # Precompute centred indicator matrices for the multiplier bootstrap.
  # Binary D: two planes (D = 0, D = 1). Multivalued D: 2L planes matching
  # the cumulative-tail structure of F_arr.
  centred <- vector("list", n_L_planes)
  for (d_idx in seq_len(n_L_planes)) {
    mat <- matrix(0, nrow = n, ncol = G)
    if (!multivalued) {
      d_val <- d_idx - 1L
      raw <- leY * (d_num == d_val)
    } else if (d_idx <= L) {
      # lower plane: D <= d_vals[ell] where ell = d_idx
      raw <- leY * (d_num <= d_vals[d_idx])
    } else {
      # upper plane: D >= d_vals[ell + 1] where ell = d_idx - L
      ell <- d_idx - L
      raw <- leY * (d_num >= d_vals[ell + 1L])
    }
    for (k in seq_len(K)) {
      idx <- idx_by_z[[k]]
      if (length(idx) == 0L) next
      mat[idx, ] <- w[idx] * (raw[idx, , drop = FALSE] -
        matrix(F_arr[, k, d_idx], nrow = length(idx), ncol = G, byrow = TRUE))
    }
    centred[[d_idx]] <- mat
  }

  one_boot <- function() {
    W <- sample(c(-1, 1), n, replace = TRUE)
    F_star <- array(0, dim = c(G, K, n_L_planes))
    for (d_idx in seq_len(n_L_planes)) {
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
    weighting = weighting,
    n_treatment_levels = L + 1L,
    multivalued = multivalued
  )
}
