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
                               y_grid_size = 50L,
                               d_labels = NULL,
                               treatment_order = c("ordered", "unordered"),
                               monotonicity_set = NULL,
                               multiplier = c("rademacher", "gaussian",
                                              "mammen")) {
  weighting <- match.arg(weighting)
  treatment_order <- match.arg(treatment_order)
  multiplier <- match.arg(multiplier)
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
  unordered <- multivalued && treatment_order == "unordered"

  # Unordered multivalued requires a user-specified monotonicity set C
  # (Sun 2023 Assumption 2.4(iii)). C is a data frame with columns
  # d, z_from, z_to; for each row the test uses the inequality
  # P(Y in B, D = d | Z = z_to) <= P(Y in B, D = d | Z = z_from),
  # corresponding to a one-sided monotonicity restriction that the
  # analyst is willing to impose. Without C we cannot proceed.
  if (unordered) {
    if (is.null(monotonicity_set)) {
      cli::cli_abort(c(
        "{.arg monotonicity_set} is required when {.arg treatment_order} is {.val unordered}.",
        i = "Supply a {.cls data.frame} with columns {.val d}, {.val z_from}, {.val z_to} listing the triples for which {.code 1{{D_{{z_to}} = d}} <= 1{{D_{{z_from}} = d}}} holds almost surely (Sun 2023 Assumption 2.4(iii))."
      ))
    }
    required <- c("d", "z_from", "z_to")
    if (!all(required %in% names(monotonicity_set))) {
      cli::cli_abort(
        "{.arg monotonicity_set} must have columns {.val {required}}."
      )
    }
    bad_d <- !(monotonicity_set$d %in% d_vals)
    bad_zf <- !(monotonicity_set$z_from %in% sort(unique(z_num)))
    bad_zt <- !(monotonicity_set$z_to %in% sort(unique(z_num)))
    if (any(bad_d)) {
      cli::cli_abort("{.arg monotonicity_set$d} contains values not in the observed treatment levels.")
    }
    if (any(bad_zf) || any(bad_zt)) {
      cli::cli_abort("{.arg monotonicity_set$z_from} / {.arg z_to} contain values not in the observed instrument levels.")
    }
  } else if (!is.null(monotonicity_set)) {
    cli::cli_warn(
      "{.arg monotonicity_set} is ignored when {.arg treatment_order} is {.val ordered}; the ordered case uses cumulative-tail inequalities across all adjacent pairs."
    )
  }

  # Order Z levels ascending by first-stage E[D | Z].
  raw_levels <- sort(unique(z_num))
  e_d_by_z <- vapply(raw_levels,
                     function(zk) sum(w[z_num == zk] * d_num[z_num == zk]) /
                                  sum(w[z_num == zk]),
                     numeric(1))
  ord <- order(e_d_by_z)
  z_levels <- raw_levels[ord]
  K <- length(z_levels)

  # Defensive check: the multiplier bootstrap and the plug-in SE
  # estimator can be ill-conditioned when some Z cells have very few
  # observations. Warn if the smallest cell is below 30; abort if
  # below 5.
  n_per_cell <- vapply(z_levels,
                       function(zk) sum(z_num == zk), integer(1))
  if (min(n_per_cell) < 5L) {
    cli::cli_abort(c(
      "Smallest Z cell has {min(n_per_cell)} observations; bootstrap is unreliable.",
      i = "Consider pooling neighbouring Z levels or collecting more data."
    ))
  }
  if (min(n_per_cell) < 30L) {
    cli::cli_warn(c(
      "Smallest Z cell has {min(n_per_cell)} observations (< 30).",
      i = "Bootstrap p-values and the plug-in SE may be ill-conditioned at this cell size. Results should be treated with caution."
    ))
  }

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
  # For unordered multivalued we use per-level indicators P(Y <= y, D = d | z)
  # for each d in d_vals, and the inequality direction is set pairwise
  # from the user-supplied monotonicity set C.
  n_L_planes <- if (unordered) (L + 1L) else if (multivalued) 2L * L else 2L
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
    } else if (unordered) {
      # Per-level indicators: F_arr[, k, ell] = P(Y <= y, D = d_vals[ell] | z).
      for (ell in seq_len(L + 1L)) {
        ind_eq <- (d_sub == d_vals[ell])
        F_arr[, k, ell] <- colSums(w_sub * leY_sub * ind_eq) / sw_k
      }
    } else {
      # Ordered multivalued D: cumulative-tail indicators
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

  # Marginal P(D <= c | Z) stochastic dominance (Sun 2023 equation 10,
  # second inequality). Only meaningful for the ordered multivalued
  # path: for each non-trivial cumulative threshold c and adjacent
  # pair (z_k, z_{k+1}), we test
  #   P(D <= c | z_{k+1}) - P(D <= c | z_k) <= 0.
  # For binary D and unordered multivalued D this machinery is idle.
  do_marginal_sd <- multivalued && !unordered
  marg_planes <- if (do_marginal_sd) L else 0L
  # p_marg[k, m] = P_hat(D <= d_vals[m] | Z = z_k) for m = 1..L.
  p_marg <- NULL
  se_marg <- NULL
  if (do_marginal_sd) {
    p_marg <- matrix(0, nrow = K, ncol = marg_planes)
    for (k in seq_len(K)) {
      idx <- idx_by_z[[k]]
      if (length(idx) == 0L) next
      sw_k <- sum(w[idx])
      for (m in seq_len(marg_planes)) {
        p_marg[k, m] <- sum(w[idx] * (d_num[idx] <= d_vals[m])) / sw_k
      }
    }
    # Plug-in binomial-mixture SE per adjacent pair per threshold.
    se_marg <- array(se_floor, dim = c(K - 1L, marg_planes))
    for (k in seq_len(K - 1L)) {
      for (m in seq_len(marg_planes)) {
        p_lo <- p_marg[k, m]
        p_hi <- p_marg[k + 1L, m]
        total_eff <- n_z_eff[k] + n_z_eff[k + 1L]
        w_lo <- n_z_eff[k + 1L] / total_eff
        w_hi <- n_z_eff[k] / total_eff
        v <- w_lo * p_lo * (1 - p_lo) + w_hi * p_hi * (1 - p_hi)
        se_marg[k, m] <- max(sqrt(max(v, 0)), se_floor)
      }
    }
  }

  # Precompute the triples (plane_idx, k_from, k_to) that get tested.
  # For ordered / binary this is derived from the Z ordering; for
  # unordered it is read from monotonicity_set.
  unordered_triples <- NULL
  if (unordered) {
    d_to_plane <- stats::setNames(seq_along(d_vals), d_vals)
    z_to_k <- stats::setNames(seq_along(z_levels), z_levels)
    plane_idx <- as.integer(d_to_plane[as.character(monotonicity_set$d)])
    # Labels use user's original D coding where available.
    label_d <- if (!is.null(d_labels)) {
      as.character(d_labels[plane_idx])
    } else {
      as.character(monotonicity_set$d)
    }
    unordered_triples <- data.frame(
      plane  = plane_idx,
      k_from = as.integer(z_to_k[as.character(monotonicity_set$z_from)]),
      k_to   = as.integer(z_to_k[as.character(monotonicity_set$z_to)]),
      label_d = label_d,
      stringsAsFactors = FALSE
    )
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
  # inequality and -1 for the "upper tail"). Labels use the user's
  # original D coding where available (d_labels), not the internal 0..k-1
  # remap.
  label_for <- function(i) {
    if (!is.null(d_labels) && length(d_labels) >= i) d_labels[i] else d_vals[i]
  }
  direction_info <- function() {
    if (!multivalued) {
      list(
        list(plane = 2L, label = paste0("d = ", label_for(2L)),
             direction = +1L),  # low - high >= 0
        list(plane = 1L, label = paste0("d = ", label_for(1L)),
             direction = -1L)   # high - low >= 0
      )
    } else {
      out <- vector("list", 2L * L)
      for (ell in seq_len(L)) {
        out[[ell]] <- list(plane = ell,
                           label = paste0("D <= ", label_for(ell)),
                           direction = -1L)
        out[[L + ell]] <- list(plane = L + ell,
                               label = paste0("D >= ", label_for(ell + 1L)),
                               direction = +1L)
      }
      out
    }
  }
  dirs <- direction_info()

  # Marginal SD contribution to the statistic. Takes a K x marg_planes
  # matrix of per-cell cumulative probabilities; used for both the
  # observed statistic (with p_marg) and the bootstrap (with a centred
  # perturbation). Returns the maximum contribution and, when binding,
  # a description of the worst pair/threshold.
  compute_marg_stat <- function(p_mat) {
    if (!do_marginal_sd) return(list(stat = 0, binding = NULL))
    best_m <- 0
    binding_m <- NULL
    for (k in seq_len(K - 1L)) {
      pair_scale <- sqrt(n_z_eff[k] * n_z_eff[k + 1L] /
                           (n_z_eff[k] + n_z_eff[k + 1L]))
      for (m in seq_len(marg_planes)) {
        diff <- p_mat[k + 1L, m] - p_mat[k, m]  # <= 0 under H0
        v <- pair_scale * diff / se_marg[k, m]
        if (v > best_m) {
          best_m <- v
          binding_m <- list(
            z_low = z_levels[k], z_high = z_levels[k + 1L],
            direction = paste0("D <= ", label_for(m)),
            y_lower = NA_real_, y_upper = NA_real_,
            marginal = TRUE
          )
        }
      }
    }
    list(stat = best_m, binding = binding_m)
  }

  compute_stat <- function(F_arr_local, p_marg_local = p_marg) {
    best <- 0
    binding_local <- NULL
    if (K < 2L) return(list(stat = 0, binding = NULL))
    if (unordered) {
      # One triple per row of monotonicity_set: test
      #   P(Y in B, D = d | Z = z_to) <= P(Y in B, D = d | Z = z_from)
      # which corresponds to phi(h, g) = P(.. | z_to) - P(.. | z_from) <= 0.
      for (row_idx in seq_len(nrow(unordered_triples))) {
        plane_idx <- unordered_triples$plane[row_idx]
        k_from   <- unordered_triples$k_from[row_idx]
        k_to     <- unordered_triples$k_to[row_idx]
        pair_scale <- sqrt(n_z_eff[k_from] * n_z_eff[k_to] /
                             (n_z_eff[k_from] + n_z_eff[k_to]))
        P_from <- interval_diff(F_arr_local, k_from, plane_idx)
        P_to   <- interval_diff(F_arr_local, k_to,   plane_idx)
        diff_mat <- P_to - P_from
        se_lo <- min(k_from, k_to); se_hi <- max(k_from, k_to)
        SE_m <- pair_interval_se(se_lo, se_hi, plane_idx)
        UT <- upper.tri(diff_mat, diag = TRUE)
        V <- pair_scale * diff_mat / SE_m
        V[!UT] <- -Inf
        m <- max(V)
        if (m > best) {
          best <- m
          pos <- which(V == m, arr.ind = TRUE)[1, ]
          binding_local <- list(
            z_from = z_levels[k_from],
            z_to = z_levels[k_to],
            direction = paste0("d = ", unordered_triples$label_d[row_idx]),
            y_lower = if (pos[1] == 1L) -Inf else y_grid[pos[1] - 1L],
            y_upper = y_grid[pos[2] - 1L]
          )
        }
      }
    } else for (k_low in seq_len(K - 1L)) {
      for (k_high in (k_low + 1L):K) {
        # Kitagawa (2015) eq. 2.1: the statistic for pair (z_low, z_high)
        # is sqrt(n_low * n_high / (n_low + n_high)) times the sup of
        # the positive-part difference divided by the plug-in SE.
        pair_scale <- sqrt(n_z_eff[k_low] * n_z_eff[k_high] /
                             (n_z_eff[k_low] + n_z_eff[k_high]))
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
          V <- pair_scale * diff_mat / SE_m
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
    # Add marginal P(D <= c | Z) stochastic dominance (Sun eq. 10,
    # second inequality) to the sup if applicable. Uses the supplied
    # p_marg_local matrix so bootstrap can inject its own perturbation.
    if (do_marginal_sd) {
      marg <- compute_marg_stat(p_marg_local)
      if (marg$stat > best) {
        best <- marg$stat
        binding_local <- marg$binding
      }
    }
    list(stat = best, binding = binding_local)
  }

  obs <- compute_stat(F_arr)
  # Pair-specific scaling is now applied inside compute_stat, so obs$stat
  # is already on the Kitagawa (eq. 2.1) scale for the variance-weighted
  # form. The unweighted form still gets a global sqrt(n) rescaling to
  # match eq. 2.2.
  n_eff_total <- sum(n_z_eff)
  scale_n <- if (weighting == "variance") 1 else sqrt(n_eff_total)
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
    } else if (unordered) {
      raw <- leY * (d_num == d_vals[d_idx])
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

  draw_multiplier <- switch(multiplier,
    rademacher = function(nn) sample(c(-1, 1), nn, replace = TRUE),
    gaussian   = function(nn) stats::rnorm(nn),
    mammen     = function(nn) {
      # Mammen (1993) two-point distribution: mean 0, variance 1, skewness 1.
      phi <- (1 + sqrt(5)) / 2  # golden ratio
      p   <- (sqrt(5) + 1) / (2 * sqrt(5))
      ifelse(stats::runif(nn) < p, -(phi - 1), phi)
    }
  )

  # Centred indicators for the marginal P(D <= c | Z) bootstrap:
  # marg_centred[i, m] = I{D_i <= d_vals[m]} - P_hat(D <= d_vals[m] | Z_i).
  marg_centred <- NULL
  if (do_marginal_sd) {
    marg_centred <- matrix(0, nrow = n, ncol = marg_planes)
    for (k in seq_len(K)) {
      idx <- idx_by_z[[k]]
      if (length(idx) == 0L) next
      for (m in seq_len(marg_planes)) {
        marg_centred[idx, m] <- w[idx] *
          ((d_num[idx] <= d_vals[m]) - p_marg[k, m])
      }
    }
  }

  one_boot <- function() {
    W <- draw_multiplier(n)
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
    p_marg_star <- NULL
    if (do_marginal_sd) {
      p_marg_star <- matrix(0, nrow = K, ncol = marg_planes)
      for (k in seq_len(K)) {
        idx <- idx_by_z[[k]]
        if (length(idx) == 0L) next
        p_marg_star[k, ] <-
          crossprod(W[idx], marg_centred[idx, , drop = FALSE])[1, ] /
          n_z_eff[k]
      }
    }
    compute_stat(F_star, p_marg_local = p_marg_star)$stat
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
