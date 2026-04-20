# Full CLR-style conditional moment inequality test for iv_mw with
# covariates. Replaces the bin-stratified simplification of v0.1.0-pre
# with series-regression conditional CDF estimation, plug-in
# heteroscedasticity-robust standard errors, and adaptive moment
# selection in the style of Andrews and Soares (2010).
#
# Not exported.

# K-fold cross-validated selection of the polynomial basis order for
# the series-regression conditional CDF estimator. Loss is squared
# prediction error of the indicator 1{Y <= y, D = d} on the
# polynomial(X) basis, averaged across a coarse y-grid, both d levels,
# and every Z cell. Candidates default to 2:5 (quadratic to quintic).
#
# Runtime scales with n_folds x length(candidates); at n_folds = 5 and
# four candidates this is 20 times a single-order fit. Fast for
# typical (n < 10000) designs.
#
# Not exported.
#' @noRd
mw_select_basis_order <- function(y, d_num, z_num, x_mat,
                                  candidates = 2:5,
                                  n_folds = 5L,
                                  y_grid_size = 10L) {
  n <- length(y)
  x_str <- x_mat[, 1]
  d_vals <- sort(unique(d_num))
  z_levels <- sort(unique(z_num))
  K <- length(z_levels)

  # Coarse y-grid for CV; the final test uses the finer y_grid_size.
  y_probs <- seq(0.1, 0.9, length.out = y_grid_size)
  y_grid <- sort(unique(stats::quantile(y, probs = y_probs, names = FALSE)))
  G <- length(y_grid)
  leY <- outer(y, y_grid, "<=")

  make_basis <- function(v, order) {
    out <- matrix(1, nrow = length(v), ncol = order + 1L)
    for (p in seq_len(order)) out[, p + 1L] <- v^p
    out
  }

  # CV fold assignment, stratified within each Z cell so every fold
  # contains observations from every instrument level.
  fold_id <- integer(n)
  for (k in seq_len(K)) {
    idx <- which(z_num == z_levels[k])
    if (length(idx) == 0L) next
    fold_id[idx] <- sample(rep_len(seq_len(n_folds), length(idx)))
  }

  mse_by_order <- numeric(length(candidates))
  for (p_idx in seq_along(candidates)) {
    p <- candidates[p_idx]
    total_sq_err <- 0
    n_terms <- 0L
    for (fold in seq_len(n_folds)) {
      train <- which(fold_id != fold)
      test  <- which(fold_id == fold)
      if (length(train) == 0L || length(test) == 0L) next
      for (k in seq_len(K)) {
        idx_tr <- train[z_num[train] == z_levels[k]]
        idx_te <- test[z_num[test] == z_levels[k]]
        if (length(idx_tr) <= p + 1L || length(idx_te) == 0L) next
        B_tr <- make_basis(x_str[idx_tr], p)
        B_te <- make_basis(x_str[idx_te], p)
        BtB_inv <- tryCatch(
          solve(crossprod(B_tr) + diag(1e-8, p + 1L)),
          error = function(e) NULL
        )
        if (is.null(BtB_inv)) next
        for (dv in d_vals) {
          ind_mask_tr <- d_num[idx_tr] == dv
          ind_mask_te <- d_num[idx_te] == dv
          ind_tr <- leY[idx_tr, , drop = FALSE] * ind_mask_tr
          ind_te <- leY[idx_te, , drop = FALSE] * ind_mask_te
          beta <- BtB_inv %*% (t(B_tr) %*% ind_tr)
          pred_te <- B_te %*% beta
          total_sq_err <- total_sq_err + sum((ind_te - pred_te)^2)
          n_terms <- n_terms + length(ind_te)
        }
      }
    }
    mse_by_order[p_idx] <- if (n_terms > 0L) total_sq_err / n_terms else Inf
  }
  list(
    selected = candidates[which.min(mse_by_order)],
    mse = stats::setNames(mse_by_order, candidates),
    candidates = candidates
  )
}

# Post-selection-valid iv_mw test with CV-selected basis order.
# Runs the core test under each candidate basis order and reports a
# p-value valid against any selection rule by taking the maximum of
# the bootstrap statistics across candidates. Conservative relative
# to a fixed-order test; guaranteed to control size at the nominal
# level regardless of which candidate CV selects.
#' @noRd
mw_clr_test_posi <- function(y, d_num, z_num, x_mat, n_boot, parallel,
                             candidates = 2:5,
                             x_grid_size = 20L, y_grid_size = 50L,
                             adaptive = TRUE) {
  cv <- mw_select_basis_order(y, d_num, z_num, x_mat,
                              candidates = candidates,
                              y_grid_size = 10L)
  p_hat <- cv$selected
  nP <- length(candidates)

  per_order <- vector("list", nP)
  for (i in seq_len(nP)) {
    per_order[[i]] <- mw_clr_test(y, d_num, z_num, x_mat,
                                  n_boot = n_boot, parallel = parallel,
                                  basis_order = candidates[i],
                                  x_grid_size = x_grid_size,
                                  y_grid_size = y_grid_size,
                                  adaptive = adaptive)
  }
  # Pair bootstrap replications across candidates and take the
  # maximum. Using independent Rademacher weights per candidate is
  # (conservatively) slightly larger than a shared-weights joint
  # bootstrap would give; tightening to the exact joint is a
  # v0.2.0 item.
  boot_min_len <- min(vapply(per_order, function(o) length(o$boot_stats),
                             integer(1)))
  boot_mat <- vapply(per_order, function(o) o$boot_stats[seq_len(boot_min_len)],
                     numeric(boot_min_len))
  boot_max <- apply(boot_mat, 1, max)

  sel_idx <- match(p_hat, candidates)
  observed <- per_order[[sel_idx]]$statistic
  p_value <- mean(boot_max >= observed)

  list(
    statistic = observed,
    p_value = p_value,
    boot_stats = boot_max,
    binding = per_order[[sel_idx]]$binding,
    n = per_order[[sel_idx]]$n,
    kappa_n = per_order[[sel_idx]]$kappa_n,
    basis_order = p_hat,
    basis_order_cv = cv,
    post_selection_valid = TRUE
  )
}

#' @noRd
mw_clr_test <- function(y, d_num, z_num, x_mat, n_boot, parallel,
                        basis_order = 3L,
                        x_grid_size = 20L,
                        y_grid_size = 50L,
                        adaptive = TRUE,
                        se_floor = 0.001) {
  n <- length(y)

  raw_levels <- sort(unique(z_num))
  e_d_by_z <- vapply(raw_levels,
                     function(zk) mean(d_num[z_num == zk]),
                     numeric(1))
  ord <- order(e_d_by_z)
  z_levels <- raw_levels[ord]
  K <- length(z_levels)

  x_str <- x_mat[, 1]

  make_basis <- function(v) {
    out <- matrix(1, nrow = length(v), ncol = basis_order + 1L)
    for (p in seq_len(basis_order)) {
      out[, p + 1L] <- v^p
    }
    out
  }
  B <- make_basis(x_str)
  P <- ncol(B)

  y_probs <- seq(0.02, 0.98, length.out = y_grid_size)
  y_grid <- sort(unique(stats::quantile(y, probs = y_probs, names = FALSE)))
  G <- length(y_grid)
  x_probs <- seq(0.1, 0.9, length.out = x_grid_size)
  x_grid <- sort(unique(stats::quantile(x_str, probs = x_probs, names = FALSE)))
  XG <- length(x_grid)
  B_at_xg <- make_basis(x_grid)

  betas <- array(0, dim = c(P, G, K, 2))
  residuals_by <- array(0, dim = c(n, G, 2))
  BtB_inv_by_k <- vector("list", K)
  idx_by_z <- vector("list", K)
  n_z <- integer(K)

  for (k in seq_len(K)) {
    idx <- which(z_num == z_levels[k])
    idx_by_z[[k]] <- idx
    n_z[k] <- length(idx)
    Bk <- B[idx, , drop = FALSE]
    BtB_inv <- tryCatch(
      solve(crossprod(Bk) + diag(1e-8, P)),
      error = function(e) {
        cli::cli_abort(c(
          "Series basis is rank-deficient in the Z = {z_levels[k]} cell (n = {length(idx)}, basis_order = {basis_order}).",
          i = "Reduce {.arg basis_order} or ensure {.arg x} has adequate within-cell variation.",
          i = "Fewer than {.val {P}} distinct values of {.arg x} within a {.arg z} level forces this error."
        ))
      }
    )
    BtB_inv_by_k[[k]] <- BtB_inv
    Xproj <- BtB_inv %*% t(Bk)

    for (d_idx in 1:2) {
      d_val <- d_idx - 1L
      leY_k <- outer(y[idx], y_grid, "<=")
      ind <- leY_k * (d_num[idx] == d_val)
      beta_kd <- Xproj %*% ind
      betas[, , k, d_idx] <- beta_kd
      fit <- Bk %*% beta_kd
      residuals_by[idx, , d_idx] <- ind - fit
    }
  }

  F_on_grid <- array(0, dim = c(XG, G, K, 2))
  SE_on_grid <- array(se_floor, dim = c(XG, G, K, 2))

  for (k in seq_len(K)) {
    idx <- idx_by_z[[k]]
    Bk <- B[idx, , drop = FALSE]
    BtB_inv <- BtB_inv_by_k[[k]]
    for (d_idx in 1:2) {
      F_on_grid[, , k, d_idx] <- B_at_xg %*% betas[, , k, d_idx]

      Rk <- residuals_by[idx, , d_idx]
      # Robust variance. For each y_g compute
      # Var_g = BtB_inv %*% (Bk' diag(r_g^2) Bk) %*% BtB_inv
      # and return the diagonal of B_at_xg %*% Var_g %*% t(B_at_xg).
      for (g in seq_len(G)) {
        r <- Rk[, g]
        Omega <- crossprod(Bk * r)
        V <- BtB_inv %*% Omega %*% BtB_inv
        # Quadratic form diag(B_at_xg %*% V %*% t(B_at_xg))
        qf <- rowSums((B_at_xg %*% V) * B_at_xg)
        SE_on_grid[, g, k, d_idx] <- pmax(sqrt(pmax(qf, 0)), se_floor)
      }
    }
  }

  # Observed studentised statistic xi for every moment.
  # For d = 1: xi = sqrt(n) * (F_low - F_high) / SE   (violation if xi > 0)
  # For d = 0: xi = sqrt(n) * (F_high - F_low) / SE
  # Pool SE as sqrt(SE_low^2 + SE_high^2) per pair.
  compute_stat <- function(F_array) {
    best <- -Inf
    binding_local <- NULL
    if (K < 2L) return(list(stat = 0, binding = NULL, xi_all = NULL))
    xi_list <- list()
    for (k_low in seq_len(K - 1L)) {
      for (k_high in (k_low + 1L):K) {
        for (d_idx in 1:2) {
          se_low  <- SE_on_grid[, , k_low,  d_idx]
          se_high <- SE_on_grid[, , k_high, d_idx]
          SE_pair <- sqrt(se_low^2 + se_high^2)
          F_low  <- F_array[, , k_low,  d_idx]
          F_high <- F_array[, , k_high, d_idx]
          diff <- if (d_idx == 2L) F_low - F_high else F_high - F_low
          xi <- sqrt(n) * diff / pmax(SE_pair, se_floor)
          key <- paste(k_low, k_high, d_idx, sep = "_")
          xi_list[[key]] <- xi
          m <- max(xi)
          if (m > best) {
            best <- m
            idx_max <- which(xi == m, arr.ind = TRUE)[1, ]
            binding_local <- list(
              z_low = z_levels[k_low],
              z_high = z_levels[k_high],
              d = d_idx - 1L,
              x = x_grid[idx_max[1]],
              y = y_grid[idx_max[2]]
            )
          }
        }
      }
    }
    list(stat = max(best, 0),
         binding = binding_local,
         xi_all = xi_list)
  }

  observed <- compute_stat(F_on_grid)
  T_n <- observed$stat

  kappa_n <- if (adaptive) sqrt(log(log(max(n, 3)))) else Inf

  one_boot <- function() {
    W <- sample(c(-1, 1), n, replace = TRUE)
    F_star <- array(0, dim = c(XG, G, K, 2))
    for (k in seq_len(K)) {
      idx <- idx_by_z[[k]]
      Bk <- B[idx, , drop = FALSE]
      BtB_inv <- BtB_inv_by_k[[k]]
      Wk <- W[idx]
      for (d_idx in 1:2) {
        rk <- residuals_by[idx, , d_idx]
        WR <- Wk * rk
        delta_beta <- BtB_inv %*% (t(Bk) %*% WR)
        F_star[, , k, d_idx] <- B_at_xg %*% delta_beta
      }
    }
    # Studentised bootstrap statistic with adaptive selection.
    best <- 0
    if (K < 2L) return(0)
    for (k_low in seq_len(K - 1L)) {
      for (k_high in (k_low + 1L):K) {
        for (d_idx in 1:2) {
          se_low  <- SE_on_grid[, , k_low,  d_idx]
          se_high <- SE_on_grid[, , k_high, d_idx]
          SE_pair <- pmax(sqrt(se_low^2 + se_high^2), se_floor)
          diff_star <- if (d_idx == 2L) {
            F_star[, , k_low, d_idx] - F_star[, , k_high, d_idx]
          } else {
            F_star[, , k_high, d_idx] - F_star[, , k_low, d_idx]
          }
          xi_star <- sqrt(n) * diff_star / SE_pair
          # Adaptive selection: only include moments where observed xi > -kappa_n
          key <- paste(k_low, k_high, d_idx, sep = "_")
          obs_xi <- observed$xi_all[[key]]
          keep <- obs_xi >= -kappa_n
          vals <- xi_star[keep]
          if (length(vals) > 0L) {
            m <- max(vals)
            if (m > best) best <- m
          }
        }
      }
    }
    best
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

  p_value <- mean(boot_stats >= T_n)

  list(
    statistic = T_n,
    p_value = p_value,
    boot_stats = boot_stats,
    binding = observed$binding,
    n = n,
    kappa_n = kappa_n,
    basis_order = basis_order
  )
}
