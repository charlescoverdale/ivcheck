# Full CLR-style conditional moment inequality test for iv_mw with
# covariates. Replaces the bin-stratified simplification of v0.1.0-pre
# with series-regression conditional CDF estimation, plug-in
# heteroscedasticity-robust standard errors, and adaptive moment
# selection in the style of Andrews and Soares (2010).
#
# Not exported.

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
