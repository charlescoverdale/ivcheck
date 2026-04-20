#' Frandsen-Lefgren-Leslie (2023) test for instrument validity in
#' judge-fixed-effects designs
#'
#' Jointly tests the local exclusion and monotonicity assumptions when
#' the instruments are a set of mutually exclusive dummy variables (the
#' leniency-of-assigned-judge design). Supports binary and multivalued
#' discrete treatments. Under the joint null, the per-judge mean
#' outcome `mu_j = E[Y | J = j]` must be a linear function of the
#' per-judge treatment propensities `P(D = d | J = j)`. Rejection is
#' evidence that at least one of exclusion or monotonicity fails.
#'
#' @inheritParams iv_kitagawa
#' @param z Factor, integer, or matrix of mutually exclusive dummy
#'   variables identifying the judge (or other random-assignment unit).
#' @param x Optional numeric vector, matrix, or data frame of covariates.
#'   If supplied, `y` and `d` are residualised on `x` before the per-
#'   judge means are computed.
#' @param method Reference distribution for the p-value. `"asymptotic"`
#'   (default) uses the chi-squared with `K - (basis_order + 1)`
#'   degrees of freedom. `"bootstrap"` uses the multiplier bootstrap
#'   of the restricted-model residual process. Asymptotic is fast and
#'   accurate for moderate `K`; bootstrap is preferred for small `K`
#'   or if errors are far from normal.
#' @param basis_order Order of the polynomial basis used to approximate
#'   the outcome / propensity function `phi(p)` in Frandsen-Lefgren-Leslie
#'   (2023) step 1. Default `1L` reduces to the Sargan-Hansen
#'   overidentification form, which imposes constant treatment effects.
#'   Values above 1 relax this to `phi(p) = delta_0 + delta_1 p +
#'   delta_2 p^2 + ... + delta_m p^m` and test the joint-zero
#'   restriction on judge residuals under the richer fit. Only binary
#'   treatment is supported when `basis_order > 1`. The slope-bounded
#'   moment-inequality component of the FLL test is not implemented in
#'   v0.1.0 (deferred to v0.2.0).
#'
#' @return An object of class `iv_test`; see [iv_kitagawa] for element
#'   descriptions. Additional elements:
#'   \item{n_judges}{Number of distinct judges / assignment groups.}
#'   \item{coef}{Fitted weighted-LS slope and intercept of `mu_j` on `p_j`.}
#'   \item{pairwise_late}{`K x K` matrix of pairwise Wald LATE estimates
#'     `(mu_j - mu_k) / (p_j - p_k)`. Under the null every entry
#'     estimates the common complier LATE.}
#'   \item{worst_pair}{List identifying the judge pair with the largest
#'     deviation of its Wald LATE from the fitted slope; useful for
#'     diagnosing the source of a rejection.}
#'
#' @details
#' Under the joint null, each pair of judges `(j, k)` identifies the
#' same complier LATE via the Wald estimator
#' `(mu_j - mu_k) / (p_j - p_k)`. The Frandsen-Lefgren-Leslie (2023)
#' test is the overidentification test of "all pairwise LATEs equal".
#' Under binary treatment with WLS weighting, that overidentification
#' test is algebraically the weighted sum of squared residuals from
#' the linear fit `mu_j = alpha + beta * p_j`, divided by a pooled
#' variance estimator. `iv_testjfe` computes this quadratic form and,
#' by default, compares to a chi-squared distribution with `K - 2`
#' degrees of freedom (the FLL asymptotic form). The multiplier
#' bootstrap of the restricted residual process is available via
#' `method = "bootstrap"` for small-K robustness.
#'
#' **Note on finite-sample size.** Per-judge propensities `p_j` enter
#' the test as estimated regressors. At modest per-judge sample sizes
#' (`n_j` below a few hundred), finite-sample binomial noise in
#' `hat p_j` compresses the distribution of the test statistic below
#' the asymptotic chi-squared reference, producing a test that is
#' mildly conservative at nominal 5 percent. Empirical size at
#' `K = 20`, `N = 3000` is 1.5 percent under the asymptotic method
#' and 2.5 percent under the bootstrap. Both methods sharpen toward
#' nominal as `n_j` grows. The bootstrap is recommended for
#' publication-grade p-values at modest `n_j`.
#'
#' The returned object includes `pairwise_late`, the `K x K` matrix of
#' pairwise Wald LATE estimates, and `worst_pair`, the judge pair with
#' the largest absolute deviation from the fitted slope. These are
#' diagnostic outputs in the sense of the paper's Figure 2: a pair
#' whose Wald LATE deviates far from the common slope is the first
#' place to look when investigating a rejection.
#'
#' Multivalued treatment is supported: for `D` with `M + 1` distinct
#' values (`0, 1, ..., M`), the fit becomes a multiple WLS regression
#' of `mu_j` on the `M`-vector `(P(D = 1 | J), ..., P(D = M | J))` and
#' the test statistic is compared to `chi^2_{K - M - 1}` (FLL 2023
#' section 4). `pairwise_late` and `worst_pair` are only defined for
#' binary `D` and return `NULL` otherwise.
#'
#' @references
#' Frandsen, B. R., Lefgren, L. J., and Leslie, E. C. (2023). Judging
#' Judge Fixed Effects. *American Economic Review*, 113(1), 253-277.
#' \doi{10.1257/aer.20201860}
#'
#' Imbens, G. W. and Angrist, J. D. (1994). Identification and Estimation
#' of Local Average Treatment Effects. *Econometrica*, 62(2), 467-475.
#' \doi{10.2307/2951620}
#'
#' @family iv_tests
#' @seealso [iv_kitagawa()] for the unconditional binary-treatment test,
#'   [iv_mw()] for the conditional version with covariates, and
#'   [iv_check()] for a one-shot wrapper that runs all applicable tests.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 2000
#' judge <- sample.int(20, n, replace = TRUE)
#' d <- rbinom(n, 1, 0.3 + 0.02 * judge)
#' y <- rnorm(n, mean = d)
#' iv_testjfe(y, d, judge, n_boot = 200, parallel = FALSE)
#' }
#'
#' @export
iv_testjfe <- function(object, ...) {
  UseMethod("iv_testjfe")
}

#' @rdname iv_testjfe
#' @export
iv_testjfe.default <- function(object, d, z, x = NULL, n_boot = 1000,
                               alpha = 0.05,
                               method = c("asymptotic", "bootstrap"),
                               weights = NULL,
                               basis_order = 1L,
                               parallel = TRUE, ...) {
  method <- match.arg(method)
  if (!is.numeric(basis_order) || length(basis_order) != 1L ||
      basis_order < 1L || basis_order != as.integer(basis_order)) {
    cli::cli_abort("{.arg basis_order} must be a positive integer.")
  }
  basis_order <- as.integer(basis_order)
  y <- object
  validate_numeric(y, "y")
  d_info <- validate_treatment_discrete(d, "d", max_levels = 20L)
  d_num <- d_info$d_num
  d_labels <- d_info$original_levels  # user-facing D values
  z_num <- validate_discrete(z, "z", max_levels = 500L)
  n <- as.integer(check_lengths(y, d_num, z_num))

  # Capture treatment levels BEFORE any residualisation on covariates.
  d_vals <- sort(unique(d_num))
  M <- length(d_vals) - 1L

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

  # Residualise y and the treatment-indicator columns on x if supplied.
  # Keep d_num (discrete level labels) intact so that level-based
  # bookkeeping still works; residualisation applies to the indicators
  # used in P_design, not the raw level labels.
  d_num_raw <- d_num
  if (!is.null(x)) {
    x_mat <- if (is.matrix(x) || is.data.frame(x)) as.matrix(x) else cbind(x)
    if (nrow(x_mat) != n) {
      cli::cli_abort("{.arg x} must have the same number of rows as {.arg y}.")
    }
    validate_numeric(as.vector(x_mat), "x")
    design <- cbind(1, x_mat)
    y <- as.numeric(stats::residuals(stats::lm.fit(design, y)))
    if (M == 1L) {
      d_num <- as.numeric(stats::residuals(stats::lm.fit(design, d_num)))
    }
    # For multivalued d: residualisation is applied inside P_design
    # construction (per-level). Do not overwrite d_num.
  }

  judges <- sort(unique(z_num))
  K <- length(judges)

  if (K < 3L) {
    cli::cli_abort(
      "{.fn iv_testjfe} requires at least three distinct judge levels; found {K}."
    )
  }
  if (K <= M + 1L) {
    cli::cli_abort(
      "{.fn iv_testjfe} needs K > M + 1 (K = {K} judges, M = {M} non-zero treatment levels) to have positive degrees of freedom."
    )
  }

  # Per-judge moments under weighting: n_j is effective (weighted) sample
  # size; mu_j and P_j use weighted means.
  n_j <- vapply(judges, function(j) sum(w[z_num == j]), numeric(1))
  mu_j <- vapply(judges, function(j) {
    idx <- which(z_num == j)
    sum(w[idx] * y[idx]) / sum(w[idx])
  }, numeric(1))
  if (basis_order > 1L && M > 1L) {
    cli::cli_abort(c(
      "{.arg basis_order > 1} is only supported for binary treatment in v0.1.0.",
      i = "Multivalued D with polynomial basis requires a tensor-product basis (planned v0.2.0)."
    ))
  }

  if (M == 1L) {
    p_j <- vapply(judges, function(j) {
      idx <- which(z_num == j)
      sum(w[idx] * d_num[idx]) / sum(w[idx])
    }, numeric(1))
    # Flexible basis (FLL 2023 equation 5 / 6): the fit component of the
    # FLL test generalises the linear WLS above to an arbitrary series
    # approximation phi(p) = delta_0 S_0(p) + ... + delta_m S_m(p). We
    # use polynomial bases here; b-splines could be added later. With
    # basis_order = 1 this reduces to the Sargan-Hansen form in the
    # paper's equation (6).
    if (basis_order == 1L) {
      P_design <- cbind(1, p_j)
    } else {
      P_design <- matrix(1, nrow = K, ncol = basis_order + 1L)
      for (p_idx in seq_len(basis_order)) {
        P_design[, p_idx + 1L] <- p_j^p_idx
      }
    }
    ind_resid_mat <- NULL  # binary path uses residualised d_num directly
  } else {
    P_design <- matrix(0, nrow = K, ncol = M + 1L)
    P_design[, 1L] <- 1
    # Retain the per-observation residualised indicators so the structural
    # residual u_i can use them (FWL). Using raw indicators against a
    # beta_hat estimated on residualised regressors inflates sigma^2 and
    # makes the chi^2 test conservative when covariates have within-judge
    # variation.
    ind_resid_mat <- matrix(0, nrow = n, ncol = M)
    for (m in seq_len(M)) {
      dv <- d_vals[m + 1L]
      ind <- as.numeric(d_num_raw == dv)
      if (!is.null(x)) {
        design <- cbind(1, x_mat)
        ind <- as.numeric(stats::residuals(stats::lm.fit(design, ind)))
      }
      ind_resid_mat[, m] <- ind
      P_design[, m + 1L] <- vapply(judges, function(j) {
        idx <- which(z_num == j)
        sum(w[idx] * ind[idx]) / sum(w[idx])
      }, numeric(1))
    }
  }

  # Weighted-LS fit mu_j ~ P_design with weights n_j.
  W <- diag(sqrt(n_j))
  lhs <- W %*% P_design
  rhs <- W %*% mu_j
  qr_fit <- qr(lhs)
  coef_fit <- qr.coef(qr_fit, rhs)
  alpha_hat <- coef_fit[1L]
  beta_vec <- coef_fit[-1L]
  mu_fit <- as.numeric(P_design %*% coef_fit)
  resid_j <- mu_j - mu_fit

  # Structural residuals u_i for pooled variance: remove the fitted
  # marginal effects of each treatment level. For binary D with the
  # flexible polynomial basis we evaluate phi at each observation's
  # estimated judge propensity and subtract; this is FLL (2023) step
  # 1's u_i = Y_i - phi_hat(p_hat(J_i)). For multivalued D with
  # covariates we use the per-observation residualised indicators
  # that P_design was built from, preserving the FWL identity between
  # the OLS fit and the residualised structural equation.
  if (M == 1L) {
    if (basis_order == 1L) {
      u_i <- y - alpha_hat - beta_vec * d_num
    } else {
      # Evaluate phi_hat(p_hat(J_i)) = sum_k delta_k * p_hat(J_i)^k at
      # each observation. d_num was residualised on x earlier (if any),
      # so using it as the "effective dose" matches the FWL identity.
      phi_at_i <- alpha_hat
      for (p_idx in seq_len(basis_order)) {
        phi_at_i <- phi_at_i + beta_vec[p_idx] * d_num^p_idx
      }
      u_i <- y - phi_at_i
    }
  } else if (is.null(ind_resid_mat)) {
    beta_full <- numeric(length(d_vals))
    beta_full[match(d_vals[-1L], d_vals)] <- beta_vec
    idx_map <- match(d_num_raw, d_vals)
    u_i <- y - alpha_hat - beta_full[idx_map]
  } else {
    u_i <- y - alpha_hat
    for (m in seq_len(M)) {
      u_i <- u_i - beta_vec[m] * ind_resid_mat[, m]
    }
  }
  ss_within <- 0
  for (j_idx in seq_len(K)) {
    idx <- which(z_num == judges[j_idx])
    u_sub <- u_i[idx]
    w_sub <- w[idx]
    u_mean <- sum(w_sub * u_sub) / sum(w_sub)
    ss_within <- ss_within + sum(w_sub * (u_sub - u_mean)^2)
  }
  df_within <- sum(w) - K
  sigma2_hat <- if (df_within > 0) ss_within / df_within else NA_real_

  T_n <- if (is.finite(sigma2_hat) && sigma2_hat > 0) {
    sum(n_j * resid_j^2) / sigma2_hat
  } else {
    NA_real_
  }

  df_test <- if (M == 1L && basis_order > 1L) {
    K - (basis_order + 1L)
  } else {
    K - (M + 1L)
  }
  p_value_asy <- if (is.finite(T_n) && df_test > 0L) {
    1 - stats::pchisq(T_n, df = df_test)
  } else {
    NA_real_
  }

  # Legacy scalar p_j is used by the binary pairwise diagnostics below.
  # For multivalued D and for basis_order > 1 the concept of a single
  # LATE "slope" does not apply; diagnostics are disabled in those cases.
  if (M == 1L && basis_order == 1L) {
    beta_hat <- beta_vec
  } else {
    beta_hat <- NA_real_
    if (M > 1L) p_j <- rep(NA_real_, K)
  }

  # Multiplier bootstrap over per-judge (mu, P_design) deviations. Works
  # for both the binary and multivalued case.
  one_boot <- function() {
    Wv <- sample(c(-1, 1), n, replace = TRUE)
    mu_star <- vapply(seq_len(K), function(j_idx) {
      idx <- which(z_num == judges[j_idx])
      sum(Wv[idx] * (y[idx] - mu_j[j_idx])) / length(idx)
    }, numeric(1))
    # P_star: bootstrap analogue of per-judge propensity design matrix.
    # Perturb each column of P_design by the corresponding centered
    # indicator projection through the multiplier weights.
    P_star <- matrix(0, nrow = K, ncol = ncol(P_design))
    P_star[, 1L] <- 0  # intercept column has no noise
    if (M == 1L) {
      delta_p <- vapply(seq_len(K), function(j_idx) {
        idx <- which(z_num == judges[j_idx])
        sum(Wv[idx] * (d_num[idx] - p_j[j_idx])) / length(idx)
      }, numeric(1))
      if (basis_order == 1L) {
        P_star[, 2L] <- delta_p
      } else {
        # Delta-method for polynomial bases: d(p^k)/dp = k * p^{k-1}.
        for (p_idx in seq_len(basis_order)) {
          P_star[, p_idx + 1L] <- p_idx * p_j^(p_idx - 1L) * delta_p
        }
      }
    } else {
      for (m in seq_len(M)) {
        dv <- d_vals[m + 1L]
        col_m <- P_design[, m + 1L]
        P_star[, m + 1L] <- vapply(seq_len(K), function(j_idx) {
          idx <- which(z_num == judges[j_idx])
          sum(Wv[idx] * ((d_num_raw[idx] == dv) - col_m[j_idx])) / length(idx)
        }, numeric(1))
      }
    }
    # Bootstrap WLS fit on (P_design + P_star, mu + mu_star)
    P_b <- P_design + P_star
    mu_b <- mu_j + mu_star
    Wdiag <- diag(sqrt(n_j))
    coef_b <- tryCatch(
      qr.coef(qr(Wdiag %*% P_b), Wdiag %*% mu_b),
      error = function(e) rep(NA_real_, ncol(P_b))
    )
    if (any(!is.finite(coef_b))) return(NA_real_)
    fit_b <- as.numeric(P_b %*% coef_b)
    resid_b <- mu_star - (fit_b - mu_fit)
    if (is.finite(sigma2_hat) && sigma2_hat > 0) {
      sum(n_j * resid_b^2) / sigma2_hat
    } else {
      NA_real_
    }
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

  p_value_boot <- if (length(boot_stats) > 0L && is.finite(T_n)) {
    mean(boot_stats >= T_n, na.rm = TRUE)
  } else {
    NA_real_
  }

  p_value <- switch(method,
    asymptotic = p_value_asy,
    bootstrap = p_value_boot
  )

  binding_j <- if (all(is.finite(resid_j))) {
    idx <- which.max(abs(resid_j))
    if (M == 1L) {
      list(judge = judges[idx], mu = mu_j[idx], p = p_j[idx],
           residual = resid_j[idx])
    } else {
      list(judge = judges[idx], mu = mu_j[idx],
           residual = resid_j[idx])
    }
  } else {
    NULL
  }

  # Pairwise Wald LATE matrix (binary D only). Under binary treatment,
  # each pair of judges identifies the common complier LATE via
  # (mu_j - mu_k) / (p_j - p_k). Not defined for multivalued D; NULL is
  # returned in that case.
  if (M == 1L && basis_order == 1L) {
    pairwise_late <- matrix(NA_real_, nrow = K, ncol = K,
                            dimnames = list(judges, judges))
    for (i in seq_len(K)) {
      for (k2 in seq_len(K)) {
        if (i == k2) next
        dp <- p_j[i] - p_j[k2]
        if (abs(dp) > 1e-8) {
          pairwise_late[i, k2] <- (mu_j[i] - mu_j[k2]) / dp
        }
      }
    }
    late_dev <- pairwise_late - beta_hat
    late_dev[!is.finite(late_dev)] <- NA_real_
    dp_mat <- outer(p_j, p_j, "-")
    late_dev[abs(dp_mat) < 0.01] <- NA_real_
    worst_pair <- if (any(!is.na(late_dev))) {
      pos <- which.max(abs(late_dev))
      r <- ((pos - 1L) %% K) + 1L
      c <- ((pos - 1L) %/% K) + 1L
      list(
        judge_j = judges[r], judge_k = judges[c],
        late_jk = pairwise_late[r, c],
        deviation_from_beta = late_dev[r, c]
      )
    } else {
      NULL
    }
  } else {
    pairwise_late <- NULL
    worst_pair <- NULL
  }

  structure(
    list(
      test = "Frandsen-Lefgren-Leslie (2023)",
      statistic = T_n,
      p_value = p_value,
      p_value_asymptotic = p_value_asy,
      p_value_bootstrap = p_value_boot,
      method = method,
      alpha = alpha,
      n_boot = n_boot,
      boot_stats = boot_stats,
      binding = binding_j,
      pairwise_late = pairwise_late,
      worst_pair = worst_pair,
      coef = if (M == 1L && basis_order == 1L) {
               c(intercept = alpha_hat, slope = beta_hat)
             } else if (M == 1L) {
               c(intercept = alpha_hat,
                 stats::setNames(beta_vec, paste0("delta_p", seq_len(basis_order))))
             } else {
               c(intercept = alpha_hat,
                 stats::setNames(beta_vec, paste0("beta_d", d_labels[-1L])))
             },
      n_judges = K,
      n_treatment_levels = M + 1L,
      basis_order = basis_order,
      n = n,
      call = sys.call()
    ),
    class = "iv_test"
  )
}

#' @rdname iv_testjfe
#' @export
iv_testjfe.fixest <- function(object, x = NULL, n_boot = 1000,
                              alpha = 0.05,
                              method = c("asymptotic", "bootstrap"),
                              weights = NULL,
                              basis_order = 1L,
                              parallel = TRUE, ...) {
  method <- match.arg(method)
  yz <- extract_iv_data(object)
  iv_testjfe.default(
    object = yz$y, d = yz$d, z = yz$z, x = x %||% yz$x,
    n_boot = n_boot, alpha = alpha, method = method,
    weights = weights, basis_order = basis_order,
    parallel = parallel, ...
  )
}

#' @rdname iv_testjfe
#' @export
iv_testjfe.ivreg <- function(object, x = NULL, n_boot = 1000,
                             alpha = 0.05,
                             method = c("asymptotic", "bootstrap"),
                             weights = NULL,
                             basis_order = 1L,
                             parallel = TRUE, ...) {
  method <- match.arg(method)
  yz <- extract_iv_data(object)
  iv_testjfe.default(
    object = yz$y, d = yz$d, z = yz$z, x = x %||% yz$x,
    n_boot = n_boot, alpha = alpha, method = method,
    weights = weights, basis_order = basis_order,
    parallel = parallel, ...
  )
}
