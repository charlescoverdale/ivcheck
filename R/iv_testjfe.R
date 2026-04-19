#' Frandsen-Lefgren-Leslie (2023) test for instrument validity in
#' judge-fixed-effects designs
#'
#' Jointly tests the local exclusion and monotonicity assumptions when
#' treatment is binary and the instruments are a set of mutually
#' exclusive dummy variables (the leniency-of-assigned-judge design).
#' Under the joint null, the per-judge mean outcome `mu_j = E[Y | J = j]`
#' must be a linear function of the per-judge treatment propensity
#' `p_j = E[D | J = j]`. Rejection is evidence that at least one of
#' exclusion or monotonicity fails.
#'
#' @inheritParams iv_kitagawa
#' @param z Factor, integer, or matrix of mutually exclusive dummy
#'   variables identifying the judge (or other random-assignment unit).
#' @param x Optional numeric vector, matrix, or data frame of covariates.
#'   If supplied, `y` and `d` are residualised on `x` before the per-
#'   judge means are computed.
#'
#' @return An object of class `iv_test`; see [iv_kitagawa] for element
#'   descriptions. Additional elements:
#'   \item{n_judges}{Number of distinct judges / assignment groups.}
#'   \item{coef}{Fitted weighted-LS slope and intercept of `mu_j` on `p_j`.}
#'
#' @details
#' This is a simplified first-pass implementation. It tests the linearity
#' implication of the FLL (2023) joint null using a weighted-LS fit of
#' `mu_j` on `p_j` and a chi-squared asymptotic distribution for the
#' weighted sum of squared residuals (degrees of freedom equal to the
#' number of judges minus two). The full nonparametric restricted-LS
#' test of Frandsen, Lefgren, and Leslie (2023), which handles the
#' case of many small judge cells and generalises to multivalued
#' treatment, is planned for v0.2.0. Users needing the published test
#' today should consult the Stata `testjfe` module (Frandsen, BYU,
#' 2020).
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
                               alpha = 0.05, weights = NULL,
                               parallel = TRUE, ...) {
  y <- object
  validate_numeric(y, "y")
  d_num <- validate_binary(d, "d")
  z_num <- validate_discrete(z, "z", max_levels = 500L)
  n <- check_lengths(y, d_num, z_num)

  if (!is.null(weights)) {
    cli::cli_warn(
      "The {.arg weights} argument is not yet implemented and is ignored in v0.1.0."
    )
  }

  # Residualise y and d on x if supplied
  if (!is.null(x)) {
    x_mat <- if (is.matrix(x) || is.data.frame(x)) as.matrix(x) else cbind(x)
    if (nrow(x_mat) != n) {
      cli::cli_abort("{.arg x} must have the same number of rows as {.arg y}.")
    }
    validate_numeric(as.vector(x_mat), "x")
    design <- cbind(1, x_mat)
    y <- as.numeric(stats::residuals(stats::lm.fit(design, y)))
    d_num <- as.numeric(stats::residuals(stats::lm.fit(design, d_num)))
  }

  judges <- sort(unique(z_num))
  K <- length(judges)
  if (K < 3L) {
    cli::cli_abort(
      "{.fn iv_testjfe} requires at least three distinct judge levels; found {K}."
    )
  }

  # Per-judge moments
  n_j <- vapply(judges, function(j) sum(z_num == j), integer(1))
  p_j <- vapply(judges, function(j) mean(d_num[z_num == j]), numeric(1))
  mu_j <- vapply(judges, function(j) mean(y[z_num == j]), numeric(1))

  # Weighted-LS fit mu_j ~ 1 + p_j with weights n_j
  # (minimise sum_j n_j * (mu_j - alpha - beta * p_j)^2)
  sw <- sum(n_j)
  p_bar <- sum(n_j * p_j) / sw
  mu_bar <- sum(n_j * mu_j) / sw
  s_pp <- sum(n_j * (p_j - p_bar)^2)
  s_pm <- sum(n_j * (p_j - p_bar) * (mu_j - mu_bar))
  beta_hat <- if (s_pp > 0) s_pm / s_pp else 0
  alpha_hat <- mu_bar - beta_hat * p_bar
  mu_fit <- alpha_hat + beta_hat * p_j
  resid_j <- mu_j - mu_fit

  # Pooled within-judge variance estimator, using the structural residuals
  # u_i = y_i - alpha_hat - beta_hat * d_i. This removes the first-stage
  # binomial contribution of D, so sigma2_hat matches the asymptotic
  # variance of (mu_j_hat - alpha_0 - beta_0 * p_j_hat). Under the FLL
  # null, the resulting T_n has chi^2_{K-2} distribution.
  u_i <- y - alpha_hat - beta_hat * d_num
  ss_within <- 0
  for (j_idx in seq_len(K)) {
    idx <- which(z_num == judges[j_idx])
    u_sub <- u_i[idx]
    ss_within <- ss_within + sum((u_sub - mean(u_sub))^2)
  }
  df_within <- n - K
  sigma2_hat <- if (df_within > 0) ss_within / df_within else NA_real_

  T_n <- if (is.finite(sigma2_hat) && sigma2_hat > 0) {
    sum(n_j * resid_j^2) / sigma2_hat
  } else {
    NA_real_
  }

  df_test <- K - 2L
  p_value <- if (is.finite(T_n)) 1 - stats::pchisq(T_n, df = df_test) else NA_real_

  # A bootstrap analogue of the test statistic (multiplier on within-
  # judge deviations), kept for consistency with the other tests and
  # for an alternative p-value when large-K asymptotics are suspect.
  one_boot <- function() {
    W <- sample(c(-1, 1), n, replace = TRUE)
    mu_star <- vapply(seq_len(K), function(j_idx) {
      idx <- which(z_num == judges[j_idx])
      sum(W[idx] * (y[idx] - mu_j[j_idx])) / length(idx)
    }, numeric(1))
    p_star <- vapply(seq_len(K), function(j_idx) {
      idx <- which(z_num == judges[j_idx])
      sum(W[idx] * (d_num[idx] - p_j[j_idx])) / length(idx)
    }, numeric(1))
    # Bootstrap fit: regress (mu_j + mu_star) on (p_j + p_star) with weights
    pb <- p_j + p_star
    mb <- mu_j + mu_star
    pbar_b <- sum(n_j * pb) / sw
    mbar_b <- sum(n_j * mb) / sw
    s_pp_b <- sum(n_j * (pb - pbar_b)^2)
    s_pm_b <- sum(n_j * (pb - pbar_b) * (mb - mbar_b))
    beta_b <- if (s_pp_b > 0) s_pm_b / s_pp_b else 0
    alpha_b <- mbar_b - beta_b * pbar_b
    fit_b <- alpha_b + beta_b * pb
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

  binding_j <- if (all(is.finite(resid_j))) {
    idx <- which.max(abs(resid_j))
    list(judge = judges[idx], mu = mu_j[idx], p = p_j[idx], residual = resid_j[idx])
  } else {
    NULL
  }

  structure(
    list(
      test = "Frandsen-Lefgren-Leslie (2023)",
      statistic = T_n,
      p_value = p_value,
      alpha = alpha,
      n_boot = n_boot,
      boot_stats = boot_stats,
      binding = binding_j,
      coef = c(intercept = alpha_hat, slope = beta_hat),
      n_judges = K,
      n = n,
      call = sys.call()
    ),
    class = "iv_test"
  )
}

#' @rdname iv_testjfe
#' @export
iv_testjfe.fixest <- function(object, x = NULL, n_boot = 1000,
                              alpha = 0.05, weights = NULL,
                              parallel = TRUE, ...) {
  yz <- extract_iv_data(object)
  iv_testjfe.default(
    object = yz$y, d = yz$d, z = yz$z, x = x %||% yz$x,
    n_boot = n_boot, alpha = alpha,
    weights = weights, parallel = parallel, ...
  )
}

#' @rdname iv_testjfe
#' @export
iv_testjfe.ivreg <- function(object, x = NULL, n_boot = 1000,
                             alpha = 0.05, weights = NULL,
                             parallel = TRUE, ...) {
  yz <- extract_iv_data(object)
  iv_testjfe.default(
    object = yz$y, d = yz$d, z = yz$z, x = x %||% yz$x,
    n_boot = n_boot, alpha = alpha,
    weights = weights, parallel = parallel, ...
  )
}
