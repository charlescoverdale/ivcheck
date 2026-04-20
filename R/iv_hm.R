#' Huber-Mellace (2015) test for instrument validity
#'
#' Tests the local exclusion and monotonicity assumptions via mean-based
#' moment inequalities on the implied complier outcome distribution.
#' Under the joint null, the principal-strata decomposition of the
#' observed conditional mean outcomes yields an implied mean outcome
#' for compliers that must lie inside the support of the outcome
#' variable. Rejection is evidence that at least one of exclusion or
#' monotonicity fails.
#'
#' @inheritParams iv_kitagawa
#' @param n_boot Number of multiplier-bootstrap replications.
#' @param alpha Significance level.
#' @param weights Optional non-negative numeric weights. Applied to the
#'   stratum-conditional means and to the multiplier bootstrap.
#' @param parallel Logical. Run bootstrap in parallel on POSIX.
#' @param ... Further arguments passed to methods.
#'
#' @return An object of class `iv_test` with elements:
#'   \item{test}{Character, `"Huber-Mellace (2015)"`.}
#'   \item{statistic}{Numeric test statistic, the max sqrt(n)-scaled
#'     positive-part violation of the implied-complier-mean bounds.}
#'   \item{p_value}{Bootstrap p-value.}
#'   \item{alpha}{Significance level used.}
#'   \item{implied_means}{Implied complier mean outcomes for `D = 0`
#'     and `D = 1`, with their bound-based admissibility range.}
#'   \item{boot_stats}{Vector of bootstrap test statistics.}
#'   \item{n}{Sample size.}
#'
#' @details
#' Huber and Mellace (2015) derive four inequality constraints on the
#' mean outcome of compliers that must hold under exclusion and
#' monotonicity. For binary treatment and binary instrument with
#' instrument levels ordered by first-stage propensity:
#' \deqn{E[Y \mid D = 1, Z = \mathrm{high}] \in
#'   \Bigl[ \frac{P_{AT}}{P_{AT}+P_{C}} \, y_{\min}
#'          + \frac{P_{C}}{P_{AT}+P_{C}} \, E[Y \mid D=1, Z=\mathrm{low}],
#'          \;
#'          \frac{P_{AT}}{P_{AT}+P_{C}} \, y_{\max}
#'          + \frac{P_{C}}{P_{AT}+P_{C}} \, E[Y \mid D=1, Z=\mathrm{low}]
#'   \Bigr],}
#' and analogously for the `D = 0` cells. The test statistic is the
#' maximum sqrt(n)-scaled positive-part violation across the four
#' inequalities. Critical values come from a multiplier bootstrap.
#'
#' `iv_hm` is complementary to `iv_kitagawa`. Kitagawa's distributional
#' inequalities on joint CDFs and Huber-Mellace's mean-based bounds
#' target different violations: Kitagawa has higher power against
#' shape-level failures; Huber-Mellace has higher power against
#' location-level failures and against violations that shift implied
#' complier means outside the observed outcome support.
#'
#' The package also implements related tests: `iv_kitagawa` for the CDF
#' inequalities, `iv_mw` for the covariate extension via intersection
#' bounds, and `iv_testjfe` for the judge-design case.
#'
#' @references
#' Huber, M. and Mellace, G. (2015). Testing Instrument Validity for
#' LATE Identification Based on Inequality Moment Constraints.
#' *Review of Economics and Statistics*, 97(2), 398-411.
#' \doi{10.1162/REST_a_00450}
#'
#' Imbens, G. W. and Angrist, J. D. (1994). Identification and Estimation
#' of Local Average Treatment Effects. *Econometrica*, 62(2), 467-475.
#' \doi{10.2307/2951620}
#'
#' @family iv_tests
#' @seealso [iv_kitagawa()] for the CDF-based form, [iv_mw()] for the
#'   conditional (covariate) extension, [iv_testjfe()] for judge designs,
#'   and [iv_check()] for a wrapper that runs all applicable tests.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 500
#' z <- sample(0:1, n, replace = TRUE)
#' d <- rbinom(n, 1, 0.3 + 0.4 * z)
#' y <- rnorm(n, mean = d)
#' iv_hm(y, d, z, n_boot = 200, parallel = FALSE)
#' }
#'
#' @export
iv_hm <- function(object, ...) {
  UseMethod("iv_hm")
}

#' @rdname iv_hm
#' @export
iv_hm.default <- function(object, d, z, n_boot = 1000, alpha = 0.05,
                          weights = NULL, parallel = TRUE, ...) {
  y <- object
  validate_numeric(y, "y")
  d_num <- validate_binary(d, "d")
  z_num <- validate_discrete(z, "z")
  n <- as.integer(check_lengths(y, d_num, z_num))

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

  z_vals <- sort(unique(z_num))
  if (length(z_vals) != 2L) {
    cli::cli_abort(
      "{.fn iv_hm} v0.1.0 supports binary instruments only; found {length(z_vals)} Z levels."
    )
  }

  # Order Z by first-stage E[D|Z] so z_high has higher compliance.
  p_by_z <- vapply(z_vals, function(zv) {
    idx <- which(z_num == zv)
    sum(w[idx] * d_num[idx]) / sum(w[idx])
  }, numeric(1))
  ord <- order(p_by_z)
  z_low_v  <- z_vals[ord[1]]
  z_high_v <- z_vals[ord[2]]
  p_low  <- p_by_z[ord[1]]
  p_high <- p_by_z[ord[2]]

  # Principal-strata shares under monotonicity:
  #   P(AT) = P(D=1 | Z_low)
  #   P(NT) = 1 - P(D=1 | Z_high)
  #   P(C)  = P(D=1 | Z_high) - P(D=1 | Z_low)
  P_AT <- p_low
  P_NT <- 1 - p_high
  P_C  <- p_high - p_low

  if (P_C <= 0) {
    cli::cli_abort(
      "The empirical first stage is zero or inverted (P(D=1|Z_high) <= P(D=1|Z_low)). HM test requires a positive first stage."
    )
  }

  # Cell means and outcome-support bounds.
  cell_mean <- function(d_val, z_val, y_vec = y, d_vec = d_num, z_vec = z_num,
                        w_vec = w) {
    idx <- which(d_vec == d_val & z_vec == z_val)
    if (length(idx) == 0L) return(NA_real_)
    sum(w_vec[idx] * y_vec[idx]) / sum(w_vec[idx])
  }

  mu_1_low  <- cell_mean(1, z_low_v)
  mu_1_high <- cell_mean(1, z_high_v)
  mu_0_low  <- cell_mean(0, z_low_v)
  mu_0_high <- cell_mean(0, z_high_v)

  y_min <- min(y)
  y_max <- max(y)

  # Four HM inequalities.
  # For d = 1: E[Y|D=1, Z_high] = (P_AT/(P_AT+P_C))*E[Y|AT] +
  #                                (P_C/(P_AT+P_C))*E[Y(1)|C].
  # Under exclusion E[Y|AT] = E[Y|D=1, Z_low]. Implied bound on
  # E[Y|D=1, Z_high] is inside the convex combination with E[Y(1)|C]
  # in [y_min, y_max].
  share_AT_d1 <- P_AT / (P_AT + P_C)
  share_C_d1  <- P_C  / (P_AT + P_C)
  lb_mu1_high <- share_AT_d1 * mu_1_low + share_C_d1 * y_min
  ub_mu1_high <- share_AT_d1 * mu_1_low + share_C_d1 * y_max

  # For d = 0: E[Y|D=0, Z_low] = (P_NT/(P_NT+P_C))*E[Y|NT] +
  #                              (P_C/(P_NT+P_C))*E[Y(0)|C].
  share_NT_d0 <- P_NT / (P_NT + P_C)
  share_C_d0  <- P_C  / (P_NT + P_C)
  lb_mu0_low <- share_NT_d0 * mu_0_high + share_C_d0 * y_min
  ub_mu0_low <- share_NT_d0 * mu_0_high + share_C_d0 * y_max

  # Violations: how far the observed cell means are outside the implied
  # bounds, scaled by sqrt(n).
  violations <- c(
    "d=1 upper" = max(0, mu_1_high - ub_mu1_high),
    "d=1 lower" = max(0, lb_mu1_high - mu_1_high),
    "d=0 upper" = max(0, mu_0_low  - ub_mu0_low ),
    "d=0 lower" = max(0, lb_mu0_low - mu_0_low)
  )
  obs_vec <- c(
    "d=1 upper" = mu_1_high - ub_mu1_high,
    "d=1 lower" = lb_mu1_high - mu_1_high,
    "d=0 upper" = mu_0_low  - ub_mu0_low,
    "d=0 lower" = lb_mu0_low - mu_0_low
  )

  T_n <- sqrt(n) * max(violations)
  worst_idx <- which.max(violations)
  binding <- list(
    direction = names(violations)[worst_idx],
    value = violations[worst_idx]
  )

  # Multiplier bootstrap. Recompute the four violations on a
  # Rademacher-weighted centred sample at each replication; use the
  # original plug-in shares for stratum weights (not bootstrap-perturbed)
  # to match Kitagawa-style inference.
  cell_idx <- list(
    d1_low  = which(d_num == 1 & z_num == z_low_v),
    d1_high = which(d_num == 1 & z_num == z_high_v),
    d0_low  = which(d_num == 0 & z_num == z_low_v),
    d0_high = which(d_num == 0 & z_num == z_high_v)
  )
  cell_w <- lapply(cell_idx, function(i) sum(w[i]))

  one_boot <- function() {
    W <- sample(c(-1, 1), n, replace = TRUE)
    # Bootstrap centred cell means: (1 / sum_w) * sum W_i * w_i * (y_i - mu_cell)
    mu_b <- function(cell_name, mu_cell) {
      idx <- cell_idx[[cell_name]]
      if (length(idx) == 0L) return(0)
      sum(W[idx] * w[idx] * (y[idx] - mu_cell)) / cell_w[[cell_name]]
    }
    bm1l <- mu_b("d1_low",  mu_1_low)
    bm1h <- mu_b("d1_high", mu_1_high)
    bm0l <- mu_b("d0_low",  mu_0_low)
    bm0h <- mu_b("d0_high", mu_0_high)
    v <- c(
      (mu_1_high + bm1h) - (share_AT_d1 * (mu_1_low + bm1l) + share_C_d1 * y_max) - (mu_1_high - ub_mu1_high),
      (share_AT_d1 * (mu_1_low + bm1l) + share_C_d1 * y_min) - (mu_1_high + bm1h) - (lb_mu1_high - mu_1_high),
      (mu_0_low + bm0l) - (share_NT_d0 * (mu_0_high + bm0h) + share_C_d0 * y_max) - (mu_0_low - ub_mu0_low),
      (share_NT_d0 * (mu_0_high + bm0h) + share_C_d0 * y_min) - (mu_0_low + bm0l) - (lb_mu0_low - mu_0_low)
    )
    sqrt(n) * max(pmax(v, 0))
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

  implied_means <- list(
    E_Y1_C_implied = (mu_1_high - share_AT_d1 * mu_1_low) / share_C_d1,
    E_Y0_C_implied = (mu_0_low  - share_NT_d0 * mu_0_high) / share_C_d0,
    y_support = c(y_min = y_min, y_max = y_max)
  )

  structure(
    list(
      test = "Huber-Mellace (2015)",
      statistic = T_n,
      p_value = p_value,
      alpha = alpha,
      n_boot = n_boot,
      boot_stats = boot_stats,
      binding = binding,
      violations = violations,
      observed_gaps = obs_vec,
      implied_means = implied_means,
      shares = c(AT = P_AT, NT = P_NT, C = P_C),
      n = n,
      call = sys.call()
    ),
    class = "iv_test"
  )
}

#' @rdname iv_hm
#' @export
iv_hm.fixest <- function(object, n_boot = 1000, alpha = 0.05,
                         weights = NULL, parallel = TRUE, ...) {
  yz <- extract_iv_data(object)
  iv_hm.default(
    object = yz$y, d = yz$d, z = yz$z,
    n_boot = n_boot, alpha = alpha,
    weights = weights, parallel = parallel, ...
  )
}

#' @rdname iv_hm
#' @export
iv_hm.ivreg <- function(object, n_boot = 1000, alpha = 0.05,
                        weights = NULL, parallel = TRUE, ...) {
  yz <- extract_iv_data(object)
  iv_hm.default(
    object = yz$y, d = yz$d, z = yz$z,
    n_boot = n_boot, alpha = alpha,
    weights = weights, parallel = parallel, ...
  )
}
