# Size-study generator for the ivcheck R Journal paper.
#
# Produces the empirical-size numbers quoted in the "Monte Carlo
# validation" section, and saves them as reproducible artifacts:
#   paper/size_stats.rds      iv_kitagawa size, 24 configs (se_floor = 0.15)
#   paper/mw_size_stats.rds   iv_mw conditional-path size (valid-IV-with-X null)
#
# Runtime: ~7-9 minutes across 6 cores. Invoke from the package root:
#   Rscript paper/make_size_study.R
#
# Determinism: every replication sets its own seed, so results are
# reproducible and independent of the number of cores.

suppressPackageStartupMessages({
  # Use the development source when run from the package root; fall back
  # to the installed package otherwise.
  if (requireNamespace("devtools", quietly = TRUE) && file.exists("R/iv_kitagawa.R")) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(ivcheck)
  }
  library(parallel)
})

N_REPS    <- 500       # MC replications per cell
N_BOOT    <- 199       # bootstrap replications per test
ALPHA     <- 0.05
SE_FLOOR  <- 0.15
N_CORES   <- max(1L, min(6L, detectCores() - 2L))

mc_se <- function(p, reps) sqrt(p * (1 - p) / reps) * 100

# ---------------------------------------------------------------------------
# Part A. iv_kitagawa size across 24 configurations at se_floor = 0.15.
# Valid-IV null: Y = D + N(0,1); Z shifts the first stage only.
# ---------------------------------------------------------------------------
cfg <- expand.grid(
  n           = c(300, 800, 2000),
  first_stage = c("balanced", "moderate", "card_like", "extreme"),
  z_balance   = c("equal", "skewed"),
  stringsAsFactors = FALSE
)
p_low_of  <- c(balanced = 0.5, moderate = 0.3, card_like = 0.3, extreme = 0.1)
p_high_of <- c(balanced = 0.5, moderate = 0.7, card_like = 0.4, extreme = 0.9)

kit_one <- function(i) {
  n      <- cfg$n[i]
  p_low  <- p_low_of[[cfg$first_stage[i]]]
  p_high <- p_high_of[[cfg$first_stage[i]]]
  zp1    <- if (cfg$z_balance[i] == "equal") 0.5 else 0.65
  rej <- 0L
  for (r in seq_len(N_REPS)) {
    set.seed(20000 + 1000 * i + r)
    z <- rbinom(n, 1, zp1)
    d <- rbinom(n, 1, ifelse(z == 1, p_high, p_low))
    y <- rnorm(n, mean = d)
    out <- tryCatch(
      iv_kitagawa(y, d, z, n_boot = N_BOOT, se_floor = SE_FLOOR, parallel = FALSE),
      error = function(e) NULL)
    if (!is.null(out) && is.finite(out$p_value) && out$p_value < ALPHA) rej <- rej + 1L
  }
  rej / N_REPS
}

cat(sprintf("Part A: iv_kitagawa size, %d configs x %d reps on %d cores...\n",
            nrow(cfg), N_REPS, N_CORES))
cfg$size <- unlist(mclapply(seq_len(nrow(cfg)), kit_one, mc.cores = N_CORES))
for (i in seq_len(nrow(cfg)))
  cat(sprintf("  n=%4d  fs=%-10s  zb=%-6s  size=%.3f\n",
              cfg$n[i], cfg$first_stage[i], cfg$z_balance[i], cfg$size[i]))
cat(sprintf("Max size: %.3f ; cells > nominal 0.05: %d of %d ; MC se ~%.2f pp\n",
            max(cfg$size), sum(cfg$size > ALPHA), nrow(cfg), mc_se(0.05, N_REPS)))

size_stats <- list(
  table = cfg, n_reps = N_REPS, n_boot = N_BOOT, se_floor = SE_FLOOR,
  alpha = ALPHA, mc_se_pp = mc_se(0.05, N_REPS),
  n_above_nominal = sum(cfg$size > ALPHA), max_size = max(cfg$size))
saveRDS(size_stats, "paper/size_stats.rds")

# ---------------------------------------------------------------------------
# Part B. iv_mw conditional-path size (valid IV conditional on one covariate).
# Y = 0.8 D + 0.5 X + N(0,1); D propensity monotone in Z; no direct Z->Y.
# This is the only validation of the conditional CLR series-regression path.
# ---------------------------------------------------------------------------
mw_one <- function(n, r) {
  set.seed(40000 + 1000 * n + r)
  x <- rnorm(n)
  z <- rbinom(n, 1, 0.5)
  d <- rbinom(n, 1, plogis(-0.3 + 1.0 * z + 0.5 * x))
  y <- rnorm(n, mean = 0.8 * d + 0.5 * x)
  out <- tryCatch(iv_mw(y, d, z, x = x, n_boot = N_BOOT, parallel = FALSE),
                  error = function(e) NULL)
  if (is.null(out) || !is.finite(out$p_value)) NA else as.integer(out$p_value < ALPHA)
}

cat("\nPart B: iv_mw conditional-path size (current kappa_n = sqrt(log log n))...\n")
mw_rows <- list()
for (n in c(800, 2000)) {
  rej <- unlist(mclapply(seq_len(N_REPS), function(r) mw_one(n, r), mc.cores = N_CORES))
  sz  <- mean(rej, na.rm = TRUE)
  cat(sprintf("  n=%4d  size=%.3f  (n_eff=%d, MC se ~%.2f pp)\n",
              n, sz, sum(!is.na(rej)), mc_se(0.05, sum(!is.na(rej)))))
  mw_rows[[length(mw_rows) + 1L]] <- data.frame(n = n, size = sz, n_eff = sum(!is.na(rej)))
}
mw_size_stats <- list(
  table = do.call(rbind, mw_rows), n_reps = N_REPS, n_boot = N_BOOT, alpha = ALPHA,
  kappa_rule = "sqrt(log(log(n)))")
saveRDS(mw_size_stats, "paper/mw_size_stats.rds")

cat("\n--- size study complete ---\n")
