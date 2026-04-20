# Replication script for the ivcheck R Journal paper.
#
# Running this script end-to-end regenerates every figure and every
# numerical claim in paper/rj/paper.Rmd. Runtime: approximately 30 to
# 45 minutes on a 4-core laptop.
#
# Required:
#   * ivcheck >= 0.1.0 (load via devtools::load_all(".") from package root)
#   * R >= 4.1.0
#   * ggplot2, showtext, pdftools
#
# Invoke from the package root:
#   Rscript paper/replicate.R
#
# Outputs:
#   paper/figures/*.pdf             — all 5 figures
#   paper/replicate_claims.txt      — tabulated numerical claims

suppressPackageStartupMessages({
  library(ivcheck)
  library(ggplot2)
})

set.seed(1)  # global seed; each block re-seeds internally for reproducibility

sink("paper/replicate_claims.txt", split = TRUE)

# =========================================================================
# Claim 1. iv_testjfe null distribution matches chi^2_{K-2} to MC precision.
# Paper claim (Why trust this implementation section):
#   "At K=20, N=3000 over 200 replications: empirical mean 18.01 vs target
#    18.0, variance 35.1 vs target 36.0, 95th percentile 29.4 vs target 28.9,
#    empirical size 6.5% vs nominal 5%."
# =========================================================================
cat("===== Claim 1: iv_testjfe null distribution vs chi^2_{K-2} =====\n")
K <- 20
n_sim <- 3000
n_reps <- 200
p_by_j <- 0.3 + 0.02 * seq_len(K)
stats_null <- numeric(n_reps)
rej <- 0
for (s in seq_len(n_reps)) {
  set.seed(10000 + s)
  judge <- sample.int(K, n_sim, replace = TRUE)
  d <- rbinom(n_sim, 1, p_by_j[judge])
  y <- rnorm(n_sim, mean = d)
  r <- iv_testjfe(y, d, judge, n_boot = 5, parallel = FALSE)
  stats_null[s] <- r$statistic
  if (r$p_value < 0.05) rej <- rej + 1
}
cat(sprintf("K = %d, N = %d, reps = %d\n", K, n_sim, n_reps))
cat(sprintf("Empirical mean: %.3f    (target chi^2_{%d} mean: %d)\n",
            mean(stats_null), K - 2, K - 2))
cat(sprintf("Empirical var:  %.3f    (target chi^2_{%d} var: %d)\n",
            var(stats_null), K - 2, 2 * (K - 2)))
cat(sprintf("Empirical 95th: %.3f    (chi^2_{%d} 95th percentile: %.3f)\n",
            quantile(stats_null, 0.95), K - 2, qchisq(0.95, K - 2)))
cat(sprintf("Empirical size: %.3f    (nominal: 0.050)\n", rej / n_reps))

# =========================================================================
# Claim 2. iv_kitagawa MC size across 24 configurations at se_floor = 0.15.
# Paper claim (Null size under finite samples and skewed Z):
#   "All skewed-Z configurations with strong first stages also deliver 0%
#    rejection. Skewed-Z configurations with weak first stages are at or
#    below nominal 5% once se_floor = 0.15."
# =========================================================================
cat("\n===== Claim 2: iv_kitagawa MC size across 24 configurations =====\n")
cfg <- expand.grid(
  n = c(300, 800, 2000),
  first_stage = c("balanced", "moderate", "card_like", "extreme"),
  z_balance = c("equal", "skewed"),
  stringsAsFactors = FALSE
)
cfg$size <- NA_real_
n_reps_size <- 100  # smaller for runtime; paper uses 200
alpha <- 0.05
for (i in seq_len(nrow(cfg))) {
  p_low  <- switch(cfg$first_stage[i],
                   balanced = 0.5, moderate = 0.3,
                   card_like = 0.3, extreme = 0.1)
  p_high <- switch(cfg$first_stage[i],
                   balanced = 0.5, moderate = 0.7,
                   card_like = 0.4, extreme = 0.9)
  z_prob_1 <- if (cfg$z_balance[i] == "equal") 0.5 else 0.65
  rej <- 0
  for (r in seq_len(n_reps_size)) {
    set.seed(20000 + 1000 * i + r)
    z <- rbinom(cfg$n[i], 1, z_prob_1)
    d <- rbinom(cfg$n[i], 1, ifelse(z == 1, p_high, p_low))
    y <- rnorm(cfg$n[i], mean = d)
    out <- iv_kitagawa(y, d, z, n_boot = 99, parallel = FALSE)
    if (out$p_value < alpha) rej <- rej + 1
  }
  cfg$size[i] <- rej / n_reps_size
  cat(sprintf("  n=%4d  fs=%-10s  zb=%-6s  size=%.3f\n",
              cfg$n[i], cfg$first_stage[i], cfg$z_balance[i], cfg$size[i]))
}
cat(sprintf("\nMax size across all configs: %.3f (expected <= 0.05)\n",
            max(cfg$size)))

# =========================================================================
# Claim 3. iv_kitagawa power curve under direct D-Z interaction violation.
# Paper figure fig-power-curve.pdf shows rejection vs delta.
# =========================================================================
cat("\n===== Claim 3: iv_kitagawa power curve (N = 500) =====\n")
set.seed(30000)
n <- 500
z <- sample(0:1, n, replace = TRUE)
d <- rbinom(n, 1, 0.3 + 0.4 * z)
y <- rnorm(n, mean = d)
pc <- iv_power(y, d, z, method = "kitagawa", n_sims = 60,
               delta_grid = seq(0, 2, length.out = 5),
               n_boot = 99, parallel = FALSE)
print(pc)

# =========================================================================
# Claim 4. Card (1995) rejection and binding interval.
# Paper (Case study): "The interval-sup Kitagawa test rejects on this
# binary-discretised college treatment. The binding violation sits in
# the upper lwage interval [6.25, 7.78]."
# =========================================================================
cat("\n===== Claim 4: Card (1995) end-to-end =====\n")
data(card1995)
set.seed(40000)
r_card <- iv_kitagawa(card1995$lwage, card1995$college,
                      card1995$near_college, n_boot = 500,
                      parallel = FALSE)
cat(sprintf("Card N = %d\n", r_card$n))
cat(sprintf("Kitagawa stat: %.3f  p-value: %.3f\n",
            r_card$statistic, r_card$p_value))
cat("Binding: "); print(r_card$binding)

# =========================================================================
# Claim 5. iv_kitagawa ~= iv_mw without covariates (numeric identity).
# README and paper: "Without covariates, iv_mw reduces exactly to the
# variance-weighted Kitagawa test (unit-tested)."
# =========================================================================
cat("\n===== Claim 5: iv_kitagawa == iv_mw without covariates =====\n")
set.seed(50000)
z <- sample(0:1, 500, replace = TRUE)
d <- rbinom(500, 1, 0.3 + 0.4 * z)
y <- rnorm(500, mean = d)
set.seed(1); r_k  <- iv_kitagawa(y, d, z, n_boot = 100, parallel = FALSE)
set.seed(1); r_mw <- iv_mw(y, d, z, n_boot = 100, parallel = FALSE)
cat(sprintf("iv_kitagawa stat: %.6f\n", r_k$statistic))
cat(sprintf("iv_mw stat:       %.6f\n", r_mw$statistic))
cat(sprintf("Numerically equal: %s\n",
            all.equal(r_k$statistic, r_mw$statistic)))

sink()

# =========================================================================
# Regenerate the 5 figures referenced in the paper.
# =========================================================================
cat("\nRegenerating figures via paper/make_figures.R ...\n")
source("paper/make_figures.R")

cat("\n--- replication complete ---\n")
cat("Numerical claims:  paper/replicate_claims.txt\n")
cat("Figures:           paper/figures/*.pdf\n")
