# Figure generator for the ivcheck R Journal paper.
#
# Builds five figures into paper/figures/:
#   fig-card-bootstrap.pdf   Kitagawa bootstrap distribution on Card 1995
#   fig-power-curve.pdf      Kitagawa power curve under increasing violation
#   fig-testjfe-null.pdf     Monte Carlo null vs chi^2 for iv_testjfe
#   fig-judges-scatter.pdf   mu_j vs p_j (valid vs violation)
#   fig-pairwise-late.pdf    Pairwise LATE heatmap for judge-IV violation
#
# Run from the package root after devtools::load_all(".") or with the
# installed version on the library path.

suppressPackageStartupMessages({
  library(ivcheck)
  library(ggplot2)
  library(showtext)
})

font_add("HelveticaNeue",
         regular = "/System/Library/Fonts/Helvetica.ttc",
         bold    = "/System/Library/Fonts/Helvetica.ttc",
         italic  = "/System/Library/Fonts/Helvetica.ttc")
showtext_auto()
showtext_opts(dpi = 300)

fig_dir <- "paper/figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

ok_blue   <- "#0072B2"
ok_orange <- "#E69F00"
ok_green  <- "#009E73"
ok_red    <- "#D55E00"
ok_purple <- "#CC79A7"
ok_sky    <- "#56B4E9"

fam <- "HelveticaNeue"

theme_wp <- function(base_size = 10) {
  theme_bw(base_size = base_size, base_family = fam) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      panel.border = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey85"),
      axis.line = element_line(linewidth = 0.35, colour = "grey25"),
      axis.ticks = element_line(linewidth = 0.35, colour = "grey25"),
      axis.ticks.length = unit(2.5, "pt"),
      axis.text = element_text(size = base_size, colour = "grey20"),
      axis.title = element_text(size = base_size, colour = "grey20"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1, family = fam),
      legend.key.height = unit(10, "pt"),
      legend.key.width = unit(22, "pt"),
      plot.margin = margin(6, 10, 6, 6)
    )
}

theme_wp_hmap <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = fam) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_text(size = base_size, colour = "grey20", family = fam),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = base_size - 1, family = fam),
      legend.text = element_text(size = base_size - 2, family = fam),
      plot.margin = margin(6, 10, 6, 6)
    )
}

# -----------------------------------------------------------------------------
# Figure 1: Kitagawa bootstrap distribution on Card (1995) data.
# -----------------------------------------------------------------------------
set.seed(1)
data(card1995, envir = environment())
k_card <- iv_kitagawa(
  card1995$lwage,
  card1995$college,
  card1995$near_college,
  n_boot = 1000,
  parallel = FALSE
)

df_boot <- data.frame(stat = k_card$boot_stats)
crit_95 <- quantile(df_boot$stat, 0.95, na.rm = TRUE)

p_boot <- ggplot(df_boot, aes(x = stat)) +
  geom_histogram(bins = 40, fill = ok_blue, colour = "white", alpha = 0.85) +
  geom_vline(xintercept = k_card$statistic,
             colour = ok_red, linewidth = 0.9) +
  geom_vline(xintercept = crit_95,
             colour = "grey40", linetype = "dashed", linewidth = 0.5) +
  annotate("text",
           x = k_card$statistic,
           y = 80,
           label = sprintf("observed T = %.2f\np = %.2f",
                           k_card$statistic, k_card$p_value),
           hjust = -0.1, size = 3.2, family = fam, colour = ok_red) +
  labs(x = "Bootstrap statistic", y = "Frequency") +
  theme_wp()

ggsave(file.path(fig_dir, "fig-card-bootstrap.pdf"),
       p_boot, width = 5.5, height = 3.4, device = cairo_pdf)

nreps_power_kit <- 300

# -----------------------------------------------------------------------------
# Figure 2: Power curve for iv_kitagawa under increasing violation size.
# -----------------------------------------------------------------------------
set.seed(2)
n <- 500
z <- sample(0:1, n, replace = TRUE)
d <- rbinom(n, 1, 0.3 + 0.4 * z)
y <- rnorm(n, mean = d)

power_grid <- seq(0, 2.0, length.out = 7)
pc <- iv_power(y, d, z,
               method = "kitagawa",
               n_sims = nreps_power_kit,
               delta_grid = power_grid,
               n_boot = 300,
               parallel = FALSE)

p_power <- ggplot(pc, aes(delta, power)) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             colour = "grey55", linewidth = 0.35) +
  geom_line(colour = ok_blue, linewidth = 0.8) +
  geom_point(colour = ok_blue, size = 2.2) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = scales::label_number(accuracy = 0.01)) +
  labs(x = "Violation size (delta, in units of sigma)",
       y = "Rejection probability") +
  theme_wp()

ggsave(file.path(fig_dir, "fig-power-curve.pdf"),
       p_power, width = 5.5, height = 3.4, device = cairo_pdf)

# -----------------------------------------------------------------------------
# Figure 3: Monte Carlo null distribution of iv_testjfe vs chi^2_{K-2}.
# -----------------------------------------------------------------------------
set.seed(3)
K <- 20
n_sim <- 3000
n_reps <- 1000
stats_null  <- numeric(n_reps)
p_asy_null  <- numeric(n_reps)
p_boot_null <- numeric(n_reps)
for (s in seq_len(n_reps)) {
  judge <- sample.int(K, n_sim, replace = TRUE)
  p_by_j <- 0.3 + 0.02 * seq_len(K)
  d <- rbinom(n_sim, 1, p_by_j[judge])
  y <- rnorm(n_sim, mean = d)
  r_asy  <- iv_testjfe(y, d, judge, method = "asymptotic",
                       n_boot = 1, parallel = FALSE)
  stats_null[s] <- r_asy$statistic
  p_asy_null[s] <- r_asy$p_value_asymptotic
  r_boot <- iv_testjfe(y, d, judge, method = "bootstrap",
                       n_boot = 200, parallel = FALSE)
  p_boot_null[s] <- r_boot$p_value_bootstrap
}

df_null <- data.frame(stat = stats_null)
x_grid <- seq(0, max(stats_null) * 1.05, length.out = 300)
df_chi <- data.frame(x = x_grid, y = dchisq(x_grid, df = K - 2))

p_null <- ggplot(df_null, aes(x = stat)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, fill = ok_sky, colour = "white", alpha = 0.75) +
  geom_line(data = df_chi, aes(x, y),
            colour = ok_red, linewidth = 0.9) +
  labs(x = expression(paste("Observed ", T[n])),
       y = "Density") +
  theme_wp()

ggsave(file.path(fig_dir, "fig-testjfe-null.pdf"),
       p_null, width = 5.5, height = 3.4, device = cairo_pdf)

# Null-distribution summary statistics for inline replacement in the paper.
null_stats <- list(
  nreps        = n_reps,
  mean_emp     = mean(stats_null),
  var_emp      = var(stats_null),
  q95_emp      = unname(quantile(stats_null, 0.95)),
  size_asy     = mean(p_asy_null  < 0.05, na.rm = TRUE),
  size_boot    = mean(p_boot_null < 0.05, na.rm = TRUE),
  mc_se_size5  = sqrt(0.05 * 0.95 / n_reps) * 100
)
saveRDS(null_stats, file.path("paper", "null_stats.rds"))

# -----------------------------------------------------------------------------
# Figure 4: mu_j vs p_j scatter for valid and violated judge IV designs.
# -----------------------------------------------------------------------------
set.seed(4)
K <- 20
n <- 3000
judge <- sample.int(K, n, replace = TRUE)
p_by_j <- seq(0.2, 0.7, length.out = K)
d <- rbinom(n, 1, p_by_j[judge])
y_valid <- rnorm(n, mean = d)
y_viol  <- rnorm(n, mean = d + 1.5 * sin(judge * 0.5))

per_judge <- function(yv) {
  data.frame(
    j = seq_len(K),
    p = vapply(seq_len(K), function(j) mean(d[judge == j]), numeric(1)),
    mu = vapply(seq_len(K), function(j) mean(yv[judge == j]), numeric(1)),
    n_j = vapply(seq_len(K), function(j) sum(judge == j), integer(1))
  )
}
df_valid <- transform(per_judge(y_valid), design = "Valid IV")
df_viol  <- transform(per_judge(y_viol),  design = "Exclusion violation")
df_jfe <- rbind(df_valid, df_viol)
df_jfe$design <- factor(df_jfe$design,
                        levels = c("Valid IV", "Exclusion violation"))

p_jfe <- ggplot(df_jfe, aes(p, mu, size = n_j, colour = design)) +
  geom_smooth(method = "lm", formula = y ~ x,
              se = FALSE, linewidth = 0.6,
              mapping = aes(weight = n_j)) +
  geom_point(alpha = 0.75) +
  scale_colour_manual(values = c(ok_blue, ok_red), name = NULL) +
  scale_size_area(max_size = 3.5, guide = "none") +
  facet_wrap(~ design, nrow = 1) +
  labs(x = expression(paste("Per-judge propensity ", p[j])),
       y = expression(paste("Per-judge mean outcome ", mu[j]))) +
  theme_wp() +
  theme(legend.position = "none",
        strip.text = element_text(family = fam, size = 10),
        strip.background = element_blank())

ggsave(file.path(fig_dir, "fig-judges-scatter.pdf"),
       p_jfe, width = 5.5, height = 3.4, device = cairo_pdf)

# -----------------------------------------------------------------------------
# Figure 5: Pairwise LATE heatmap under violation.
# -----------------------------------------------------------------------------
r_viol <- iv_testjfe(y_viol, d, judge, n_boot = 30, parallel = FALSE)
pw <- r_viol$pairwise_late
if (!is.null(pw)) {
  pw_long <- as.data.frame(as.table(pw))
  names(pw_long) <- c("j", "k", "late")
  pw_long$late[!is.finite(pw_long$late)] <- NA
  pw_long$late_trimmed <- pmax(pmin(pw_long$late, 3), -3)

  p_pw <- ggplot(pw_long, aes(j, k, fill = late_trimmed)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = ok_red, mid = "grey97", high = ok_blue,
                         midpoint = r_viol$coef["slope"],
                         limits = c(-3, 3), na.value = "grey90",
                         name = expression(paste(hat(beta)[jk]))) +
    coord_fixed() +
    labs(x = NULL, y = NULL) +
    theme_wp_hmap()

  ggsave(file.path(fig_dir, "fig-pairwise-late.pdf"),
         p_pw, width = 5.0, height = 3.6, device = cairo_pdf)
}

# -----------------------------------------------------------------------------
# Figure 6: Power curve for iv_testjfe under sinusoidal judge-effect exclusion
# violation. Y = D + eta * sigma * sin(0.5 * J) + eps.
# -----------------------------------------------------------------------------
set.seed(6)
K <- 20
n_sim <- 3000
nreps_power_jfe <- 200
eta_grid <- c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5)
power_jfe <- numeric(length(eta_grid))
for (i in seq_along(eta_grid)) {
  rejections <- 0
  for (s in seq_len(nreps_power_jfe)) {
    judge <- sample.int(K, n_sim, replace = TRUE)
    p_by_j <- 0.3 + 0.02 * seq_len(K)
    d <- rbinom(n_sim, 1, p_by_j[judge])
    y <- rnorm(n_sim, mean = d + eta_grid[i] * sin(0.5 * judge))
    r <- iv_testjfe(y, d, judge, method = "asymptotic",
                    n_boot = 1, parallel = FALSE)
    if (is.finite(r$p_value) && r$p_value < 0.05) rejections <- rejections + 1L
  }
  power_jfe[i] <- rejections / nreps_power_jfe
}

df_power_jfe <- data.frame(eta = eta_grid, power = power_jfe)

p_power_jfe <- ggplot(df_power_jfe, aes(eta, power)) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             colour = "grey55", linewidth = 0.35) +
  geom_line(colour = ok_blue, linewidth = 0.8) +
  geom_point(colour = ok_blue, size = 2.2) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = scales::label_number(accuracy = 0.01)) +
  labs(x = expression(paste("Violation amplitude (", eta, ", units of ", sigma, ")")),
       y = "Rejection probability") +
  theme_wp()

ggsave(file.path(fig_dir, "fig-power-testjfe.pdf"),
       p_power_jfe, width = 5.5, height = 3.4, device = cairo_pdf)

# Save the MC rep counts for inline substitution in the paper.
saveRDS(
  list(nreps_power_kit = nreps_power_kit,
       nreps_power_jfe = nreps_power_jfe),
  file.path("paper", "power_stats.rds")
)

cat("\n--- figures built ---\n")
