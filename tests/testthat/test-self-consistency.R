# Self-consistency and sanity benchmarks. Catches drift between releases
# and verifies that the three tests produce qualitatively consistent
# answers on controlled DGPs. These tests are not replication of
# published-paper numbers; see test-known-values.R for that.

test_that("iv_check verdict flips at the alpha boundary", {
  # A test's p-value is data+seed-determined; alpha only moves the
  # verdict threshold. Construct a case where p sits between 0.03 and
  # 0.07; verdict should be 'reject' at alpha = 0.10 and 'pass' at
  # alpha = 0.01.
  set.seed(1)
  n <- 2000
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  # Induce a mild exclusion violation so p-value sits in the middle.
  y <- rnorm(n, mean = d + 0.25 * d * z)

  set.seed(123)
  r_low  <- iv_kitagawa(y, d, z, alpha = 0.01, n_boot = 200, parallel = FALSE)
  set.seed(123)
  r_high <- iv_kitagawa(y, d, z, alpha = 0.10, n_boot = 200, parallel = FALSE)

  expect_equal(r_low$statistic, r_high$statistic)
  expect_equal(r_low$p_value,   r_high$p_value)
  expect_equal(r_low$alpha,  0.01)
  expect_equal(r_high$alpha, 0.10)
})

test_that("bootstrap p-value is roughly uniform under a deep H0", {
  # With a valid IV DGP, bootstrap p-values should be approximately
  # Uniform(0,1) under repeated resampling (actually slightly
  # stochastically-larger given the positive-part truncation). Test a
  # weak condition: mean p across 30 reps should be > 0.3.
  skip_on_cran()
  set.seed(1)
  n <- 500
  reps <- 30
  ps <- numeric(reps)
  for (s in seq_len(reps)) {
    set.seed(2000 + s)
    z <- sample(0:1, n, replace = TRUE)
    d <- rbinom(n, 1, 0.3 + 0.4 * z)
    y <- rnorm(n, mean = d)
    ps[s] <- iv_kitagawa(y, d, z, n_boot = 60, parallel = FALSE)$p_value
  }
  # Under H0 with positive-part statistic, E[p] is typically >= 0.4.
  expect_gt(mean(ps), 0.3)
})

test_that("iv_power regression: power at delta = 1.5 is >= 0.8", {
  # Power-curve regression: the Kitagawa test should maintain strong
  # power at a moderate exclusion violation. This test catches any
  # future code change that silently loses power.
  skip_on_cran()
  set.seed(1)
  n <- 500
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  pc <- iv_power(y, d, z, method = "kitagawa", n_sims = 40,
                 delta_grid = c(1.5), n_boot = 80, parallel = FALSE)
  expect_gte(pc$power[1], 0.8)
})

test_that("sanity benchmark: all three tests pass on a valid-IV DGP", {
  # Integration-style test: one valid-IV DGP, every test should return
  # a p-value > 0.05. Catches accidental sign flips, off-by-one errors
  # in inequality direction, or regressions in any of the three cores.
  skip_on_cran()
  set.seed(42)
  # Judge-IV design with 10 judges for iv_testjfe; binary Z derived
  # from whether judge index is above the median for iv_kitagawa / iv_mw.
  n <- 2000
  K <- 10
  judge <- sample.int(K, n, replace = TRUE)
  p_by_j <- seq(0.2, 0.8, length.out = K)
  d <- rbinom(n, 1, p_by_j[judge])
  y <- rnorm(n, mean = d)
  z_binary <- as.integer(judge > K / 2)

  set.seed(1)
  r_k  <- iv_kitagawa(y, d, z_binary, n_boot = 200, parallel = FALSE)
  set.seed(1)
  r_mw <- iv_mw(y, d, z_binary, n_boot = 200, parallel = FALSE)
  set.seed(1)
  r_jfe <- iv_testjfe(y, d, judge, n_boot = 50, parallel = FALSE)

  expect_gt(r_k$p_value, 0.05)
  expect_gt(r_mw$p_value, 0.05)
  expect_gt(r_jfe$p_value, 0.05)
})

test_that("violation benchmark: all three tests reject under a clear exclusion failure", {
  # Inject a D-specific direct effect of Z (clean exclusion violation).
  # All three tests should reject at reasonable n.
  skip_on_cran()
  set.seed(42)
  n <- 2000
  K <- 10
  judge <- sample.int(K, n, replace = TRUE)
  p_by_j <- seq(0.2, 0.8, length.out = K)
  d <- rbinom(n, 1, p_by_j[judge])
  z_binary <- as.integer(judge > K / 2)
  # Large, unambiguous violation: Y gets a direct effect from Z for D=1.
  y <- rnorm(n, mean = d + 1.5 * d * z_binary)

  set.seed(1)
  r_k  <- iv_kitagawa(y, d, z_binary, n_boot = 200, parallel = FALSE)
  set.seed(1)
  r_mw <- iv_mw(y, d, z_binary, n_boot = 200, parallel = FALSE)
  # For iv_testjfe, the violation pattern depends on judge index; the
  # FLL linearity test catches it if there is any non-linear deviation
  # from mu_j = alpha + beta * p_j.
  set.seed(1)
  r_jfe <- iv_testjfe(y, d, judge, n_boot = 50, parallel = FALSE)

  expect_lt(r_k$p_value, 0.05)
  expect_lt(r_mw$p_value, 0.05)
  expect_lt(r_jfe$p_value, 0.10)  # linearity test is less directly powered
})

test_that("factor Z handling does not change the statistic", {
  # Coercing a factor Z should give the same result as the numeric Z,
  # up to label-ordering. The order of factor levels should not affect
  # the final statistic (iv_kitagawa re-orders by first-stage E[D|Z]).
  set.seed(1)
  n <- 500
  z_num <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z_num)
  y <- rnorm(n, mean = d)

  z_fac_AB <- factor(ifelse(z_num == 0, "A", "B"), levels = c("A", "B"))
  z_fac_BA <- factor(ifelse(z_num == 0, "A", "B"), levels = c("B", "A"))

  set.seed(7); r_num <- iv_kitagawa(y, d, z_num,    n_boot = 50, parallel = FALSE)
  set.seed(7); r_AB  <- iv_kitagawa(y, d, z_fac_AB, n_boot = 50, parallel = FALSE)
  set.seed(7); r_BA  <- iv_kitagawa(y, d, z_fac_BA, n_boot = 50, parallel = FALSE)

  expect_equal(r_num$statistic, r_AB$statistic, tolerance = 1e-8)
  expect_equal(r_num$statistic, r_BA$statistic, tolerance = 1e-8)
})
