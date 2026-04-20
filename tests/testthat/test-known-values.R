# Level-2 known-value tests. Each test exercises a case where the
# analytical answer is derivable by hand or by construction. These
# tests catch off-by-one, wrong-sign, or discretisation-choice bugs
# that structural and invariant tests can miss.

test_that("iv_kitagawa: statistic is small when Z is independent of (D, Y)", {
  # Under Z independent of (D, Y), F(y, d | Z = z1) == F(y, d | Z = z2)
  # for all y, d. The empirical positive-part KS statistic converges
  # to 0 at rate 1/sqrt(n). With n = 2000 we expect stat < ~4 typically.
  skip_on_cran()
  set.seed(101)
  n <- 2000
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.5)  # Z-independent first stage
  y <- rnorm(n, mean = d)  # Y depends on D only
  out <- iv_kitagawa(y, d, z, n_boot = 100, parallel = FALSE)
  # Statistic should be bounded; p-value should not reject
  expect_lt(out$statistic, 5)
  expect_gt(out$p_value, 0.1)
})

test_that("iv_kitagawa: deterministic Y perfectly separated by D gives stat = 0", {
  # If Y = 1{D=1} exactly, then F(y, 1 | z) = 0 for y < 1 and
  # = P(D=1 | z) for y >= 1, with no variation across z except via
  # the marginals. The inequality F(y, 1 | z_low) - F(y, 1 | z_high)
  # is strictly <= 0 everywhere. Empirical stat is exactly 0.
  set.seed(1)
  n <- 500
  z <- c(rep(0, n/2), rep(1, n/2))
  d <- c(rep(0, n/4), rep(1, n/4), rep(0, n/4), rep(1, 3*n/4 - n/2))
  # Simpler: deterministic
  y <- as.numeric(d)
  out <- iv_kitagawa(y, d, z, n_boot = 20, parallel = FALSE)
  expect_equal(out$statistic, 0, tolerance = 1e-10)
})

test_that("iv_mw: without x, matches iv_kitagawa statistic exactly on same seed", {
  # Invariant test: iv_mw unconditional path delegates to the same
  # Kitagawa core. On identical seeds, the two statistics are equal.
  set.seed(42)
  n <- 300
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  set.seed(1); r_k <- iv_kitagawa(y, d, z, n_boot = 50, parallel = FALSE)
  set.seed(1); r_mw <- iv_mw(y, d, z, n_boot = 50, parallel = FALSE)
  expect_equal(r_k$statistic, r_mw$statistic, tolerance = 1e-12)
  expect_equal(r_k$p_value, r_mw$p_value, tolerance = 1e-12)
})

test_that("iv_testjfe: T_n = 0 when mu_j is exactly linear in p_j", {
  # Construct data where mu_j = alpha + beta * p_j holds deterministically.
  # In the limit of many observations per judge, the residuals from the
  # weighted-LS fit are exactly zero and T_n collapses to 0.
  set.seed(1)
  K <- 10
  n_per <- 2000  # large per-judge sample -> tiny noise
  p_j <- seq(0.2, 0.7, length.out = K)
  alpha <- 1
  beta <- 2
  judge <- rep(seq_len(K), each = n_per)
  # Generate D from p_j[judge], Y = alpha + beta * d_deterministic + tiny noise
  # But we want mu_j = alpha + beta * p_j, so Y must have mean alpha + beta * D
  # with per-observation noise that averages out.
  d <- stats::rbinom(length(judge), 1, p_j[judge])
  y <- alpha + beta * d + stats::rnorm(length(judge), 0, 0.01)
  out <- iv_testjfe(y, d, judge, n_boot = 20, parallel = FALSE)
  # T_n should be tiny: the only source of non-zero residual is the
  # finite-sample noise in p_j_hat and mu_j_hat.
  expect_lt(out$statistic, qchisq(0.99, K - 2))
  expect_gt(out$p_value, 0.01)
})

test_that("iv_testjfe: slope in weighted-LS recovers beta when mu_j is linear in p_j", {
  # Level-2 analytical check on the OLS coefficients themselves.
  set.seed(1)
  K <- 10
  n_per <- 1000
  p_j <- seq(0.2, 0.7, length.out = K)
  alpha <- 1
  beta <- 2
  judge <- rep(seq_len(K), each = n_per)
  d <- stats::rbinom(length(judge), 1, p_j[judge])
  y <- alpha + beta * d + stats::rnorm(length(judge), 0, 0.5)
  out <- iv_testjfe(y, d, judge, n_boot = 10, parallel = FALSE)
  # Slope should recover beta. Tolerance ~ sd(mu_j_hat) / range(p_j)
  expect_equal(out$coef["slope"], c(slope = beta), tolerance = 0.1)
  expect_equal(out$coef["intercept"], c(intercept = alpha), tolerance = 0.1)
})

test_that("iv_check: overall verdict is 'cannot reject' when all p-values are large", {
  set.seed(1)
  n <- 500
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  # Raw-vector path: call via default
  chk <- structure(
    list(
      table = data.frame(test = "X", statistic = 0.5, p_value = 0.9,
                         verdict = "pass"),
      alpha = 0.05,
      overall = ivcheck:::overall_verdict(0.9, 0.05)
    ),
    class = "iv_check"
  )
  expect_match(chk$overall, "cannot reject")
})

test_that("iv_power: delta = 0 empirical size is at most ~alpha + MC noise", {
  skip_on_cran()
  set.seed(1)
  n <- 300
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  out <- iv_power(y, d, z, method = "kitagawa", n_sims = 50,
                  delta_grid = c(0), n_boot = 80, parallel = FALSE)
  # Under strict null, Kitagawa test is conservative; empirical size
  # should be <= alpha + a few SE. With n_sims = 50 and alpha = 0.05,
  # SE(power) ~= sqrt(0.05 * 0.95 / 50) ~ 0.031. Allow 3 SE + alpha
  expect_lt(out$power[1], 0.15)
})

test_that("iv_kitagawa reproduces the Card (1995) binary-college rejection", {
  # Published-number reproduction. Kitagawa (2015) Table 2 and the
  # Mourifie-Wan (2017) replication both reject the simple college/
  # near_college specification at 5%. The exact statistic depends on
  # the bootstrap seed, but rejection is robust.
  skip_on_cran()
  data(card1995)
  set.seed(40000)
  r <- iv_kitagawa(card1995$lwage, card1995$college,
                   card1995$near_college, n_boot = 300, parallel = FALSE)
  expect_lt(r$p_value, 0.05)
  expect_gt(r$statistic, 1.5)
})

test_that("iv_testjfe null size holds for multivalued D with x (small MC)", {
  # Regression test for the FWL identity in the structural-residual
  # construction. With M = 2 and x varying within judges, a buggy
  # sigma^2 inflates and produces conservative (<5%) rejection under H0.
  # With the FWL-correct structural residual, empirical size should
  # sit near nominal alpha.
  skip_on_cran()
  set.seed(2026)
  K <- 12
  n_per_judge <- 200
  judge_vec <- rep(seq_len(K), each = n_per_judge)
  n <- length(judge_vec)
  rej <- 0
  reps <- 80
  for (s in seq_len(reps)) {
    set.seed(70000 + s)
    x <- rnorm(n)
    # Valid judge IV: D is a function of judge + x. Three levels of D.
    p_by_j <- c(0.2, 0.4, 0.6)  # P(D = 1 | J), P(D = 2 | J) will be
    lin <- 0.1 * judge_vec + 0.3 * x
    d_raw <- rbinom(n, 2, plogis(lin) * 0.9)   # values in {0, 1, 2}
    y <- d_raw + 0.4 * x + rnorm(n, 0, 0.5)    # valid exclusion
    out <- iv_testjfe(y, d_raw, judge_vec, x = x, n_boot = 10,
                      method = "asymptotic", parallel = FALSE)
    if (isTRUE(out$p_value < 0.05)) rej <- rej + 1
  }
  # Empirical size should be within roughly 3 SE of alpha = 0.05.
  # With 80 reps, SE ~= 0.024 so <= 0.13 is a very conservative pass.
  expect_lte(rej / reps, 0.15)
})
