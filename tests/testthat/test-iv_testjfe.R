test_that("iv_testjfe returns an iv_test object with the right structure", {
  set.seed(1)
  n <- 500
  judge <- sample.int(10, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.02 * judge)
  y <- rnorm(n, mean = d)
  out <- iv_testjfe(y, d, judge, n_boot = 50, parallel = FALSE)
  expect_s3_class(out, "iv_test")
  expect_identical(out$test, "Frandsen-Lefgren-Leslie (2023)")
  expect_equal(out$n_judges, 10L)
  expect_named(out$coef, c("intercept", "slope"))
})

test_that("iv_testjfe basis_order > 1 uses polynomial basis and adjusts df", {
  set.seed(1)
  K <- 15
  n_per <- 150
  judge <- rep(seq_len(K), each = n_per)
  p_j <- seq(0.2, 0.8, length.out = K)
  d <- rbinom(length(judge), 1, p_j[judge])
  y <- rnorm(length(judge), mean = d)
  r1 <- iv_testjfe(y, d, judge, n_boot = 30, parallel = FALSE,
                   basis_order = 1L)
  r2 <- iv_testjfe(y, d, judge, n_boot = 30, parallel = FALSE,
                   basis_order = 2L)
  r3 <- iv_testjfe(y, d, judge, n_boot = 30, parallel = FALSE,
                   basis_order = 3L)
  expect_equal(r1$basis_order, 1L)
  expect_equal(r2$basis_order, 2L)
  expect_equal(r3$basis_order, 3L)
  expect_true("slope" %in% names(r1$coef))
  expect_true(all(c("delta_p1", "delta_p2") %in% names(r2$coef)))
  expect_true(all(c("delta_p1", "delta_p2", "delta_p3") %in% names(r3$coef)))
})

test_that("iv_testjfe basis_order must be a positive integer", {
  set.seed(1)
  j <- sample(1:10, 200, replace = TRUE)
  d <- rbinom(200, 1, 0.5)
  y <- rnorm(200)
  expect_error(iv_testjfe(y, d, j, basis_order = "banana",
                          n_boot = 10, parallel = FALSE),
               "positive integer")
  expect_error(iv_testjfe(y, d, j, basis_order = 1.5,
                          n_boot = 10, parallel = FALSE),
               "positive integer")
})

test_that("iv_testjfe basis_order > 1 aborts for multivalued D", {
  set.seed(1)
  j <- sample(1:10, 300, replace = TRUE)
  d <- sample(0:2, 300, replace = TRUE)
  y <- rnorm(300)
  expect_error(iv_testjfe(y, d, j, basis_order = 2L,
                          n_boot = 10, parallel = FALSE),
               "binary")
})

test_that("iv_testjfe requires at least three judges", {
  set.seed(1)
  n <- 100
  judge <- sample(1:2, n, replace = TRUE)
  d <- rbinom(n, 1, 0.4)
  y <- rnorm(n)
  expect_error(iv_testjfe(y, d, judge, n_boot = 10, parallel = FALSE),
               "at least three")
})

test_that("iv_testjfe null distribution is approximately chi^2_{K-2}", {
  skip_on_cran()
  K <- 20
  stats <- numeric(60)
  p_by_j <- 0.3 + 0.02 * seq_len(K)
  for (s in seq_along(stats)) {
    set.seed(1000 + s)
    n <- 2000
    judge <- sample.int(K, n, replace = TRUE)
    d <- rbinom(n, 1, p_by_j[judge])
    y <- rnorm(n, mean = d)
    r <- iv_testjfe(y, d, judge, n_boot = 10, parallel = FALSE)
    stats[s] <- r$statistic
  }
  # Under chi^2_{K-2 = 18}, mean ~ 18 and sd ~ 6. Allow generous slack.
  expect_gt(mean(stats), 14)
  expect_lt(mean(stats), 22)
})

test_that("iv_testjfe rejects a large exclusion violation", {
  skip_on_cran()
  set.seed(2)
  K <- 20
  n <- 3000
  judge <- sample.int(K, n, replace = TRUE)
  p_by_j <- 0.3 + 0.02 * seq_len(K)
  d <- rbinom(n, 1, p_by_j[judge])
  # Strong non-linear judge effect on Y
  y <- rnorm(n, mean = d + 1.5 * sin(judge * 0.5))
  out <- iv_testjfe(y, d, judge, n_boot = 10, parallel = FALSE)
  expect_lt(out$p_value, 0.05)
})

test_that("iv_testjfe handles multivalued treatment", {
  set.seed(1)
  K <- 12
  n <- 2000
  judge <- sample.int(K, n, replace = TRUE)
  # Multivalued D: 0, 1, 2
  d <- vapply(judge, function(j) {
    u <- runif(1)
    p_hi <- 0.1 + 0.03 * j
    p_md <- 0.4
    if (u < p_hi) 2L else if (u < p_hi + p_md) 1L else 0L
  }, integer(1))
  y <- rnorm(n, mean = d)
  out <- iv_testjfe(y, d, judge, n_boot = 30, parallel = FALSE)
  expect_s3_class(out, "iv_test")
  expect_equal(out$n_treatment_levels, 3L)
  expect_null(out$pairwise_late)
  expect_null(out$worst_pair)
  expect_named(out$coef, c("intercept", "beta_d1", "beta_d2"))
  expect_true(is.finite(out$p_value))
  expect_gt(out$p_value, 0.05)
})

test_that("iv_testjfe multivalued rejects a large violation", {
  skip_on_cran()
  set.seed(2)
  K <- 15
  n <- 3000
  judge <- sample.int(K, n, replace = TRUE)
  d <- vapply(judge, function(j) {
    u <- runif(1)
    p_hi <- 0.1 + 0.03 * j
    p_md <- 0.4
    if (u < p_hi) 2L else if (u < p_hi + p_md) 1L else 0L
  }, integer(1))
  # Direct non-linear judge effect
  y <- rnorm(n, mean = d + 1.0 * sin(judge))
  out <- iv_testjfe(y, d, judge, n_boot = 30, parallel = FALSE)
  expect_lt(out$p_value, 0.05)
})

test_that("iv_testjfe accepts covariates and residualises internally", {
  set.seed(1)
  K <- 10
  n <- 1000
  judge <- sample.int(K, n, replace = TRUE)
  p_by_j <- 0.3 + 0.04 * seq_len(K)
  d <- rbinom(n, 1, p_by_j[judge])
  x <- rnorm(n)
  y <- rnorm(n, mean = d + 0.5 * x)
  out <- iv_testjfe(y, d, judge, x = x, n_boot = 50, parallel = FALSE)
  expect_s3_class(out, "iv_test")
  expect_equal(out$n_judges, K)
})
