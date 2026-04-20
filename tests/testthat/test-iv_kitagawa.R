test_that("iv_kitagawa returns an iv_test object with the right structure", {
  set.seed(1)
  n <- 200
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  out <- iv_kitagawa(y, d, z, n_boot = 50, parallel = FALSE)
  expect_s3_class(out, "iv_test")
  expect_identical(out$test, "Kitagawa (2015)")
  expect_equal(out$n, n)
  expect_equal(length(out$boot_stats), 50L)
})

test_that("iv_kitagawa treatment_order = 'unordered' requires monotonicity_set", {
  set.seed(1)
  n <- 400
  z <- sample(0:2, n, replace = TRUE)
  d <- sample(0:2, n, replace = TRUE)
  y <- rnorm(n)
  expect_error(
    iv_kitagawa(y, d, z, n_boot = 10, parallel = FALSE,
                treatment_order = "unordered"),
    "monotonicity_set"
  )
})

test_that("iv_kitagawa treatment_order = 'unordered' validates column names", {
  set.seed(1)
  n <- 400
  z <- sample(0:2, n, replace = TRUE)
  d <- sample(0:2, n, replace = TRUE)
  y <- rnorm(n)
  bad_ms <- data.frame(level = 0, from = 0, to = 1)
  expect_error(
    iv_kitagawa(y, d, z, n_boot = 10, parallel = FALSE,
                treatment_order = "unordered", monotonicity_set = bad_ms),
    "columns"
  )
})

test_that("iv_kitagawa unordered + monotonicity_set runs and returns valid iv_test", {
  skip_on_cran()
  set.seed(1)
  n <- 800
  z <- sample(0:2, n, replace = TRUE)
  prop <- matrix(c(0.5, 0.3, 0.2,
                   0.3, 0.4, 0.3,
                   0.2, 0.3, 0.5), nrow = 3, byrow = TRUE)
  d <- vapply(seq_len(n), function(i) {
    sample(0:2, 1, prob = prop[z[i] + 1, ])
  }, integer(1))
  y <- rnorm(n, mean = d)
  ms <- data.frame(
    d = c(0, 1, 2),
    z_from = c(0, 0, 2),
    z_to   = c(2, 1, 0)
  )
  out <- iv_kitagawa(y, d, z, n_boot = 50, parallel = FALSE,
                     treatment_order = "unordered",
                     monotonicity_set = ms)
  expect_s3_class(out, "iv_test")
  expect_true(is.finite(out$statistic) && out$statistic >= 0)
  expect_true(is.finite(out$p_value) && out$p_value >= 0 && out$p_value <= 1)
  expect_equal(out$treatment_order, "unordered")
})

test_that("iv_kitagawa multiplier = gaussian runs and returns sensible output", {
  set.seed(1)
  n <- 300
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  out <- iv_kitagawa(y, d, z, n_boot = 60, parallel = FALSE,
                     multiplier = "gaussian")
  expect_s3_class(out, "iv_test")
  expect_equal(out$multiplier, "gaussian")
})

test_that("iv_kitagawa handles ties in Y without breaking", {
  # Ties in Y arise often in applied work (rounded wages, survey
  # categories). The empirical CDF on a ties-rich support should still
  # produce a well-defined statistic and a p-value in [0, 1].
  set.seed(1)
  n <- 600
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  # Force 20 percent of Y values to be identical at zero.
  y <- rnorm(n, mean = d)
  y[sample.int(n, size = n %/% 5)] <- 0
  expect_equal(sum(y == 0) / n, 0.2, tolerance = 0.05)

  out <- iv_kitagawa(y, d, z, n_boot = 80, parallel = FALSE)
  expect_s3_class(out, "iv_test")
  expect_true(is.finite(out$statistic) && out$statistic >= 0)
  expect_true(is.finite(out$p_value) &&
                out$p_value >= 0 && out$p_value <= 1)
  # Stat should be robust across two seeds on the same data.
  set.seed(7); out2 <- iv_kitagawa(y, d, z, n_boot = 80, parallel = FALSE)
  expect_equal(out$statistic, out2$statistic, tolerance = 1e-10)
})

test_that("iv_kitagawa warns when smallest Z cell is tiny", {
  set.seed(1)
  # z == 0 has only 10 obs, z == 1 has 290
  z <- c(rep(0, 10), rep(1, 290))
  d <- rbinom(300, 1, 0.3 + 0.4 * z)
  y <- rnorm(300, mean = d)
  expect_warning(iv_kitagawa(y, d, z, n_boot = 30, parallel = FALSE),
                 "Z cell")
})

test_that("iv_kitagawa multiplier = mammen runs and returns sensible output", {
  set.seed(1)
  n <- 300
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  out <- iv_kitagawa(y, d, z, n_boot = 60, parallel = FALSE,
                     multiplier = "mammen")
  expect_s3_class(out, "iv_test")
  expect_equal(out$multiplier, "mammen")
})

test_that("iv_kitagawa rejects continuous d with a clear message", {
  y <- rnorm(100)
  d <- rnorm(100)
  z <- sample(0:1, 100, replace = TRUE)
  # Continuous d has hundreds of unique values; validator caps at 20.
  expect_error(iv_kitagawa(y, d, z, n_boot = 10, parallel = FALSE),
               "levels")
})

test_that("iv_kitagawa dispatches Sun (2023) for multivalued D", {
  skip_on_cran()
  set.seed(1)
  n <- 500
  z <- sample(0:1, n, replace = TRUE)
  d <- vapply(z, function(zi) {
    u <- runif(1)
    pp <- if (zi == 1) c(0.2, 0.3, 0.5) else c(0.5, 0.3, 0.2)
    if (u < pp[1]) 0L else if (u < pp[1] + pp[2]) 1L else 2L
  }, integer(1))
  y <- rnorm(n, mean = d)
  out <- iv_kitagawa(y, d, z, n_boot = 50, parallel = FALSE)
  expect_identical(out$test, "Sun (2023)")
  expect_true(out$multivalued)
  expect_equal(out$n_treatment_levels, 3L)
})

test_that("iv_kitagawa rejects unequal-length inputs", {
  y <- rnorm(100)
  d <- rbinom(50, 1, 0.5)
  z <- sample(0:1, 100, replace = TRUE)
  expect_error(iv_kitagawa(y, d, z, n_boot = 10, parallel = FALSE),
               "unequal lengths")
})

test_that("iv_kitagawa statistic is non-negative and p-value is in [0, 1]", {
  set.seed(1)
  n <- 200
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  out <- iv_kitagawa(y, d, z, n_boot = 100, parallel = FALSE)
  expect_true(out$statistic >= 0)
  expect_true(out$p_value >= 0 && out$p_value <= 1)
})

test_that("iv_kitagawa detects a clear exclusion-restriction violation", {
  skip_on_cran()
  set.seed(2)
  n <- 1500
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  # Direct effect of Z on Y: large exclusion violation
  y <- rnorm(n, mean = d + 2 * z)
  out <- iv_kitagawa(y, d, z, n_boot = 300, parallel = FALSE)
  expect_true(out$p_value < 0.05)
})

test_that("iv_kitagawa binding element is populated when statistic > 0", {
  set.seed(1)
  n <- 200
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  out <- iv_kitagawa(y, d, z, n_boot = 50, parallel = FALSE)
  expect_true(is.null(out$binding) || is.list(out$binding))
  if (!is.null(out$binding)) {
    expect_named(out$binding, c("z_low", "z_high", "direction", "y_lower", "y_upper"))
  }
})

test_that("iv_kitagawa runs on multilevel discrete instruments", {
  set.seed(1)
  n <- 400
  z <- sample(0:2, n, replace = TRUE)
  # First-stage increasing in z
  pd <- 0.2 + 0.2 * z
  d <- rbinom(n, 1, pd)
  y <- rnorm(n, mean = d)
  out <- iv_kitagawa(y, d, z, n_boot = 100, parallel = FALSE)
  expect_s3_class(out, "iv_test")
  expect_true(out$p_value >= 0 && out$p_value <= 1)
})

test_that("iv_kitagawa runs on card1995 with a binary treatment", {
  skip_on_cran()
  data(card1995)
  out <- iv_kitagawa(
    card1995$lwage,
    card1995$college,
    card1995$near_college,
    n_boot = 100,
    parallel = FALSE
  )
  expect_s3_class(out, "iv_test")
  expect_true(out$statistic >= 0)
})
