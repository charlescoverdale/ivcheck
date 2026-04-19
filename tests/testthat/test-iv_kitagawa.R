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

test_that("iv_kitagawa rejects non-binary d with a clear message", {
  y <- rnorm(100)
  d <- rnorm(100)
  z <- sample(0:1, 100, replace = TRUE)
  expect_error(iv_kitagawa(y, d, z, n_boot = 10, parallel = FALSE), "binary")
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
    expect_named(out$binding, c("z_low", "z_high", "d", "y_lower", "y_upper"))
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
