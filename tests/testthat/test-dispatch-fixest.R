test_that("extract_iv_data.fixest returns y, d, z from a fitted feols IV model", {
  skip_if_not_installed("fixest")
  set.seed(1)
  n <- 400
  df <- data.frame(
    z = sample(0:1, n, replace = TRUE),
    x = rnorm(n)
  )
  df$d <- rbinom(n, 1, 0.3 + 0.4 * df$z)
  df$y <- rnorm(n, mean = df$d + 0.2 * df$x)
  m <- fixest::feols(y ~ x | d ~ z, data = df)
  yz <- ivcheck:::extract_iv_data(m)
  expect_equal(length(yz$y), n)
  expect_equal(length(yz$d), n)
  expect_equal(length(yz$z), n)
})

test_that("iv_kitagawa dispatches on a fitted fixest IV model", {
  skip_if_not_installed("fixest")
  set.seed(1)
  n <- 400
  df <- data.frame(z = sample(0:1, n, replace = TRUE))
  df$d <- rbinom(n, 1, 0.3 + 0.4 * df$z)
  df$y <- rnorm(n, mean = df$d)
  m <- fixest::feols(y ~ 1 | d ~ z, data = df)
  out <- iv_kitagawa(m, n_boot = 50, parallel = FALSE)
  expect_s3_class(out, "iv_test")
  expect_equal(out$n, n)
})

test_that("iv_check dispatches on a fitted fixest IV model and returns a tidy table", {
  skip_if_not_installed("fixest")
  set.seed(1)
  n <- 400
  df <- data.frame(z = sample(0:1, n, replace = TRUE))
  df$d <- rbinom(n, 1, 0.3 + 0.4 * df$z)
  df$y <- rnorm(n, mean = df$d)
  m <- fixest::feols(y ~ 1 | d ~ z, data = df)
  chk <- iv_check(m, n_boot = 30, parallel = FALSE)
  expect_s3_class(chk, "iv_check")
  expect_true(nrow(chk$table) >= 1L)
  expect_named(chk$table, c("test", "statistic", "p_value", "verdict"))
})

test_that("extract_iv_data.fixest errors on non-IV fixest models", {
  skip_if_not_installed("fixest")
  set.seed(1)
  n <- 200
  df <- data.frame(y = rnorm(n), x = rnorm(n))
  m <- fixest::feols(y ~ x, data = df)
  expect_error(ivcheck:::extract_iv_data(m), "not an IV model")
})
