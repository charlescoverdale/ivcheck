test_that("extract_iv_data.ivreg returns y, d, z from a fitted ivreg model", {
  skip_if_not_installed("ivreg")
  set.seed(1)
  n <- 400
  df <- data.frame(
    z = sample(0:1, n, replace = TRUE),
    x = rnorm(n)
  )
  df$d <- rbinom(n, 1, 0.3 + 0.4 * df$z)
  df$y <- rnorm(n, mean = df$d + 0.2 * df$x)
  m <- ivreg::ivreg(y ~ x + d | x + z, data = df)
  yz <- ivcheck:::extract_iv_data(m)
  expect_equal(length(yz$y), n)
  expect_equal(length(yz$d), n)
  expect_equal(length(yz$z), n)
})

test_that("iv_kitagawa dispatches on a fitted ivreg model", {
  skip_if_not_installed("ivreg")
  set.seed(1)
  n <- 400
  df <- data.frame(z = sample(0:1, n, replace = TRUE))
  df$d <- rbinom(n, 1, 0.3 + 0.4 * df$z)
  df$y <- rnorm(n, mean = df$d)
  m <- ivreg::ivreg(y ~ d | z, data = df)
  out <- iv_kitagawa(m, n_boot = 50, parallel = FALSE)
  expect_s3_class(out, "iv_test")
  expect_equal(out$n, n)
})
