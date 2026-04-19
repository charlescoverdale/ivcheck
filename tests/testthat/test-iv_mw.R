test_that("iv_mw returns an iv_test object with the right structure", {
  set.seed(1)
  n <- 200
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  out <- iv_mw(y, d, z, n_boot = 50)
  expect_s3_class(out, "iv_test")
  expect_identical(out$test, "Mourifie-Wan (2017)")
  expect_false(out$conditional)
})

test_that("iv_mw sets conditional = TRUE when covariates supplied", {
  set.seed(1)
  n <- 200
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  x <- rnorm(n)
  out <- iv_mw(y, d, z, x = x, n_boot = 50)
  expect_true(out$conditional)
})
