test_that("iv_mw returns an iv_test object with the right structure", {
  set.seed(1)
  n <- 200
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  out <- iv_mw(y, d, z, n_boot = 50, parallel = FALSE)
  expect_s3_class(out, "iv_test")
  expect_identical(out$test, "Mourifie-Wan (2017)")
  expect_false(out$conditional)
  expect_true(is.na(out$kappa_n))
})

test_that("iv_mw with x flags conditional = TRUE and reports kappa_n", {
  set.seed(1)
  n <- 300
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  x <- rnorm(n)
  out <- iv_mw(y, d, z, x = x, n_boot = 50, parallel = FALSE)
  expect_true(out$conditional)
  expect_true(is.numeric(out$kappa_n) && is.finite(out$kappa_n))
  expect_equal(out$basis_order, 3L)
})

test_that("iv_mw without x matches iv_kitagawa closely on the same data", {
  set.seed(1)
  n <- 300
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  r_k  <- iv_kitagawa(y, d, z, n_boot = 100, parallel = FALSE)
  r_mw <- iv_mw(y, d, z, n_boot = 100, parallel = FALSE)
  expect_equal(r_k$statistic, r_mw$statistic)
})

test_that("iv_mw errors when x rows do not match y", {
  set.seed(1)
  n <- 100
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.5)
  y <- rnorm(n)
  x <- rnorm(50)
  expect_error(iv_mw(y, d, z, x = x, n_boot = 10, parallel = FALSE),
               "same number of rows")
})

test_that("iv_mw statistic and p-value are in valid ranges", {
  set.seed(1)
  n <- 200
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  out <- iv_mw(y, d, z, n_boot = 50, parallel = FALSE)
  expect_true(out$statistic >= 0)
  expect_true(out$p_value >= 0 && out$p_value <= 1)
})
