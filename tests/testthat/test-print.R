test_that("print.iv_test runs without error on stub object", {
  set.seed(1)
  n <- 100
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.5)
  y <- rnorm(n)
  out <- iv_kitagawa(y, d, z, n_boot = 10)
  expect_no_error(print(out))
  expect_invisible(print(out))
})

test_that("format.iv_test returns a character string", {
  set.seed(1)
  n <- 100
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.5)
  y <- rnorm(n)
  out <- iv_kitagawa(y, d, z, n_boot = 10)
  expect_type(format(out), "character")
  expect_length(format(out), 1L)
})
