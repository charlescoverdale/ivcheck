test_that("iv_power returns a data frame with delta and power columns", {
  set.seed(1)
  n <- 200
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  out <- iv_power(y, d, z, method = "kitagawa", n_sims = 5,
                  delta_grid = c(0, 1), n_boot = 30, parallel = FALSE)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("delta", "power"))
  expect_equal(nrow(out), 2L)
  expect_true(all(out$power >= 0 & out$power <= 1))
})

test_that("iv_power method validation accepts only registered tests", {
  set.seed(1)
  n <- 100
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.5)
  y <- rnorm(n)
  expect_error(iv_power(y, d, z, method = "banana",
                        n_sims = 5, delta_grid = c(0),
                        n_boot = 20, parallel = FALSE))
})

test_that("iv_power approaches 1 for large deltas (slow test)", {
  skip_on_cran()
  set.seed(1)
  n <- 500
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  out <- iv_power(y, d, z, method = "kitagawa", n_sims = 30,
                  delta_grid = c(0, 2), n_boot = 100, parallel = FALSE)
  expect_lt(out$power[out$delta == 0], 0.2)
  expect_gt(out$power[out$delta == 2], 0.8)
})
