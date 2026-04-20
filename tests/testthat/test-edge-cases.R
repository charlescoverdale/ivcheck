# Edge-case tests. Hostile inputs a PhD student might try:
# empty vectors, length-1, NA, constant treatment, constant instrument,
# negative weights, misaligned lengths, factors with unused levels.

test_that("iv_kitagawa rejects empty input", {
  expect_error(iv_kitagawa(numeric(0), integer(0), integer(0),
                           n_boot = 10, parallel = FALSE))
})

test_that("iv_kitagawa rejects length-1 input", {
  expect_error(iv_kitagawa(1, 1L, 1L, n_boot = 10, parallel = FALSE))
})

test_that("iv_kitagawa rejects length mismatches", {
  expect_error(
    iv_kitagawa(rnorm(100), rbinom(50, 1, 0.5), rbinom(100, 1, 0.5),
                n_boot = 10, parallel = FALSE),
    "unequal lengths"
  )
})

test_that("iv_kitagawa rejects NA in y, d, or z", {
  set.seed(1)
  y <- rnorm(100); d <- rbinom(100, 1, 0.5); z <- rbinom(100, 1, 0.5)
  y_na <- y; y_na[10] <- NA
  d_na <- d; d_na[5] <- NA
  z_na <- z; z_na[2] <- NA
  expect_error(iv_kitagawa(y_na, d, z, n_boot = 10, parallel = FALSE), "NA")
  expect_error(iv_kitagawa(y, d_na, z, n_boot = 10, parallel = FALSE), "NA")
  expect_error(iv_kitagawa(y, d, z_na, n_boot = 10, parallel = FALSE), "NA")
})

test_that("iv_kitagawa rejects constant Z", {
  y <- rnorm(100)
  d <- rbinom(100, 1, 0.5)
  z <- rep(1, 100)
  expect_error(iv_kitagawa(y, d, z, n_boot = 10, parallel = FALSE),
               "distinct")
})

test_that("iv_kitagawa rejects non-binary d with continuous values", {
  y <- rnorm(100)
  d <- rnorm(100)  # continuous, not discrete
  z <- rbinom(100, 1, 0.5)
  expect_error(iv_kitagawa(y, d, z, n_boot = 10, parallel = FALSE),
               "levels")
})

test_that("iv_kitagawa handles factor Z by coercing", {
  set.seed(1)
  y <- rnorm(200)
  d <- rbinom(200, 1, 0.5)
  z_fac <- factor(sample(c("A", "B"), 200, replace = TRUE))
  out <- iv_kitagawa(y, d, z_fac, n_boot = 30, parallel = FALSE)
  expect_s3_class(out, "iv_test")
})

test_that("iv_kitagawa rejects negative weights", {
  set.seed(1)
  y <- rnorm(100); d <- rbinom(100, 1, 0.5); z <- rbinom(100, 1, 0.5)
  w <- c(-1, rep(1, 99))
  expect_error(
    iv_kitagawa(y, d, z, n_boot = 10, parallel = FALSE, weights = w),
    "non-negative"
  )
})

test_that("iv_kitagawa rejects weights of wrong length", {
  set.seed(1)
  y <- rnorm(100); d <- rbinom(100, 1, 0.5); z <- rbinom(100, 1, 0.5)
  w <- rep(1, 50)
  expect_error(
    iv_kitagawa(y, d, z, n_boot = 10, parallel = FALSE, weights = w),
    "same number|length"
  )
})

test_that("iv_kitagawa returns iv_test with stat in [0, inf) and p in [0, 1]", {
  set.seed(1)
  for (seed in 1:5) {
    set.seed(seed)
    n <- 200 + seed * 100
    y <- rnorm(n)
    d <- rbinom(n, 1, 0.5)
    z <- rbinom(n, 1, 0.5)
    out <- iv_kitagawa(y, d, z, n_boot = 30, parallel = FALSE)
    expect_true(is.finite(out$statistic))
    expect_gte(out$statistic, 0)
    expect_true(is.finite(out$p_value))
    expect_gte(out$p_value, 0)
    expect_lte(out$p_value, 1)
  }
})

test_that("iv_testjfe rejects fewer than 3 judges", {
  set.seed(1)
  n <- 100
  j <- sample(1:2, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.02 * j)
  y <- rnorm(n, mean = d)
  expect_error(iv_testjfe(y, d, j, n_boot = 10, parallel = FALSE),
               "at least three")
})

test_that("iv_testjfe rejects continuous d beyond max levels", {
  set.seed(1)
  n <- 200
  j <- sample(1:10, n, replace = TRUE)
  d <- rnorm(n)  # continuous
  y <- rnorm(n)
  expect_error(iv_testjfe(y, d, j, n_boot = 10, parallel = FALSE),
               "levels")
})

test_that("iv_mw rejects x with wrong number of rows", {
  set.seed(1)
  n <- 200
  y <- rnorm(n); d <- rbinom(n, 1, 0.5); z <- rbinom(n, 1, 0.5)
  x <- rnorm(50)
  expect_error(iv_mw(y, d, z, x = x, n_boot = 10, parallel = FALSE),
               "rows")
})

test_that("iv_power validates method", {
  set.seed(1)
  n <- 100
  y <- rnorm(n); d <- rbinom(n, 1, 0.5); z <- rbinom(n, 1, 0.5)
  expect_error(
    iv_power(y, d, z, method = "not_a_test", n_sims = 5,
             delta_grid = c(0), n_boot = 10, parallel = FALSE)
  )
})

test_that("iv_check on a non-IV model returns empty", {
  skip_if_not_installed("fixest")
  set.seed(1)
  df <- data.frame(y = rnorm(100), x = rnorm(100))
  m <- fixest::feols(y ~ x, data = df)
  expect_error(iv_check(m, n_boot = 10, parallel = FALSE),
               "not an IV model|Could not")
})
