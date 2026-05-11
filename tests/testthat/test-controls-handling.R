# Tests for v0.1.2 behaviour change: iv_kitagawa() refuses fitted models
# with exogenous controls, and iv_check() filters its applicable-tests
# list to reflect that.

test_that("iv_kitagawa errors on a fixest model with one exogenous control", {
  skip_if_not_installed("fixest")
  set.seed(1)
  n <- 300
  df <- data.frame(
    z = sample(0:1, n, replace = TRUE),
    x = rnorm(n)
  )
  df$d <- rbinom(n, 1, 0.3 + 0.4 * df$z)
  df$y <- rnorm(n, mean = df$d + 0.2 * df$x)
  m <- fixest::feols(y ~ x | d ~ z, data = df)
  expect_error(
    iv_kitagawa(m, n_boot = 20, parallel = FALSE),
    "unconditional"
  )
})

test_that("iv_kitagawa errors on a fixest model with multiple exogenous controls", {
  skip_if_not_installed("fixest")
  set.seed(1)
  n <- 300
  df <- data.frame(
    z = sample(0:1, n, replace = TRUE),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  df$d <- rbinom(n, 1, 0.3 + 0.4 * df$z)
  df$y <- rnorm(n, mean = df$d + 0.2 * df$x1 + 0.1 * df$x2)
  m <- fixest::feols(y ~ x1 + x2 | d ~ z, data = df)
  expect_error(
    iv_kitagawa(m, n_boot = 20, parallel = FALSE),
    "unconditional"
  )
})

test_that("iv_kitagawa errors on an ivreg model with exogenous controls", {
  skip_if_not_installed("ivreg")
  set.seed(1)
  n <- 300
  df <- data.frame(
    z = sample(0:1, n, replace = TRUE),
    x = rnorm(n)
  )
  df$d <- rbinom(n, 1, 0.3 + 0.4 * df$z)
  df$y <- rnorm(n, mean = df$d + 0.2 * df$x)
  m <- ivreg::ivreg(y ~ x + d | x + z, data = df)
  expect_error(
    iv_kitagawa(m, n_boot = 20, parallel = FALSE),
    "unconditional"
  )
})

test_that("iv_kitagawa raw-vector path is unchanged regardless of separate X data", {
  set.seed(1)
  n <- 300
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  out <- iv_kitagawa(y, d, z, n_boot = 50, parallel = FALSE)
  expect_s3_class(out, "iv_test")
  expect_identical(out$test, "Kitagawa (2015)")
})

test_that("iv_check skips Kitagawa when fixest model has one exogenous control", {
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
  chk <- NULL
  expect_message(
    chk <- iv_check(m, n_boot = 30, parallel = FALSE),
    "Kitagawa test skipped"
  )
  expect_s3_class(chk, "iv_check")
  expect_false("Kitagawa (2015)" %in% chk$table$test)
  expect_true(any(grepl("Mourifie", chk$table$test)))
})

test_that("iv_check skips both Kitagawa and MW when fixest model has multivariate X", {
  skip_if_not_installed("fixest")
  set.seed(1)
  n <- 400
  df <- data.frame(
    z = sample(0:1, n, replace = TRUE),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  df$d <- rbinom(n, 1, 0.3 + 0.4 * df$z)
  df$y <- rnorm(n, mean = df$d + 0.2 * df$x1 + 0.1 * df$x2)
  m <- fixest::feols(y ~ x1 + x2 | d ~ z, data = df)
  # Both Kitagawa and MW should be skipped, and there's no FLL applicable
  # (binary z), so iv_check has nothing to run and aborts.
  expect_error(
    iv_check(m, n_boot = 30, parallel = FALSE),
    "Could not detect an applicable"
  )
})

test_that("iv_check is silent (no informational messages) when no X is present", {
  skip_if_not_installed("fixest")
  set.seed(1)
  n <- 400
  df <- data.frame(z = sample(0:1, n, replace = TRUE))
  df$d <- rbinom(n, 1, 0.3 + 0.4 * df$z)
  df$y <- rnorm(n, mean = df$d)
  m <- fixest::feols(y ~ 1 | d ~ z, data = df)
  expect_no_message(
    iv_check(m, n_boot = 30, parallel = FALSE),
    message = "skipped"
  )
})

test_that("detect_applicable_tests filters Kitagawa on X presence and MW on multivariate X", {
  yz_no_x   <- list(y = rnorm(100), d = rbinom(100, 1, 0.5),
                    z = sample(0:1, 100, replace = TRUE), x = NULL)
  yz_one_x  <- list(y = rnorm(100), d = rbinom(100, 1, 0.5),
                    z = sample(0:1, 100, replace = TRUE),
                    x = matrix(rnorm(100), ncol = 1L))
  yz_two_x  <- list(y = rnorm(100), d = rbinom(100, 1, 0.5),
                    z = sample(0:1, 100, replace = TRUE),
                    x = matrix(rnorm(200), ncol = 2L))
  yz_judge  <- list(y = rnorm(100), d = rbinom(100, 1, 0.5),
                    z = sample(1:5, 100, replace = TRUE),
                    x = NULL)
  expect_setequal(
    ivcheck:::detect_applicable_tests(NULL, yz = yz_no_x),
    c("kitagawa", "mw")
  )
  expect_setequal(
    ivcheck:::detect_applicable_tests(NULL, yz = yz_one_x),
    c("mw")
  )
  expect_identical(
    ivcheck:::detect_applicable_tests(NULL, yz = yz_two_x),
    character(0)
  )
  expect_setequal(
    ivcheck:::detect_applicable_tests(NULL, yz = yz_judge),
    c("kitagawa", "mw", "testjfe")
  )
})
