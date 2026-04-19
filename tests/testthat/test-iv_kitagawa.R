test_that("iv_kitagawa returns an iv_test object with the right structure", {
  set.seed(1)
  n <- 200
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d)
  out <- iv_kitagawa(y, d, z, n_boot = 50)
  expect_s3_class(out, "iv_test")
  expect_identical(out$test, "Kitagawa (2015)")
  expect_equal(out$n, n)
  expect_equal(length(out$boot_stats), 50L)
})

test_that("iv_kitagawa rejects non-binary d with a clear message", {
  y <- rnorm(100)
  d <- rnorm(100)
  z <- sample(0:1, 100, replace = TRUE)
  expect_error(iv_kitagawa(y, d, z, n_boot = 10), "binary")
})

test_that("iv_kitagawa rejects unequal-length inputs", {
  y <- rnorm(100)
  d <- rbinom(50, 1, 0.5)
  z <- sample(0:1, 100, replace = TRUE)
  expect_error(iv_kitagawa(y, d, z, n_boot = 10), "unequal lengths")
})
