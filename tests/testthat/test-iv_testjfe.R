test_that("iv_testjfe returns an iv_test object with the right structure", {
  set.seed(1)
  n <- 500
  judge <- sample.int(10, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.02 * judge)
  y <- rnorm(n, mean = d)
  out <- iv_testjfe(y, d, judge, n_boot = 50)
  expect_s3_class(out, "iv_test")
  expect_identical(out$test, "Frandsen-Lefgren-Leslie (2023)")
  expect_equal(out$n_judges, 10L)
})
