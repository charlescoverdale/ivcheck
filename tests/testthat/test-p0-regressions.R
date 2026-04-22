# Regression tests for the three P0 behaviours flagged in DEEP_REVIEW.md:
#   P0.1  fixest / ivreg dispatch must abort on fitted FE-IV models.
#   P0.2  iv_testjfe on multivalued D with within-judge-varying x must use
#         FWL-residualised indicators in the structural residual; if it does,
#         the asymptotic null size stays within a couple MC SEs of nominal.
#   P0.3  binding-direction labels must reference the user's original D
#         coding (both numeric levels like {10, 20, 30} and factor levels
#         like c("low", "med", "high")), not the internal 0..(k-1) remap.

test_that("P0.1: iv_kitagawa.fixest aborts on models with fixed effects", {
  skip_if_not_installed("fixest")
  set.seed(1)
  n <- 400
  firm <- sample.int(15, n, replace = TRUE)
  z <- sample(0:1, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.4 * z)
  y <- rnorm(n, mean = d + 0.3 * firm)
  df <- data.frame(y = y, d = d, z = z, firm = factor(firm))
  m <- fixest::feols(y ~ 1 | firm | d ~ z, data = df)
  expect_error(
    iv_kitagawa(m, n_boot = 20, parallel = FALSE),
    "does not support models with fixed effects"
  )
})

test_that("P0.1: iv_testjfe.fixest aborts on models with fixed effects", {
  skip_if_not_installed("fixest")
  set.seed(2)
  n <- 400
  firm <- sample.int(15, n, replace = TRUE)
  z <- sample.int(12, n, replace = TRUE)
  d <- rbinom(n, 1, 0.3 + 0.02 * z)
  y <- rnorm(n, mean = d)
  df <- data.frame(y = y, d = d, z = factor(z), firm = factor(firm))
  m <- fixest::feols(y ~ 1 | firm | d ~ z, data = df)
  expect_error(
    iv_testjfe(m, n_boot = 20, parallel = FALSE),
    "does not support models with fixed effects"
  )
})

test_that("P0.2: iv_testjfe M = 2 with within-judge-varying x keeps size near nominal", {
  skip_on_cran()
  set.seed(9001)
  K <- 15
  n <- 3000
  n_mc <- 120
  pvals <- numeric(n_mc)
  for (s in seq_len(n_mc)) {
    judge <- sample.int(K, n, replace = TRUE)
    x <- rnorm(n)
    p1 <- plogis(-0.5 + 0.05 * judge)
    p2 <- plogis(-1.0 + 0.03 * judge + 0.3 * x)
    u <- runif(n)
    d <- as.integer(u > p1) + as.integer(u > p1 + (1 - p1) * p2)
    y <- 0.5 * d + 0.2 * x + rnorm(n)
    pvals[s] <- iv_testjfe(y, d, judge, x = x, n_boot = 1,
                           parallel = FALSE)$p_value_asymptotic
  }
  size5 <- mean(pvals < 0.05, na.rm = TRUE)
  # MC SE at nominal 5% with 120 draws is sqrt(.05*.95/120) ~ 0.0199, so +/- 2 SE
  # is about [0.01, 0.09]. The bug inflates sigma^2 and pushes size well below
  # this band (often < 0.01); the fix restores it.
  expect_gt(size5, 0.005)
  expect_lt(size5, 0.12)
})

test_that("P0.3: binding$direction preserves numeric D values", {
  set.seed(101)
  n <- 800
  z <- sample(0:2, n, replace = TRUE)
  d_user <- sample(c(10L, 20L, 30L), n, replace = TRUE,
                   prob = c(0.3, 0.4, 0.3))
  d_user[z == 2] <- sample(c(20L, 30L), sum(z == 2), replace = TRUE,
                           prob = c(0.3, 0.7))
  y <- rnorm(n, mean = d_user * 0.02)
  r <- iv_kitagawa(y, d_user, z, n_boot = 30, parallel = FALSE)
  # Binding direction string must reference one of the user's actual D values
  # (10, 20, or 30), not the internal 0/1/2 remap.
  expect_match(r$binding$direction, "10|20|30")
  expect_false(grepl("(^|[^0-9])0([^0-9]|$)", r$binding$direction))
})

test_that("P0.3: binding$direction preserves factor D labels", {
  set.seed(202)
  n <- 800
  z <- sample(0:2, n, replace = TRUE)
  d_fac <- factor(sample(c("low", "med", "high"), n, replace = TRUE,
                         prob = c(0.3, 0.4, 0.3)),
                  levels = c("low", "med", "high"))
  d_fac[z == 2] <- factor(sample(c("med", "high"), sum(z == 2), replace = TRUE,
                                 prob = c(0.3, 0.7)),
                          levels = c("low", "med", "high"))
  y <- rnorm(n, mean = as.integer(d_fac) * 0.5)
  r <- iv_kitagawa(y, d_fac, z, n_boot = 30, parallel = FALSE)
  expect_match(r$binding$direction, "low|med|high")
})

test_that("P0.3: binding$direction preserves binary D values in {5, 7}", {
  set.seed(303)
  n <- 600
  z <- sample(0:1, n, replace = TRUE)
  d_bin <- sample(c(5L, 7L), n, replace = TRUE, prob = c(0.4, 0.6))
  d_bin[z == 1] <- sample(c(5L, 7L), sum(z == 1), replace = TRUE,
                          prob = c(0.2, 0.8))
  y <- rnorm(n, mean = d_bin * 0.3)
  r <- iv_kitagawa(y, d_bin, z, n_boot = 30, parallel = FALSE)
  expect_match(r$binding$direction, "5|7")
})
