test_that("iv_check rejects unknown test names", {
  expect_error(iv_check(NULL, tests = "banana"), "Unknown test name")
})

test_that("verdict_label returns expected strings", {
  expect_identical(ivcheck:::verdict_label(0.01, 0.05), "reject")
  expect_identical(ivcheck:::verdict_label(0.50, 0.05), "pass")
  expect_true(is.na(ivcheck:::verdict_label(NA_real_, 0.05)))
})

test_that("overall_verdict handles all-NA and mixed p-values", {
  expect_match(
    ivcheck:::overall_verdict(c(NA_real_, NA_real_), 0.05),
    "inconclusive"
  )
  expect_match(
    ivcheck:::overall_verdict(c(0.01, 0.5), 0.05),
    "reject"
  )
  expect_match(
    ivcheck:::overall_verdict(c(0.3, 0.5), 0.05),
    "cannot reject"
  )
})
