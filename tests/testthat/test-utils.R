test_that("validate_binary accepts 0/1 vectors and rejects others", {
  expect_invisible(ivcheck:::validate_binary(c(0, 1, 1, 0)))
  expect_error(ivcheck:::validate_binary(c(0.5, 1)), "binary")
  expect_error(ivcheck:::validate_binary(c(1, 2, 3)), "binary")
})

test_that("validate_binary rejects NA", {
  expect_error(ivcheck:::validate_binary(c(0, 1, NA, 1)), "NA")
})

test_that("validate_discrete enforces minimum two levels", {
  expect_error(ivcheck:::validate_discrete(rep(1, 10)), "two distinct")
})

test_that("validate_discrete rejects NA", {
  expect_error(ivcheck:::validate_discrete(c(0, 1, NA, 2)), "NA")
})

test_that("check_lengths flags unequal lengths", {
  expect_error(
    ivcheck:::check_lengths(1:3, 1:5),
    "unequal lengths"
  )
})

test_that("%||% returns fallback on NULL", {
  expect_identical(ivcheck:::`%||%`(NULL, 7), 7)
  expect_identical(ivcheck:::`%||%`(5, 7), 5)
})
