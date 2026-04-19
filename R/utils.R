# Internal helpers (not exported)

`%||%` <- function(a, b) if (is.null(a)) b else a

iv_env <- new.env(parent = emptyenv())

validate_binary <- function(d, arg = "d") {
  if (!is.numeric(d) && !is.logical(d) && !is.integer(d)) {
    cli::cli_abort("{.arg {arg}} must be numeric, integer, or logical.")
  }
  if (anyNA(d)) {
    cli::cli_abort("{.arg {arg}} contains {.val NA} values.")
  }
  vals <- unique(as.numeric(d))
  if (length(vals) != 2L || !all(vals %in% c(0, 1))) {
    cli::cli_abort(
      "{.arg {arg}} must be a binary 0/1 vector. Found values: {.val {vals}}."
    )
  }
  invisible(as.numeric(d))
}

validate_discrete <- function(z, arg = "z", max_levels = 50L) {
  if (is.factor(z)) {
    z <- as.integer(z)
  }
  if (!is.numeric(z) && !is.integer(z)) {
    cli::cli_abort("{.arg {arg}} must be numeric, integer, or factor.")
  }
  if (anyNA(z)) {
    cli::cli_abort("{.arg {arg}} contains {.val NA} values.")
  }
  k <- length(unique(z))
  if (k < 2L) {
    cli::cli_abort("{.arg {arg}} must take at least two distinct values.")
  }
  if (k > max_levels) {
    cli::cli_abort(
      "{.arg {arg}} has {k} levels; this test expects a discrete instrument."
    )
  }
  invisible(as.numeric(z))
}

validate_numeric <- function(x, arg = "x", allow_na = FALSE) {
  if (!is.numeric(x)) {
    cli::cli_abort("{.arg {arg}} must be a numeric vector.")
  }
  if (!allow_na && anyNA(x)) {
    cli::cli_abort("{.arg {arg}} contains {.val NA} values.")
  }
  invisible(x)
}

check_lengths <- function(...) {
  args <- list(...)
  lens <- vapply(args, length, integer(1))
  if (length(unique(lens)) != 1L) {
    cli::cli_abort(
      "Input vectors have unequal lengths: {.val {lens}}."
    )
  }
  invisible(lens[1])
}
