#' Print method for an IV-validity test
#'
#' @param x An object of class `iv_test`.
#' @param digits Number of significant digits to display.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.iv_test <- function(x, digits = 3L, ...) {
  cli::cli_h1(x$test)
  cli::cli_text("Sample size: {.val {x$n}}")
  cli::cli_text(
    "Statistic: {.val {format(x$statistic, digits = digits)}}, ",
    "p-value: {.val {format(x$p_value, digits = digits)}}"
  )
  verdict <- if (is.na(x$p_value)) {
    "not yet implemented (v0.1.0 stub)"
  } else if (x$p_value < x$alpha) {
    sprintf("reject IV validity at %.2f", x$alpha)
  } else {
    sprintf("cannot reject IV validity at %.2f", x$alpha)
  }
  cli::cli_text("Verdict: {verdict}")
  invisible(x)
}

#' Summary method for an IV-validity test
#'
#' @param object An object of class `iv_test`.
#' @param ... Ignored.
#' @return Invisibly returns `object`.
#' @export
summary.iv_test <- function(object, ...) {
  print(object)
  if (!is.null(object$binding)) {
    cli::cli_text("Binding pair: {.val {object$binding}}")
  }
  if (!is.null(object$n_judges)) {
    cli::cli_text("Number of judges: {.val {object$n_judges}}")
  }
  cli::cli_text("Bootstrap replications: {.val {object$n_boot}}")
  invisible(object)
}

#' Format method for an IV-validity test
#'
#' Used when an `iv_test` object is included as a column of a data frame
#' or tibble.
#'
#' @param x An object of class `iv_test`.
#' @param ... Ignored.
#' @return A one-line character summary.
#' @export
format.iv_test <- function(x, ...) {
  sprintf(
    "%s: stat = %.3f, p = %.3f",
    x$test, x$statistic %||% NA_real_, x$p_value %||% NA_real_
  )
}

#' Print method for an iv_check result
#'
#' @param x An object of class `iv_check`.
#' @param digits Number of significant digits.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.iv_check <- function(x, digits = 3L, ...) {
  cli::cli_h1("IV validity diagnostic")
  tab <- x$table
  for (i in seq_len(nrow(tab))) {
    cli::cli_text(
      "  {.strong {tab$test[i]}}: ",
      "stat = {.val {format(tab$statistic[i], digits = digits)}}, ",
      "p = {.val {format(tab$p_value[i], digits = digits)}}, ",
      "{tab$verdict[i]}"
    )
  }
  cli::cli_text("Overall: {x$overall}")
  invisible(x)
}
