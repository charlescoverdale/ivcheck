# Shared helpers for extracting (y, d, z, x) from fitted IV model objects.
# Internal only; S3 dispatch within the package namespace.

#' @noRd
extract_iv_data <- function(object) {
  UseMethod("extract_iv_data")
}

#' @noRd
extract_iv_data.default <- function(object) {
  cli::cli_abort(
    c(
      "No {.fn extract_iv_data} method for objects of class {.cls {class(object)[1]}}.",
      i = "Supported classes: {.cls fixest}, {.cls ivreg}. Or use the default method with raw vectors."
    )
  )
}

#' @noRd
extract_iv_data.fixest <- function(object) {
  if (!requireNamespace("fixest", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg fixest} is required for this method.")
  }
  if (!isTRUE(object$iv)) {
    cli::cli_abort(
      "Model is not an IV model. Use {.code feols(y ~ x | d ~ z, data)} syntax."
    )
  }
  # Placeholder extraction. Full implementation wired up once we have
  # a test case; fixest exposes endogenous/instrument model frames via
  # its accessor methods.
  list(
    y = NULL,
    d = NULL,
    z = NULL,
    x = NULL
  )
}

#' @noRd
extract_iv_data.ivreg <- function(object) {
  if (!requireNamespace("ivreg", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ivreg} is required for this method.")
  }
  # Placeholder extraction.
  list(
    y = NULL,
    d = NULL,
    z = NULL,
    x = NULL
  )
}
