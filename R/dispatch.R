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
      c(
        "Model is not an IV model.",
        i = "Use {.code feols(y ~ x | d ~ z, data)} syntax to fit an IV model."
      )
    )
  }
  if (length(object$fixef_vars) > 0L) {
    cli::cli_abort(c(
      "{.fn ivcheck} v0.1.0 does not support models with fixed effects.",
      i = "The discrete-Z tests of Kitagawa / Mourifie-Wan / FLL operate on the raw (Y, D, Z) joint distribution; residualising on fixed effects destroys the discrete structure of Z.",
      i = "Workaround: pre-demean {.var Y} and {.var D} within each FE cell, keep {.var Z} discrete, and call the default method on the resulting vectors."
    ))
  }
  endo_names <- object$iv_endo_names
  inst_names <- object$iv_inst_names
  if (length(endo_names) != 1L) {
    cli::cli_abort(
      "{.fn ivcheck} v0.1.0 supports exactly one endogenous treatment; model has {length(endo_names)}."
    )
  }
  # Recover y, d, z without relying on the call environment.
  y <- as.numeric(stats::model.matrix(object, type = "lhs"))
  fs <- object$iv_first_stage[[1]]
  d <- as.numeric(stats::fitted(fs) + stats::residuals(fs))
  mm_fs <- stats::model.matrix(fs)
  # The first stage's model.matrix has exogenous X columns plus the
  # instruments. Extract the instruments by name.
  z_cols <- intersect(inst_names, colnames(mm_fs))
  if (length(z_cols) == 0L) {
    cli::cli_abort(
      "Could not locate instrument {.val {inst_names}} in the first-stage design matrix."
    )
  }
  z <- if (length(z_cols) == 1L) mm_fs[, z_cols] else mm_fs[, z_cols]
  x_cols <- setdiff(colnames(mm_fs), c("(Intercept)", inst_names))
  list(
    y = y,
    d = d,
    z = z,
    x = if (length(x_cols) > 0L) mm_fs[, x_cols, drop = FALSE] else NULL
  )
}

#' @noRd
extract_iv_data.ivreg <- function(object) {
  if (!requireNamespace("ivreg", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ivreg} is required for this method.")
  }
  endo_names <- names(object$endogenous)
  inst_names <- names(object$instruments)
  exog_names <- setdiff(names(object$exogenous), "(Intercept)")
  if (length(endo_names) != 1L) {
    cli::cli_abort(
      "{.fn ivcheck} v0.1.0 supports exactly one endogenous treatment; model has {length(endo_names)}."
    )
  }
  mf <- object$model
  list(
    y = as.numeric(object$y),
    d = mf[[endo_names]],
    z = if (length(inst_names) == 1L) mf[[inst_names]] else mf[, inst_names],
    x = if (length(exog_names) > 0L) as.matrix(mf[, exog_names, drop = FALSE]) else NULL
  )
}
