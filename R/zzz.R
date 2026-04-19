.onLoad <- function(libname, pkgname) {
  # Conditional registration of broom::glance methods for modelsummary.
  # If broom or modelsummary is not installed, this is a no-op.
  if (requireNamespace("broom", quietly = TRUE)) {
    registerS3method("glance", "iv_test",  glance_iv_test,  envir = asNamespace("broom"))
    registerS3method("glance", "iv_check", glance_iv_check, envir = asNamespace("broom"))
  }
  invisible()
}
