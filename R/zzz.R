.onLoad <- function(libname, pkgname) {
  # Internal S3 method registration. These are not @export'd because
  # extract_iv_data is internal; registerS3method makes S3 dispatch find
  # them at runtime.
  registerS3method("extract_iv_data", "default", extract_iv_data.default)
  registerS3method("extract_iv_data", "fixest",  extract_iv_data.fixest)
  registerS3method("extract_iv_data", "ivreg",   extract_iv_data.ivreg)

  # Conditional registration of broom::glance methods for modelsummary.
  if (requireNamespace("broom", quietly = TRUE)) {
    registerS3method("glance", "iv_test",  glance_iv_test,  envir = asNamespace("broom"))
    registerS3method("glance", "iv_check", glance_iv_check, envir = asNamespace("broom"))
  }
  invisible()
}
