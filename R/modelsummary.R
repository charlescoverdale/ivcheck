# modelsummary integration. Registered in .onLoad via zzz.R only when
# the modelsummary package is installed. No hard dependency.

glance_iv_test <- function(x, ...) {
  data.frame(
    test = x$test,
    statistic = x$statistic %||% NA_real_,
    p.value = x$p_value %||% NA_real_,
    n = x$n %||% NA_integer_
  )
}

glance_iv_check <- function(x, ...) {
  x$table
}
