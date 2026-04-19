#' Plot method for an IV-validity test
#'
#' Plots the bootstrap distribution of the test statistic with the
#' observed statistic and the rejection region highlighted.
#'
#' @param x An object of class `iv_test`.
#' @param ... Further graphical arguments passed to [graphics::hist].
#' @return Invisibly returns `x`.
#' @export
plot.iv_test <- function(x, ...) {
  if (is.null(x$boot_stats) || all(is.na(x$boot_stats))) {
    cli::cli_warn(
      "No bootstrap statistics stored; cannot plot. This is expected for the v0.1.0 stub."
    )
    return(invisible(x))
  }
  crit <- stats::quantile(x$boot_stats, 1 - x$alpha, na.rm = TRUE)
  graphics::hist(
    x$boot_stats,
    main = sprintf("%s: bootstrap distribution", x$test),
    xlab = "Test statistic",
    ...
  )
  graphics::abline(v = x$statistic, lwd = 2)
  graphics::abline(v = crit, lty = 2)
  invisible(x)
}
