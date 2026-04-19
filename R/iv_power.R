#' Monte Carlo power curve for IV-validity tests
#'
#' Simulates data under a user-specified deviation from validity and
#' estimates the rejection probability of the chosen test at each
#' deviation size. Useful for sample-size planning and for benchmarking
#' different tests on the same DGP.
#'
#' The deviation is parametrised as the share of the population that
#' violates monotonicity (defiers) relative to the share of compliers.
#' A deviation of 0 corresponds to the null; larger values correspond to
#' greater violations.
#'
#' @param y,d,z Observed data used to anchor the DGP (marginal
#'   distributions, cell counts).
#' @param method Which test to benchmark. One of `"kitagawa"`, `"mw"`,
#'   `"testjfe"`.
#' @param alpha Significance level.
#' @param n_sims Number of Monte Carlo simulations per deviation.
#' @param delta_grid Numeric vector of deviation sizes to evaluate.
#'   If `NULL`, defaults to `seq(0, 0.3, by = 0.05)`.
#' @param ... Further arguments passed to the underlying test.
#'
#' @return A data frame with columns `delta` and `power` (estimated
#'   rejection probability).
#'
#' @examples
#' \dontrun{
#' # Headline power curve for a small-N design
#' set.seed(1)
#' y <- rnorm(300); d <- rbinom(300, 1, 0.5); z <- sample(0:1, 300, TRUE)
#' iv_power(y, d, z, method = "kitagawa", n_sims = 50)
#' }
#'
#' @export
iv_power <- function(y, d, z, method = c("kitagawa", "mw", "testjfe"),
                     alpha = 0.05, n_sims = 500, delta_grid = NULL, ...) {
  method <- match.arg(method)
  if (is.null(delta_grid)) {
    delta_grid <- seq(0, 0.3, by = 0.05)
  }
  # TODO v0.1.0: implement DGP perturbation + Monte Carlo loop.
  data.frame(
    delta = delta_grid,
    power = rep(NA_real_, length(delta_grid))
  )
}
