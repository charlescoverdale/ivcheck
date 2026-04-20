#' Run all applicable IV-validity tests on a fitted model
#'
#' Detects which tests are applicable from the structure of the fitted
#' instrumental variable model and runs them. Returns a tidy summary with
#' a one-line verdict.
#'
#' Applicability is determined by:
#' * Kitagawa (2015) applies to any binary treatment with a discrete instrument.
#' * Mourifie-Wan (2017) applies to the same case, and additionally
#'   supports covariates.
#' * Frandsen-Lefgren-Leslie (2023) applies when the instrument is a set
#'   of mutually exclusive dummy variables (judge-IV / group design).
#'
#' @param model A fitted IV model from [fixest::feols] or `ivreg::ivreg()`.
#' @param tests Character vector of test names to run, or `"all"` (the
#'   default) to run every applicable test.
#' @param alpha Significance level for the verdict. Default 0.05.
#' @param n_boot Number of bootstrap replications. Default 1000.
#' @param ... Further arguments passed to each underlying test.
#'
#' @return An object of class `iv_check` containing a data frame with one
#'   row per test (test name, statistic, p-value, verdict) plus an
#'   overall verdict string.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("fixest", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 500
#'   df <- data.frame(
#'     z = sample(0:1, n, replace = TRUE),
#'     x = rnorm(n)
#'   )
#'   df$d <- rbinom(n, 1, 0.3 + 0.4 * df$z)
#'   df$y <- rnorm(n, mean = df$d + 0.5 * df$x)
#'   m <- fixest::feols(y ~ x | d ~ z, data = df)
#'   iv_check(m, n_boot = 200)
#' }
#' }
#'
#' @export
iv_check <- function(model, tests = "all", alpha = 0.05,
                     n_boot = 1000, ...) {
  # Validate test names before probing the model, so a mistyped name
  # produces a clean "Unknown test name" error even for inputs that
  # cannot be interrogated.
  if (!identical(tests, "all")) {
    unknown <- setdiff(tests, c("kitagawa", "mw", "testjfe"))
    if (length(unknown) > 0L) {
      cli::cli_abort("Unknown test name{?s}: {.val {unknown}}.")
    }
  }

  applicable <- detect_applicable_tests(model)

  if (length(applicable) == 0L) {
    cli::cli_abort(c(
      "Could not detect an applicable IV-validity test for this model.",
      "i" = "The model does not appear to be an IV model, or the treatment is not binary.",
      "i" = "Supported: {.fn fixest::feols} with formula {.code y ~ x | d ~ z} or {.fn ivreg::ivreg}."
    ))
  }

  if (identical(tests, "all")) {
    tests <- applicable
  } else {
    dropped <- setdiff(tests, applicable)
    if (length(dropped) > 0L) {
      cli::cli_warn(
        "Requested test{?s} {.val {dropped}} not applicable to this model; skipping."
      )
      tests <- intersect(tests, applicable)
    }
  }

  results <- list()
  if ("kitagawa" %in% tests) {
    results$kitagawa <- iv_kitagawa(model, n_boot = n_boot,
                                    alpha = alpha, ...)
  }
  if ("mw" %in% tests) {
    # iv_mw conditional path only supports a single covariate in v0.1.0.
    # When iv_check is handed a fitted model with multiple exogenous
    # covariates, run iv_mw in unconditional mode by calling the default
    # method on extracted vectors (the fitted-model dispatch would pull
    # the exogenous x matrix and trigger the multivariate-x abort). The
    # conditional test is available to users who invoke iv_mw() directly
    # with their preferred single covariate.
    used_unconditional <- FALSE
    if (inherits(model, c("fixest", "ivreg"))) {
      yz <- tryCatch(extract_iv_data(model), error = function(e) NULL)
      if (!is.null(yz) && !is.null(yz$x) &&
          is.matrix(yz$x) && ncol(yz$x) > 1L) {
        results$mw <- iv_mw(yz$y, d = yz$d, z = yz$z,
                            n_boot = n_boot, alpha = alpha, ...)
        used_unconditional <- TRUE
      }
    }
    if (!used_unconditional) {
      results$mw <- iv_mw(model, n_boot = n_boot, alpha = alpha, ...)
    }
  }
  if ("testjfe" %in% tests && "testjfe" %in% applicable) {
    results$testjfe <- iv_testjfe(model, n_boot = n_boot,
                                  alpha = alpha, ...)
  }

  tab <- data.frame(
    test      = vapply(results, function(r) r$test, character(1)),
    statistic = vapply(results, function(r) r$statistic, numeric(1)),
    p_value   = vapply(results, function(r) r$p_value, numeric(1)),
    verdict   = vapply(results, function(r) {
      verdict_label(r$p_value, alpha)
    }, character(1)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  overall <- overall_verdict(tab$p_value, alpha)

  structure(
    list(
      table = tab,
      results = results,
      alpha = alpha,
      overall = overall,
      call = sys.call()
    ),
    class = "iv_check"
  )
}

# Internal: decide which tests make sense for a given fitted model.
detect_applicable_tests <- function(model) {
  yz <- tryCatch(extract_iv_data(model), error = function(e) NULL)
  if (is.null(yz)) return(character(0))
  is_binary_d <- length(unique(stats::na.omit(yz$d))) == 2L
  is_binary_z <- !is.null(yz$z) &&
    length(unique(stats::na.omit(yz$z))) == 2L
  is_judge    <- !is.null(yz$z) &&
    length(unique(stats::na.omit(yz$z))) > 2L &&
    is_binary_d
  out <- character(0)
  if (is_binary_d) out <- c(out, "kitagawa", "mw")
  if (is_judge) out <- c(out, "testjfe")
  out
}

verdict_label <- function(p, alpha) {
  if (is.na(p)) return(NA_character_)
  if (p < alpha) "reject" else "pass"
}

overall_verdict <- function(p_values, alpha) {
  if (all(is.na(p_values))) return("inconclusive (no test implemented)")
  any_reject <- any(p_values < alpha, na.rm = TRUE)
  if (any_reject) {
    sprintf("at least one test rejects IV validity at %.2f.", alpha)
  } else {
    sprintf("cannot reject IV validity at %.2f.", alpha)
  }
}
