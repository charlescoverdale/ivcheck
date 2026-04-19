# Bootstrap helpers used by iv_kitagawa, iv_mw, iv_testjfe.
# Not exported.

# Multiplier bootstrap weights.
#
# Rademacher weights used by default; Mammen two-point weights
# available as an option. See Kitagawa (2015) section 3.2 for
# the multiplier bootstrap applied to the KS statistic.
multiplier_weights <- function(n, kind = c("rademacher", "mammen")) {
  kind <- match.arg(kind)
  switch(kind,
    rademacher = sample(c(-1, 1), size = n, replace = TRUE),
    mammen = {
      phi <- (1 + sqrt(5)) / 2
      w1 <- -(sqrt(5) - 1) / 2
      w2 <-  (sqrt(5) + 1) / 2
      p1 <- (sqrt(5) + 1) / (2 * sqrt(5))
      u  <- stats::runif(n)
      ifelse(u < p1, w1, w2)
    }
  )
}

# Run B bootstrap replications of a statistic fn(data, weights).
# Parallelised on POSIX; serial on Windows.
run_bootstrap <- function(fn, n, B, parallel = TRUE, mc_cores = NULL, ...) {
  do_one <- function(b) {
    w <- multiplier_weights(n)
    fn(w, ...)
  }
  if (parallel && .Platform$OS.type == "unix") {
    cores <- mc_cores %||% max(1L, parallel::detectCores() - 1L)
    out <- parallel::mclapply(seq_len(B), do_one, mc.cores = cores)
  } else {
    out <- lapply(seq_len(B), do_one)
  }
  unlist(out, use.names = FALSE)
}
