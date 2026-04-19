# Build the bundled card1995 dataset.
#
# Source: Card, D. (1995). "Using Geographic Variation in College Proximity
# to Estimate the Return to Schooling." In Aspects of Labour Market
# Behaviour: Essays in Honour of John Vanderkamp, ed. L. N. Christofides,
# E. K. Grant, and R. Swidinsky. Toronto: University of Toronto Press,
# pp. 201-222.
#
# The raw data are from the National Longitudinal Survey of Young Men,
# 1966-1976. We obtain the cleaned extract from the `wooldridge` R
# package (Wooldridge, J. M., 2020, `wooldridge: 115 Data Sets from
# "Introductory Econometrics: A Modern Approach"`, CRAN).
#
# We keep the subset of columns needed for IV-validity illustrations and
# add a binary `college` indicator (educ >= 16) for use with tests that
# require a binary treatment, notably Kitagawa (2015).

library(wooldridge)

raw <- wooldridge::card

card1995 <- data.frame(
  id           = seq_len(nrow(raw)),
  lwage        = raw$lwage,
  educ         = raw$educ,
  college      = as.integer(raw$educ >= 16),
  near_college = as.integer(raw$nearc4),
  age          = raw$age,
  exper        = raw$exper,
  black        = as.integer(raw$black),
  south        = as.integer(raw$south),
  smsa         = as.integer(raw$smsa),
  married      = as.integer(raw$married)
)

# Drop rows with any NA in the retained columns
card1995 <- card1995[stats::complete.cases(card1995), ]
rownames(card1995) <- NULL

# Sanity checks
stopifnot(
  nrow(card1995) > 2900,
  all(card1995$near_college %in% c(0, 1)),
  all(card1995$college %in% c(0, 1)),
  all(!is.na(card1995$lwage))
)

# Save to data/ (use version 2 for R >= 3.5 compatibility per CRAN)
usethis::use_data(card1995, overwrite = TRUE, version = 2, compress = "xz")
