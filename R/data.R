#' Card (1995) proximity-to-college extract
#'
#' A data extract from the National Longitudinal Survey of Young Men,
#' as used in Card (1995) to estimate the return to schooling using
#' proximity to a four-year college as an instrument for years of
#' schooling. The extract adds a binary `college` indicator (16+ years
#' of schooling) so the data can be used with IV-validity tests that
#' require a binary treatment.
#'
#' @format A data frame with 2991 rows and 11 variables:
#' \describe{
#'   \item{id}{Integer row identifier.}
#'   \item{lwage}{Log hourly wage in 1976 (outcome in Card's specification).}
#'   \item{educ}{Years of completed schooling (continuous; Card's
#'     endogenous regressor).}
#'   \item{college}{Integer 0/1 indicator for `educ >= 16`. Use this
#'     when a test requires a binary treatment.}
#'   \item{near_college}{Integer 0/1 indicator for growing up near a
#'     four-year college (Card's instrument).}
#'   \item{age}{Age in 1976.}
#'   \item{exper}{Years of potential labour-market experience (age minus
#'     schooling minus six).}
#'   \item{black}{Integer 0/1 indicator for black respondents.}
#'   \item{south}{Integer 0/1 indicator for residence in the US south.}
#'   \item{smsa}{Integer 0/1 indicator for residence in a Standard
#'     Metropolitan Statistical Area.}
#'   \item{married}{Integer 0/1 indicator for married respondents.}
#' }
#'
#' @source Card, D. (1995). Using Geographic Variation in College
#'   Proximity to Estimate the Return to Schooling. In *Aspects of
#'   Labour Market Behaviour: Essays in Honour of John Vanderkamp*, ed.
#'   L. N. Christofides, E. K. Grant, and R. Swidinsky, 201-222.
#'   University of Toronto Press. Original data from the 1966-1976
#'   National Longitudinal Survey of Young Men. Cleaned extract via
#'   the `wooldridge` package on CRAN.
#'
#' @references
#' Card, D. (1995). Using Geographic Variation in College Proximity to
#' Estimate the Return to Schooling. In Christofides, Grant, and
#' Swidinsky (eds.), *Aspects of Labour Market Behaviour: Essays in
#' Honour of John Vanderkamp*, 201-222.
#'
#' Wooldridge, J. M. (2020). *wooldridge: 115 Data Sets from
#' "Introductory Econometrics: A Modern Approach"*. R package.
#'
#' @examples
#' data(card1995)
#' summary(card1995$lwage)
#' table(near_college = card1995$near_college,
#'       college      = card1995$college)
"card1995"
