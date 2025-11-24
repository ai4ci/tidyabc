# ---
# repo: ai4ci/tidyabc
# file: standalone-sigmoid.R
# last-updated: '2025-11-13'
# license: https://unlicense.org
# ---

#' Sigmoid functions
#'
#' @param n a vector of lengths
#' @param n_50 the size at 50 percent
#' @param n_01 the size at 1 percent, and the limit of the returned `p` value
#'   such that \eqn{p \times n \ge n_{01}}
#' @param n_100 an offset under which the value is 1 (default is 0)
#' @param p_max the maximum value at `n <= n_100`
#' @keywords internal
#' @name sigmoid_params
NULL

#' Proportions from a vector of lengths
#'
#' Typically useful in finding alpha values for plotting, or bandwidths for
#' smoothing data this is a decreasing sigmoid with tunable parameters.
#'
#' This is using a composite of a
#' \href{https://en.wikipedia.org/wiki/Hill_equation_(biochemistry)}{Hill function}
#' and a reciprocal:
#'
#' \deqn{
#' f(n) = \frac{1}{1 + \left( \dfrac{n - n_{100}}{n_{50}} \right)^b} \times p_{max}\\
#' b = \frac{\ln(99)}{\ln(n_{01}/n_{50})}
#' }
#'
#' If \eqn{n-n_{100} \le 0} then the maximum value is returned
#'
#' @inheritParams sigmoid_params
#' @returns a vector of proportions for each length
#' @keywords internal
#'
#' @unit
#'
#' testthat::expect_equal(.p_from_n(c(0, 5, 100), 5, 100), c(1, 0.5, 0.01))
#'
#' testthat::expect_error(
#'   {
#'     .p_from_n(0, 50, 1)
#'   },
#'   "inconsistent parameters, n_01 must be larger than n_50.",
#'   fixed = TRUE
#' )
#'
#' testthat::expect_equal(.p_from_n(c(0:10), 5, 10, 0.1), c(
#'   0.1,
#'   0.0999976758268449,
#'   0.0997704297474442,
#'   0.0967278219151741,
#'   0.081446655140766,
#'   0.05,
#'   0.0229935645768342,
#'   0.0097036542590398,
#'   0.00424593238101868,
#'   0.00199056073397113,
#'   0.001
#' ))
#'
#' # does not go below
#' testthat::expect_equal(
#'   .p_from_n(c(1000, 10000, 100000), 5, 50),
#'   c(5e-04, 5e-05, 5e-06)
#' )
#'
.p_from_n = function(n, n_50, n_01, p_max = 1, n_100 = 0) {
  if (n_01 < n_50 || n_50 < n_100) {
    stop("inconsistent parameters, n_01 must be larger than n_50.")
  }

  n = n - n_100
  n_50 = n_50 - n_100
  n_01 = n_01 - n_100
  b = log(99) / (log(n_01) - log(n_50))

  ifelse(n < 0, 1, ifelse(n < n_01, 1 / (1 + (n / n_50)^b), 0.01 * n_01 / n)) *
    p_max
}


#' Simple sigmoid tending towards a minimum proportion of input
#'
#' @inheritParams sigmoid_params
#' @returns a vector of proportions for each length
#' @keywords internal
#'
#' @unit
#' testthat::expect_equal(.p_from_n2(0:5, 1), c(
#'   1,
#'   0.707106781186547,
#'   0.447213595499958,
#'   0.316227766016838,
#'   0.242535625036333,
#'   0.196116135138184
#' ))
.p_from_n2 = function(n, n_inf, p_max = 1, n_100 = 0) {
  n = n - n_100
  return(
    ifelse(n <= 0, 1, n_inf / sqrt(n^2 + n_inf^2)) * p_max
  )
}


#' Scale probabilities by odds ratios
#'
#' \deqn{
#'   O = \frac{p}{1-p} \\
#'   p = \frac{O}{1+O} \\
#'   p_k = \frac{\frac{kp}{1-p}}{( 1 + \frac{kp}{1-p})} \\
#'   p_k = \frac{kp}{1 - p + kp}
#' }
#'
#' @param p A vector of probabilities
#' @param odds_ratio a vector of odds ratios
#'
#' @returns an adjusted vector of probabilities
#' @keywords internal
#' @unit
#' scale_probability(0.5, c(0.25,0.5,1,2,4))
scale_probability = function(p, odds_ratio) {
  # O = p/(1-p)
  # (1-p)O = p; O-Op = p; p = O/(1+O)
  # p_k = kp/(1-p) / ( 1 + kp/(1-p))
  # p_k = kp/(1-p) / (( 1 - p + kp)/(1-p))
  # p_k = kp / ( 1 - p + kp)
  return(odds_ratio * p / (1 - p + odds_ratio * p))
  # changep(c(0.25,0.5,1,2,4),0.5)
}
