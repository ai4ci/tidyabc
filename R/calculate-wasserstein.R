#' Calculate a wasserstein distance
#'
#' This function takes simulation and observed data and calculates a normalised
#' wasserstein distance.
#'
#' In the comparison unequal lengths of the data can be accommodated. The
#' simulated data is recycled, and sampled, until the same length as the
#' observed data before the comparison.
#'
#' @param sim A vector of simulated data points
#' @param obs A vector of observed data points
#' @param debias Should the simulations be shifted to match the mean of the
#'   observed data
#'
#' @returns a length normalised wasserstein distance. This is the average distance an
#'   individual simulated data point must be shifted to match the observed data normalised
#'   by the average distance of the observed data from the mean.
#' @export
#' @concept workflow
#'
#' @unit
#'
#' testthat::expect_equal(calculate_wasserstein(0:10, 10:0), 0)
#'
#' # zero if no distance
#' testthat::expect_equal(calculate_wasserstein(0:10, 0:10), 0)
#' testthat::expect_equal(calculate_wasserstein(10:0, 0:10), 0)
#'
#' # normalised so that all mass at mean = 1
#' testthat::expect_equal(calculate_wasserstein(rep(5, 11), 0:10), 1)
#'
#' # smaller sample recycled and normalises to same value
#' testthat::expect_equal(calculate_wasserstein(rep(5, 5), 0:10), 1)
#'
#' # should be ((0+1+0+1+0+0+0+1+0+1+0) / 11) / ((5+4+3+2+1+0+1+2+3+4+5) / 11) = 0.1333...
#' testthat::expect_equal(
#'   calculate_wasserstein(
#'     c(0, 0, 2, 2, 4, 5, 6, 8, 8, 10, 10),
#'     0:10
#'  ),
#'   0.133333333333333
#' )
#'
#' withr::with_seed(100, {
#'    ref = rnorm(1000)
#'    cmp1 = rnorm(1000)
#'    cmp2 = rnorm(1000, sd=2)
#'    cmp3 = rnorm(100)
#' })
#'
#' testthat::expect_equal(calculate_wasserstein(cmp1,ref), 0.0576417503974498)
#' testthat::expect_equal(calculate_wasserstein(cmp2,ref), 1.04950621760385)
#' testthat::expect_equal(calculate_wasserstein(cmp3,ref), 0.212977231452674)
#'
calculate_wasserstein = function(sim, obs, debias = FALSE) {
  obs = sort(obs)
  av = mean(obs)
  # What is the average distance of points in the reference data from the mean?
  # we use this to normalise the result
  norm = mean(abs(obs - av))
  size = length(obs)

  sim = .match_size(sim, size)
  sim = sort(sim)
  if (debias) {
    sim = sim - mean(sim) + av
  }
  return(mean(abs(sim - obs)) / norm)
}


#' Generate a function to calculate a wasserstein distance
#'
#' `r lifecycle::badge('experimental')`
#'
#' This function takes reference data in the forms of individual observations in
#' the form of for example event times and returns a crated function to
#' calculate the wasserstein distance from simulated data to the observed data.
#' For large simulations this will be quicker than `calculate_wasserstein` but
#' must be used with a crated `scorer_fn`
#'
#' In the comparison unequal lengths of the data can be accommodated. The
#' simulated data is recycled, and sampled, until the same length as the
#' observed data before the comparison.
#'
#' @param obs A vector of observations
#' @param debias Should the simulations be shifted to match the mean of the
#'   observed data
#'
#' @returns a function that takes parameter `sim` and returns a
#'   length normalised wasserstein distance. This is the average distance an
#'   individual data point must be shifted to match the reference data normalised
#'   by the average distance of the reference data from the mean. The function
#'   also takes a `p` parameter which can be a `progressr` progress bar which must
#'   be named.
#' @export
#'
#' @unit
#' tmp = wasserstein_calculator(0:10)
#'
#' # zero if no distance
#' testthat::expect_equal(tmp(0:10), 0)
#' testthat::expect_equal(tmp(10:0), 0)
#'
#' # normalised so that all mass at mean = 1
#' testthat::expect_equal(tmp(rep(5, 11)), 1)
#'
#' # smaller sample recycled and normalises to same value
#' testthat::expect_equal(tmp(rep(5, 5)), 1)
#'
#' # should be ((0+1+0+1+0+0+0+1+0+1+0) / 11) / ((5+4+3+2+1+0+1+2+3+4+5) / 11) = 0.1333...
#' testthat::expect_equal(
#'   tmp(c(0, 0, 2, 2, 4, 5, 6, 8, 8, 10, 10)),
#'   0.133333333333333
#' )
#'
#' withr::with_seed(100, {
#'    ref = rnorm(1000)
#'    cmp1 = rnorm(1000)
#'    cmp2 = rnorm(1000, sd=2)
#'    cmp3 = rnorm(100)
#' })
#'
#' tmp2 = wasserstein_calculator(ref)
#'
#' testthat::expect_equal(tmp2(cmp1), 0.0576417503974498)
#' testthat::expect_equal(tmp2(cmp2), 1.04950621760385)
#' testthat::expect_equal(tmp2(cmp3), 0.212977231452674)
#'
wasserstein_calculator = function(obs, debias = FALSE) {
  # TODO: bootstrapping

  # Data that can be precalculated:
  obs = sort(obs)
  av = mean(obs)
  # What is the average distance of points in the reference data from the mean?
  # we use this to normalise the result
  # these are inlined into the crated function
  norm = mean(abs(obs - av))
  size = length(obs)

  return(
    carrier::crate(
      function(sim, ...) {
        sim = .match_size(sim, !!size)
        sim = sort(sim)
        if (!!debias) {
          sim = sim - mean(sim) + !!av
        }
        return(mean(abs(sim - obs)) / !!norm)
      },
      obs = obs,
      .match_size = .match_size
    )
  )
}


# This recycles the input to match the size of comparison data by repeating it.
.match_size = carrier::crate(function(x, size) {
  if (length(x) < size) {
    x = c(
      rep(x, size %/% length(x)),
      sample(x, size %% length(x), replace = FALSE)
    )
  } else if (length(x) > size) {
    x = sample(x, size, replace = FALSE)
  }
  return(x)
})
