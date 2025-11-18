#' Calculate a wasserstein distance
#'
#' This function takes simulation and observed data and calculates a normalised
#' wasserstein distance.
#'
#' In the comparison unequal lengths of the data can be accommodated. The
#' simulated data is sorted and linearly interpolated to the same length as the
#' observed data before the comparison.
#'
#' @param sim A vector of simulated data points
#' @param obs A vector of observed data points
#' @param debias Should the simulations be shifted to match the mean of the
#'   observed data
#' @param bootstraps Randomly resample from the simulated data points to match
#'   the observed size this many times and combine the output by averaging. The
#'   alternative, when this is 1 (the default) matches the sizes by selecting
#'   and/or repeating the simulated data points in order (deterministically)
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
#' calculate_wasserstein(cmp1,ref,bootstraps=10)
#'
#' gen = function(n, mean=0, sd=1) {
#'   sample(pnorm(seq(-1,1,length.out = n - n%%2 + 1), mean, sd))
#' }
#'
#' # there should be approximately zero
#' calculate_wasserstein(gen(1000), gen(1000))
#' calculate_wasserstein(gen(1000), gen(100))
#' calculate_wasserstein(gen(100), gen(1000))
#' calculate_wasserstein(gen(200), gen(1000))
#'
#' # these should be approximately equal:
#' calculate_wasserstein(gen(100,0.1), gen(1000))
#' calculate_wasserstein(gen(200,0.1), gen(1000))
#' calculate_wasserstein(gen(1000,0.1), gen(1000))
#' calculate_wasserstein(gen(1000, 0.1), gen(200))
#' calculate_wasserstein(gen(1000, 0.1), gen(100))
#'
calculate_wasserstein = function(sim, obs, debias = FALSE, bootstraps = 1) {
  sim = sim[!is.na(sim)]
  obs = obs[!is.na(obs)]

  obs = sort(obs)
  av = mean(obs)
  # What is the average distance of points in the reference data from the mean?
  # we use this to normalise the result
  norm = stats::sd(obs) # mean(abs(obs - av))
  size = length(obs)

  # no non zero items
  if (length(sim) == 0) {
    return(Inf)
  }

  out = sapply(seq_len(bootstraps), function(i) {
    if (bootstraps == 1) {
      sim2 = .match_size(sim, size)
    } else {
      sim2 = sample(sim, size, replace = TRUE)
      sim2 = sort(sim2)
    }

    if (debias) {
      sim2 = sim2 - mean(sim2) + av
    }
    return(mean(abs(sim2 - obs)) / norm)
  })

  return(mean(out))
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
#' @param bootstraps Randomly resample from the simulated data points to match
#'   the observed size this many times and combine the output by averaging. The
#'   alternative, when this is 1 (the default) matches the sizes by selecting
#'   and/or repeating the simulated data points in order (deterministically)
#'
#' @returns a function that takes parameter `sim` and returns a
#'   length normalised wasserstein distance. This is the average distance an
#'   individual data point must be shifted to match the reference data normalised
#'   by the average distance of the reference data from the mean. The function
#'   also takes a `p` parameter which can be a `progressr` progress bar which must
#'   be named.
#' @export
#' @concept workflow
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
wasserstein_calculator = function(obs, debias = FALSE, bootstraps = 1) {
  obs = obs[!is.na(obs)]

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
        sim = sim[!is.na(sim)]
        if (length(sim) == 0) {
          return(Inf)
        }

        out = sapply(seq_len(!!bootstraps), function(i) {
          if (!!bootstraps == 1) {
            sim2 = .match_size(sim, !!size)
          } else {
            sim2 = sample(sim, !!size, replace = TRUE)
            sim2 = sort(sim2)
          }

          if (!!debias) {
            sim2 = sim2 - mean(sim2) + !!av
          }
          return(mean(abs(sim2 - obs)) / !!norm)
        })

        return(mean(out))
      },
      obs = obs,
      .match_size = .match_size
    )
  )
}


# Deterministically match size of sorted vector.
# linear interpolation between values
.match_size = carrier::crate(function(x, size) {
  if (length(x) == 0) {
    stop("Cannot recycle object of length zero.")
  }

  idx = seq(1, length(x), length.out = size)
  p = idx - floor(idx)
  idx0 = floor(idx)
  idx1 = ceiling(idx)
  x = sort(x)

  out = ifelse(
    idx0 == idx1,
    x[idx0],
    (1 - p) * x[idx0] + p * x[idx1]
  )

  return(out)
})
