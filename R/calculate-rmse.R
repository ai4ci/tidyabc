#' Generate a function to calculate a Root Mean Squared Error (RMSE)
#'
#' This function takes reference data in the form, for example of count data, and
#' returns a crated function to calculate the mean squared error from
#' simulated data to the observed data.
#' \deqn{
#' \text{RMSE} = \sqrt{\frac{1}{N} \sum_{i=1}^{N} (x_i - y_i)^2}
#' }
#' where \eqn{x_i} are the simulated data points, \eqn{y_i} are the observed
#' data points, and \eqn{N} is the number of data points. Both input vectors
#' must be of the same length. Missing values (\code{NA}) are removed
#' pairwise before calculation. Returns \code{NA} if the input vectors are
#' not of equal length.
#'
#' @param sim A vector of simulated counts
#' @param obs A vector of observed counts
#'
#' @returns The square root of the average squared
#'   distance between simulation and observation. Simulation and observation
#'   must be the same length.
#' @export
#' @concept workflow
#'
#' @unit
#'
#' # zero if no distance
#' testthat::expect_equal(calculate_rmse(0:10, 0:10), 0)
#'
#' testthat::expect_equal(
#'   calculate_rmse(rep(5, 11), 0:10),
#'   3.16227766016838
#' )
#'
#' withr::with_seed(100, {
#'    ref = rnorm(1000)
#'    cmp1 = rnorm(1000)
#'    cmp2 = rnorm(1000, sd=2)
#' })
#'
#' breaks = seq(-2,2,length.out=8)
#' nref = table(cut(ref,breaks))
#' ncmp1 = table(cut(cmp1,breaks))
#' ncmp2 = table(cut(cmp2,breaks))
#'
#' tmp1 = calculate_rmse(ncmp1, nref)
#' tmp2 = calculate_rmse(ncmp2, nref)
#'
#' testthat::expect_equal(tmp1, 10.3578817470424)
#' testthat::expect_equal(tmp2, 68.8403951179829)
#'
#' @examples
#'
#' # example case counts from an exponential growth process
#' sim = table(floor(rexpgrowth(1000, 0.05, 40, 0)))
#' obs = table(floor(rexpgrowth(1000, 0.075, 40, 0)))
#' obs2 = table(floor(rexpgrowth(1000, 0.05, 40, 0)))
#'
#' # obs is a different distribution to sim (larger growth)
#' calculate_rmse(sim, obs)
#'
#' # obs2 is from the same distribution as sim so the RMSE should be lower:
#' calculate_rmse(sim, obs2)
#'
#'
calculate_rmse = function(sim, obs) {
  if (length(sim) != length(obs)) {
    return(NA)
  }

  invalid = is.na(sim) | is.na(obs)
  sim = sim[!invalid]
  obs = obs[!invalid]

  return(sqrt(mean((sim - obs)^2)))
}
