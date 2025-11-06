#' Generate a function to calculate a root mean squared error
#'
#' This function takes reference data in the form, for example of count data, and
#' returns a crated function to calculate the mean squared error from
#' simulated data to the observed data.
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
calculate_rmse = function(sim, obs) {
  if (length(sim) != length(obs)) {
    return(NA)
  }

  return(sqrt(mean((sim - obs)^2)))
}
