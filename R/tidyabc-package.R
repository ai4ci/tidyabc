#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
#' @importFrom splines interpSpline
NULL


#' Common workflow options
#'
#' @param obsdata The observational data. The data in this will typically
#'   be a named list, but could be anything, e.g. dataframe. It is the reference
#'   data that the simulation model is aiming to replicate.
#' @param priors_list a named list of priors specified as a `abc_prior` S3 object
#'   (see `priors()`), this can include derived values as unnamed 2-sided
#'   formulae, where the LHS of the formula will be assigned to the value of the
#'   RHS, plus optionally a set of constraints as one sided formulae where the
#'   RHS of the formulae will resolve to a boolean value.
#' @param proposal_list a named list of statistical or empirical distributions
#'   resulting from previous waves specified as a `abc_prior` S3 object (see
#'   \code{\link{priors()}}).
#' @param sim_fn a user defined function that takes a set of parameters named
#'   the same as `priors_list`. It must return a simulated data set in the
#'   same format as `obsdata`, or that can be compared to `simdata` by
#'   `scorer_fn`. This function must not refer to global parameters, and will be
#'   automatically crated with `carrier`.
#' @param simdata The simulated data. The simulated data must be in exactly
#'   the same format as `obsdata`, or that can be compared to `simdata` by
#'   `scorer_fn`.
#' @param scorer_fn a user supplied function that matches the following
#'   signature `scorer_fn(simdata, obsdata, ....)`, i.e. it takes data in the
#'   format of `simdata` paired with the original `obsdata` and returns a named
#'   list of component scores per simulation. This function can make use of the
#'   `calculate_*()` set of functions to compare components of the simulation to
#'   the original data. This function must not refer to global parameters, and
#'   will be automatically crated with `carrier`. If this is a purrr style
#'   function then `.x` will refer to simulation output and `.y` to original
#'   observation data.
#' @param obsscores Summary scores for the observational data. This will
#'   be a named list, and is equivalent to the output of `scorer_fn`,
#'   on the observed data. If not given typically it will be assumed to be all
#'   zeros.
#' @param n_sims The number of simulations to run per wave (for SMC and Adaptive)
#'   or overall (for Rejection). For rejection sampling a large number is
#'   recommended, for the others sma
#' @param simscores Scores for the simulated data, as output by `scorer_fn`.
#' @param distance_method what metric is used to combine `simscores` and `obsscores`
#'   and is one of `"euclidean"`, `"manhattan"`, or `"mahalanobis"`.
#' @param keep_simulations keep the individual simulation results in the output
#'   of an ABC workflow. This can have large implications for the size of the
#'   result. It may also not be what you want and it is probably worth considering
#'   resampling the posteriors rather than keeping the simulations.
#' @param acceptance_rate What proportion of simulations to keep in ABC rejection
#'   or hard ABC parts of the algorithms.
#' @param posteriors_df a dataframe of posteriors that have been selected by ABC
#'   this may include columns for scores, weight and/or simulation outputs (
#'   `abc_component_score`, `abc_summary_distance`, `abc_weight`, `abc_simulation`
#'   ) as well as columns matching the `priors` input specification.
#' @param converged_fn a function that takes a `summary` and `per_param` input
#'   and generates a logical indicator that the function has converged
#' @param parallel parallelise the simulation? If this is set to true then the
#'   simulation step will be parallelised using `furrr`. For this to make any
#'   difference it must have been set up with the following:
#'   `future::plan(future::multisession, workers = parallel::detectCores()-2)`
#' @param fit A S3 `abc_fit` object as output by the `abc_XXX` functions
#' @param p a `progressr` progress bar
#' @param seed an optional random seed
#' @param knots the number of knots to model the CDF with. Optional, and will be
#'   typically inferred from the data size. Small numbers tend to work better if
#'   we expect the distribution to be unimodal.
#' @param bw for Adaptive ABC data distributions are smoothed before modelling
#'   the CDF. Over smoothing can reduce convergence, under-smoothing may result
#'   in noisy posterior estimates. This is in units of the ESS and defaults to
#'   0.1.
#' @param max_time the maximum time in seconds to spend in ABC waves before admitting
#'   defeat. This time may not be all used if the algorithm converges.
#' @param allow_continue if SMC or adaptive algorithms have not converged after
#' `max_time` allow the algorithm to interactively prompt the user to continue.
#' @param max_recover if the effective sample size of SMC or adaptive algorithms
#'  drops below 200, the algorithm will retry the wave with double the sample
#'  size to try and recover the shape of the distribution, up to a maximum of
#' `max_recover` times.
#' @param debug_errors Errors that crop up in `sim_fn` during a simulation due
#'   to anomolous value combinations are hard to debug. If this flag is set,
#'   whenever a `sim_fn` or `scorer_fn` throws an error an interactive debugging
#'   session is started with the failing parameter combinations. This is not
#'   compatible with running in parallel.
#' @param scoreweights A named vector with names matching output of `scorer_fn`
#'   that defines the importance of this component of the scoring in the overall
#'   distance and weighting of any given simulation. This can be used to assign
#'   more weight on certain parts of the model output.
#' @param kernel one of `"epanechnikov"`, `"uniform"`, `"triangular"`, `"biweight"`,
#'   or `"gaussian"`. The kernel defines how the distance metric translates into
#'   the importance weight that decides whether a given simulation and associated
#'   parameters should be rejected or held for the next round.
#' @name tidyabc_common
#' @keywords internal
#' @concept workflow
NULL
