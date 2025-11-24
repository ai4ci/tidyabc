#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
#' @importFrom splines interpSpline
#' @importFrom rlang `:=`
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
#'   `priors()`).
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
#' @param distance_method what metric is used to combine `simscores` and `obsscores`.
#'   One of `"euclidean"`, `"normalised"`, `"manhattan"`, or `"mahalanobis"`.
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
#'
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
#'   more weight on certain parts of the model output. For `euclidean` and `manhattan`
#'   distance methods these weights multiply the output of `scorer_fn` directly.
#'   For the other 2 distance methods some degree of normalisation is done first
#'   on the first wave scores to make different components have approximately the
#'   same relevance to the overall score.
#' @param kernel one of `"epanechnikov"` (default), `"uniform"`, `"triangular"`,
#'   `"biweight"`, or `"gaussian"`. The kernel defines how the distance metric
#'   translates into the importance weight that decides whether a given
#'   simulation and associated parameters should be rejected or held for the
#'   next round. All kernels except `gaussian` have a hard cut-off outside of
#'   which the probability of acceptance of a particle is zero. Use of
#'   `gaussian` kernels can result in poor convergence.
#' @param distfit one of `"empirical"` or `"analytical"` determines what kind of
#'   distribution the ABC adaptive algorithm will fit for the posteriors.
#' @param bw for Adaptive ABC data distributions are smoothed before modelling
#'   empirical CDF. Over smoothing can reduce convergence rate, under-smoothing
#'   may result in noisy posterior estimates, and appearance of local modes.
#' @param widen_by change the dispersion of the empirical proposal distribution in ABC
#'   adaptive, preserving the median. This is akin to a nonlinear,
#'   heteroscedastic random walk in the quantile space, and can help address
#'   over-fitting or local modes in the ABC adaptive waves. `widen_by` is an odds
#'   ratio and describes how much further from the median any given part of the
#'   distribution is after transformation. E.g. if the median of a distribution
#'   is zero, and the `widen_by` is 2 then the 0.75 quantile will move to the
#'   position of the 0.9 quantile. The distribution will stay within the
#'   support of the prior. This is by default 1.05 which allows for some
#'   additional variability in proposals.
#' @param use_proposal_correlation When calculating the weight of a particle the
#'   proposal correlation structure is available, to help determine how unusual
#'   or otherwise a particle is.
#' @param ess_limit a numeric vector of length 2 which for ABC adaptive, defines
#'   the limits which rate at which the algorithm will converge in terms of
#'   effective sample size. If for example the algorithm is converging too quickly
#'   and some high weight particles are dominating then the ESS will drop
#'   below the lower limit. In this case more particles will be accepted to try
#'   and offset this. On the other hand if the algorithm is converging too slowly
#'   low probability particles in proposal space are not filtered out quickly
#'   enough and this can lead to too much importance being given to unlikely
#'   proposals and wide bi-modal peaked posteriors.
#'
#' @name tidyabc_common
#' @keywords internal
#' @concept workflow
NULL
