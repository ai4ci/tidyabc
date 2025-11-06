#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom lifecycle deprecated
## usethis namespace: end
NULL


#' Common workflow options
#'
#' @param obsdata The observational data. The data in this will typically
#'   be a named list, but could be anything, e.g. dataframe. It is the reference
#'   data that the simulation model is aiming to replicate.
#' @param priors_list a named list of priors specified as `dist_fns` (see
#'   `as.dist_fns()`), plus optionally, derived values as unnamed 2-sided
#'   formulae, where the LHS of the formula will be assigned to the value of the
#'   RHS, plus optionally a set of constraints as one sided formulae where the
#'   RHS of the formulae will resolve to a boolean value.
#' @param sim_fn a user defined function that takes a set of parameters named
#'   the same as the list `priors`. It must return a simulated data set in the
#'   same format as `obsdata`, or that can be compared to `simdata` by
#'   `scorer_fn`. This function must not refer to global parameters, and will be
#'   automatically crated with `carrier`.
#' @param simdata The simulated data. The simulated data must be in exactly
#'   the same format as `obsdata`, or that can be compared to `simdata` by
#'   `scorer_fn`.
#' @param scorer_fn a user supplied function that matches the following
#'   signature `scorer_fn(simdata, obsdata, ....)`, i.e. it takes data in the
#'   format of `simdata` paired with the original `obsdata` and returns a named
#'   list of scores per simulation. This function can make use of the
#'   `calculate_*()` set of functions to compare components of the simulation to
#'   the original data. This function must not refer to global parameters, and
#'   will be automatically crated with `carrier`. If this is a purrr style
#'   function then `.x` will refer to simulation output and `.y` to original
#'   observation data.
#' @param obsscores Summary scores for the observational data. This will
#'   be a named list, and is equivalent to the output of `scorer_fn`,
#'   on the observed data. If not given typically it will be assumed to be all
#'   zeros.
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
#' @param kernel_t A kernel bandwidth parameter for proposals. This controls the
#'   amount of noise that particles are perturbed by (and hence the spread of
#'   the PDF of the proposal distribution), and 1 is approximately 34% of the
#'   proposal distribution at the centre. Smaller values (default is 0.2) give a
#'   smaller step size in generating new proposals, and proposals will be closer
#'   to currently accepted particles.
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
#' @param allow_continue if SMC or adaptive algorithms have not converged after
#' `max_waves` allow the algorithm to interactively prompt the user to continue.
#' @name tidyabc_common
#' @keywords internal
#' @concept workflow
NULL
