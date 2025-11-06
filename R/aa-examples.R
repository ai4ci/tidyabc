# example data cache
exdata = new.env(parent = emptyenv())

#' Example generators
#'
#' These are a set of internally cached functions to support
#' examples. They are exported as internal functions so that the examples can
#' run correctly and cache their output to prevent excessive repetition of the
#' code examples.
#'
#' @keywords internal
#' @name example_fns
NULL


#' @describeIn example_fns Example simulation function
#' @param mean the example mean
#' @param sd1 the example normal sd (A)
#' @param sd2 the example gamma sd (B)
#' @returns a list of A and B of samples from a normal and gamma
#' @export
#' @examples
#' example_sim_fn(3,2,1) %>% lapply(head,10)
example_sim_fn = function(mean, sd1, sd2) {
  return(list(
    A = rnorm(1000, mean, sd1),
    B = rgamma2(1000, mean, sd2)
  ))
}

#' @describeIn example_fns Example simulation function
#' @param obsdata e.g. example_obsdata()
#' @param simdata output of example_sim_fn()
#'
#' @returns a list of scores
#' @export
#' @examples
#' example_scorer_fn(
#'   example_obsdata(),
#'   example_sim_fn(3,2,1)
#' )
example_scorer_fn = function(obsdata, simdata) {
  return(list(
    A = calculate_wasserstein(obsdata$A, simdata$A),
    B = calculate_wasserstein(obsdata$B, simdata$B)
  ))
}


#' @describeIn example_fns Example output of `test_simulation()`
#' @returns output of `test_simulation()`
#' @export
#' @examples
#' example_obs() %>% rapply(head,n=10, how="replace")
example_obs = function() {
  if (is.null(exdata$obs)) {
    exdata$obs = test_simulation(
      example_sim_fn,
      example_scorer_fn,
      mean = 5,
      sd1 = 2,
      sd2 = 1,
      seed = 123
    )
  }
  return(exdata$obs)
}

#' @describeIn example_fns Example for `truth` parameter
#' @returns example `truth` parameter
#' @export
#' @examples
#' example_truth()
example_truth = function() {
  return(example_obs()$truth)
}

#' @describeIn example_fns Example for `obsdata` parameter
#' @returns example `obsdata` parameter
#' @export
#' @examples
#' example_obsdata() %>% lapply(head,10)
example_obsdata = function() {
  return(example_obs()$obsdata)
}

#' @describeIn example_fns Example for `priors_list` parameter
#' @returns example `priors_list` parameter
#' @export
#' @examples
#' example_priors_list()
example_priors_list = function() {
  return(
    list(
      mean = as.dist_fns(runif, 0, 10),
      sd1 = as.dist_fns(runif, 0, 5),
      sd2 = as.dist_fns(runif, 0, 5)
    )
  )
}

#' @describeIn example_fns Example `abc_rejection` output
#' @returns example `abc_rejection` output
#' @export
#' @examples
#' example_rejection_fit()
example_rejection_fit = function() {
  if (is.null(exdata$rejection_fit)) {
    exdata$rejection_fit = abc_rejection(
      obsdata = example_obsdata(),
      priors_list = example_priors_list(),
      sim_fn = example_sim_fn,
      scorer_fn = example_scorer_fn,
      n_sims = 10000,
      acceptance_rate = 0.01,
      parallel = FALSE
    )
  }
  return(exdata$rejection_fit)
}

#' @describeIn example_fns Example `abc_smc` output
#' @returns example `abc_smc` output
#' @export
#' @examples
#' example_smc_fit()
example_smc_fit = function() {
  if (is.null(exdata$smc_fit)) {
    exdata$smc_fit = abc_smc(
      obsdata = example_obsdata(),
      priors_list = example_priors_list(),
      sim_fn = example_sim_fn,
      scorer_fn = example_scorer_fn,
      n_sims = 1000,
      acceptance_rate = 0.25,
      max_waves = 5,
      parallel = FALSE,
      allow_continue = FALSE
    )
  }
  return(exdata$smc_fit)
}

#' @describeIn example_fns Example `abc_adaptive` output
#' @returns example `abc_adaptive` output
#' @export
#' @examples
#' example_adaptive_fit()
example_adaptive_fit = function() {
  if (is.null(exdata$adaptive_fit)) {
    exdata$adaptive_fit = abc_adaptive(
      obsdata = example_obsdata(),
      priors_list = example_priors_list(),
      sim_fn = example_sim_fn,
      scorer_fn = example_scorer_fn,
      n_sims = 1000,
      acceptance_rate = 0.25,
      max_waves = 5,
      parallel = FALSE,
      allow_continue = FALSE
    )
  }
  return(exdata$adaptive_fit)
}
