# Perform ABC sequential Monte Carlo fitting

This function will execute a simulation for a random selection of
parameters. Based on the `acceptance_rate` it will reject a proportion
of the results. The remaining results are weighted (using a kernel with
a tolerance equivalent to half the acceptance rate). Weighted parameter
particles generate proposals for further waves but a particle
perturbation. Waves are executed until a maximum is reached or the
results converge sufficiently that the changes between waves are small.
A relatively small number of simulations may be attempted with a high
acceptance rate, over multiple waves.

## Usage

``` r
abc_smc(
  obsdata,
  priors_list,
  sim_fn,
  scorer_fn,
  n_sims,
  acceptance_rate,
  ...,
  max_time = 5 * 60,
  converged_fn = default_termination_fn(),
  obsscores = NULL,
  distance_method = "euclidean",
  seed = NULL,
  parallel = FALSE,
  allow_continue = interactive(),
  debug_errors = FALSE,
  kernel = "epanechnikov",
  scoreweights = NULL
)
```

## Arguments

- obsdata:

  The observational data. The data in this will typically be a named
  list, but could be anything, e.g. dataframe. It is the reference data
  that the simulation model is aiming to replicate.

- priors_list:

  a named list of priors specified as a `abc_prior` S3 object (see
  [`priors()`](https://ai4ci.github.io/tidyabc/reference/priors.md)),
  this can include derived values as unnamed 2-sided formulae, where the
  LHS of the formula will be assigned to the value of the RHS, plus
  optionally a set of constraints as one sided formulae where the RHS of
  the formulae will resolve to a boolean value.

- sim_fn:

  a user defined function that takes a set of parameters named the same
  as the list `priors`. It must return a simulated data set in the same
  format as `obsdata`, or that can be compared to `simdata` by
  `scorer_fn`. This function must not refer to global parameters, and
  will be automatically crated with `carrier`.

- scorer_fn:

  a user supplied function that matches the following signature
  `scorer_fn(simdata, obsdata, ....)`, i.e. it takes data in the format
  of `simdata` paired with the original `obsdata` and returns a named
  list of component scores per simulation. This function can make use of
  the `calculate_*()` set of functions to compare components of the
  simulation to the original data. This function must not refer to
  global parameters, and will be automatically crated with `carrier`. If
  this is a purrr style function then `.x` will refer to simulation
  output and `.y` to original observation data.

- n_sims:

  The number of simulations to run per wave (for SMC and Adaptive) or
  overall (for Rejection). For rejection sampling a large number is
  recommended, for the others sma

- acceptance_rate:

  What proportion of simulations to keep in ABC rejection or hard ABC
  parts of the algorithms.

- ...:

  must be empty

- max_time:

  the maximum time in seconds to spend in ABC waves before admitting
  defeat. This time may not be all used if the algorithm converges.

- converged_fn:

  a function that takes a `summary` and `per_param` input and generates
  a logical indicator that the function has converged

- obsscores:

  Summary scores for the observational data. This will be a named list,
  and is equivalent to the output of `scorer_fn`, on the observed data.
  If not given typically it will be assumed to be all zeros.

- distance_method:

  what metric is used to combine `simscores` and `obsscores` and is one
  of `"euclidean"`, `"manhattan"`, or `"mahalanobis"`.

- seed:

  an optional random seed

- parallel:

  parallelise the simulation? If this is set to true then the simulation
  step will be parallelised using `furrr`. For this to make any
  difference it must have been set up with the following:
  `future::plan(future::multisession, workers = parallel::detectCores()-2)`

- allow_continue:

  if SMC or adaptive algorithms have not converged after `max_time`
  allow the algorithm to interactively prompt the user to continue.

- debug_errors:

  Errors that crop up in `sim_fn` during a simulation due to anomolous
  value combinations are hard to debug. If this flag is set, whenever a
  `sim_fn` or `scorer_fn` throws an error an interactive debugging
  session is started with the failing parameter combinations. This is
  not compatible with running in parallel.

- scoreweights:

  A named vector with names matching output of `scorer_fn` that defines
  the importance of this component of the scoring in the overall
  distance and weighting of any given simulation. This can be used to
  assign more weight on certain parts of the model output.

## Value

an S3 object of class `abc_fit` this contains the following:

- type: the type of ABC algorithm

- iterations: number of completed iterations

- converged: boolean - did the result meet convergence criteria

- waves: a list of dataframes of wave convergence metrics

- summary: a dataframe with the summary of the parameter fits after each
  wave.

- priors: the priors for the fit as a `abc_prior` S3 object

- posteriors: the final wave posteriors

## Examples

``` r
fit = abc_smc(
  obsdata = example_obsdata(),
  priors_list = example_priors_list(),
  sim_fn = example_sim_fn,
  scorer_fn = example_scorer_fn,
  n_sims = 1000,
  acceptance_rate = 0.25,
  max_time = 5, # 5 seconds to fit within examples limit
  parallel = FALSE,
  allow_continue = FALSE
)
#> ABC-SMC
#> SMC waves:  ■■■■■■■■■■                        30% | wave 2 ETA:  4s
#> SMC waves:  ■■■■■■■■■■■■■                     42% | wave 3 ETA:  3s
#> Converged on wave: 8
#> SMC waves:  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      90% | wave 7 ETA:  1s

summary(fit)
#> ABC SMC fit: 8 waves - (converged)
#> Parameter estimates:
#> # A tibble: 3 × 4
#> # Groups:   param [3]
#>   param mean_sd       median_95_CrI           ESS
#>   <chr> <chr>         <chr>                 <dbl>
#> 1 mean  4.981 ± 0.031 4.978 [4.911 — 5.074]  259.
#> 2 sd1   1.992 ± 0.070 1.991 [1.830 — 2.190]  259.
#> 3 sd2   0.988 ± 0.038 0.985 [0.895 — 1.090]  259.
```
