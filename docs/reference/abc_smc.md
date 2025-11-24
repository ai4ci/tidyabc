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
  as `priors_list`. It must return a simulated data set in the same
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

  what metric is used to combine `simscores` and `obsscores`. One of
  `"euclidean"`, `"normalised"`, `"manhattan"`, or `"mahalanobis"`.

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

- kernel:

  one of `"epanechnikov"` (default), `"uniform"`, `"triangular"`,
  `"biweight"`, or `"gaussian"`. The kernel defines how the distance
  metric translates into the importance weight that decides whether a
  given simulation and associated parameters should be rejected or held
  for the next round. All kernels except `gaussian` have a hard cut-off
  outside of which the probability of acceptance of a particle is zero.
  Use of `gaussian` kernels can result in poor convergence.

- scoreweights:

  A named vector with names matching output of `scorer_fn` that defines
  the importance of this component of the scoring in the overall
  distance and weighting of any given simulation. This can be used to
  assign more weight on certain parts of the model output. For
  `euclidean` and `manhattan` distance methods these weights multiply
  the output of `scorer_fn` directly. For the other 2 distance methods
  some degree of normalisation is done first on the first wave scores to
  make different components have approximately the same relevance to the
  overall score.

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

## Details

Performs the ABC Sequential Monte Carlo (SMC) algorithm. This iterative
method refines parameter estimates across multiple waves.

1.  **Initialization (Wave 1):** Parameters \\\theta^{(i)}\\ are sampled
    from the prior \\P(\theta)\\. Simulations \\D_s^{(i)} =
    M(\theta^{(i)})\\ are run, summaries \\S_s^{(i)}\\ are computed, and
    distances \\d^{(i)} = d(S_s^{(i)}, S_o)\\ are calculated. A
    tolerance threshold \\\epsilon_1\\ is set as the \\\alpha =
    \texttt{acceptance\\rate}\\ quantile of these initial distances.
    Unnormalized weights \\\tilde{w}^{(i)}\_1\\ are calculated using a
    kernel \\K\_{\epsilon_1}(d^{(i)})\\.

2.  **Subsequent Waves (\\t \> 1\\):**

    - **Proposal Generation:** A particle \\\theta^{(j)}\_{t-1}\\ is
      selected from the previous wave's accepted particles with
      probability proportional to its weight \\w^{(j)}\_{t-1}\\. The
      particle is then perturbed in a transformed MVN space using a
      multivariate normal kernel with covariance \\\Sigma_t =
      \frac{\kappa_t^2}{d} \text{Cov}\_{w\_{t-1}}(\theta\_{t-1})\\,
      where \\\text{Cov}\_{w\_{t-1}}\\ is the weighted covariance from
      wave \\t-1\\, \\d\\ is the parameter dimension, and \\\kappa_t\\
      is the `kernel_t` parameter. The new proposal \\\theta^{(i)}\_t\\
      is generated as: \$\$ \theta^{(i)}\_t = \theta^{(j)}\_{t-1} +
      \zeta, \quad \zeta \sim \mathcal{N}(0, \Sigma_t) \$\$ This
      proposal is mapped back to the original parameter space using the
      prior's copula transformation (from MVN space defined by prior
      CDFs).

    - **Simulation and Weighting:** Simulations \\D_s^{(i)} =
      M(\theta^{(i)}\_t)\\ are run for the new proposals. Distances
      \\d^{(i)}\_t\\ are computed. The tolerance \\\epsilon_t\\ is set
      as the \\\alpha\\-quantile of distances from the *current* wave's
      simulations. The unnormalized weight for particle \\i\\ in wave
      \\t\\ is calculated as: \$\$ \tilde{w}^{(i)}\_t =
      \frac{P(\theta^{(i)}\_t)
      K\_{\epsilon_t}(d^{(i)}\_t)}{q_t(\theta^{(i)}\_t)} \$\$ where
      \\P(\theta^{(i)}\_t)\\ is the prior density, \\K\_{\epsilon_t}\\
      is the ABC kernel, and \\q_t(\theta^{(i)}\_t)\\ is the proposal
      density from the previous wave's weighted particles (calculated
      using the perturbation kernel). This proposal density is computed
      as a weighted sum: \$\$ q_t(\theta^{(i)}\_t) = \sum_j
      w^{(j)}\_{t-1} \phi(\theta^{(i)}\_t; \theta^{(j)}\_{t-1},
      \Sigma_t) \$\$ where \\\phi(\cdot; \mu, \Sigma)\\ is the PDF of a
      multivariate normal with mean \\\mu\\ and covariance \\\Sigma\\.

3.  **Normalization:** Weights \\w^{(i)}\_t\\ are normalized to sum to
    one. Particles with negligible weights are typically filtered out.

4.  **Termination:** The process repeats until a maximum number of waves
    or time is reached, or convergence criteria are met based on changes
    in parameter estimates or effective sample size (ESS).

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
#> SMC waves:  ■■■■■■■                           21% | wave 1 ETA:  4s
#> SMC waves:  ■■■■■■■■■■■■■■■■■■■               60% | wave 4 ETA:  2s
#> Converged on wave: 7
#> SMC waves:  ■■■■■■■■■■■■■■■■■■■■■■■■■■■       88% | wave 6 ETA:  1s

summary(fit)
#> ABC SMC fit: 7 waves - (converged)
#> Parameter estimates:
#> # A tibble: 3 × 4
#> # Groups:   param [3]
#>   param mean_sd       median_95_CrI           ESS
#>   <chr> <chr>         <chr>                 <dbl>
#> 1 mean  4.978 ± 0.036 4.979 [4.883 — 5.097]  274.
#> 2 sd1   1.998 ± 0.097 1.996 [1.769 — 2.310]  274.
#> 3 sd2   0.990 ± 0.045 0.990 [0.879 — 1.102]  274.
```
