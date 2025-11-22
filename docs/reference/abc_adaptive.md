# Perform ABC sequential adaptive fitting

This function will execute a simulation for a random selection of
parameters. Based on the `acceptance_rate` it will reject a proportion
of the results. The remaining results are weighted (using a kernel with
a tolerance equivalent to half the acceptance rate). Empirical
distributions are fitted to weighted parameter particles and from these
proposals are generated for further waves by fresh sampling. Waves are
executed until a maximum is reached or the results converge sufficiently
that the changes between waves are small. A relatively small number of
simulations may be attempted with a high acceptance rate, over multiple
waves.

## Usage

``` r
abc_adaptive(
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
  knots = NULL,
  parallel = FALSE,
  max_recover = 3,
  allow_continue = interactive(),
  debug_errors = FALSE,
  kernel = "epanechnikov",
  bw = 0.1,
  widen_by = 1.05,
  scoreweights = NULL,
  use_proposal_correlation = TRUE
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

- knots:

  the number of knots to model the CDF with. Optional, and will be
  typically inferred from the data size. Small numbers tend to work
  better if we expect the distribution to be unimodal.

- parallel:

  parallelise the simulation? If this is set to true then the simulation
  step will be parallelised using `furrr`. For this to make any
  difference it must have been set up with the following:
  `future::plan(future::multisession, workers = parallel::detectCores()-2)`

- max_recover:

  if the effective sample size of SMC or adaptive algorithms drops below
  200, the algorithm will retry the wave with double the sample size to
  try and recover the shape of the distribution, up to a maximum of
  `max_recover` times.

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

- bw:

  for Adaptive ABC data distributions are smoothed before modelling the
  CDF. Over smoothing can reduce convergence rate, under-smoothing may
  result in noisy posterior estimates, and appearance of local modes.
  This is a proportion of the ESS and defaults to 0.1.

- widen_by:

  change the dispersion of proposal distribution in ABC adaptive,
  preserving the median. This is akin to a nonlinear, heteroscedastic
  random walk in the quantile space, and can help address over-fitting
  or local modes in the ABC adaptive waves. `widen_by` is an odds ratio
  and describes how much further from the median any given part of the
  distribution is after transformation. E.g. if the median of a
  distribution is zero, and the `widen_by` is 2 then the 0.75 quantile
  will move to the position of the 0.9 quantile. The distribution will
  stay within the support of the prior. This is by default 1.05 which
  allows for some additional variability in proposals.

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

- use_proposal_correlation:

  When calculating the weight of a particle the proposal correlation
  structure is available, to help determine how unusual or otherwise a
  particle is.

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

Performs the ABC Adaptive algorithm. This iterative method refines
parameter estimates across waves by fitting empirical proposal
distributions to the weighted posterior samples from the previous wave.
Unlike ABC-SMC, which uses a fixed perturbation kernel, `abc_adaptive`
constructs a new proposal distribution \\Q_t(\theta)\\ at each wave
\\t\\.

1.  **Initialization (Wave 1):** Parameters \\\theta^{(i)}\\ are sampled
    from the prior \\P(\theta)\\. Simulations are run, summary
    statistics \\S_s^{(i)}\\ are computed, and distances \\d^{(i)} =
    d(S_s^{(i)}, S_o)\\ are calculated. A tolerance threshold
    \\\epsilon_1\\ is set as the \\\alpha = \texttt{acceptance\\rate}\\
    quantile of these distances. Unnormalized weights
    \\\tilde{w}^{(i)}\_1\\ are calculated using a kernel
    \\K\_{\epsilon_1}(d^{(i)})\\.

2.  **Subsequent Waves (\\t \> 1\\):**

    - **Proposal Generation:** An empirical joint proposal distribution
      \\Q_t(\theta)\\ is constructed from the weighted posterior sample
      \\\\(\theta^{(i)}\_{t-1}, w^{(i)}\_{t-1})\\\\ of the previous
      wave. This is done by fitting marginal empirical distributions
      \\Q\_{t,j}(\theta_j)\\ to each parameter \\\theta_j\\, using the
      [`empirical()`](https://ai4ci.github.io/tidyabc/reference/empirical.md)
      function with the prior \\P_j(\theta_j)\\ as a link to enforce
      support constraints. These marginals are assumed independent, but
      their weighted correlation matrix \\R_t\\ is retained as an
      attribute and used to induce dependence in the MVN sampling space.
      New proposals \\\theta^{(i)}\_t\\ are generated by:

      - Sampling a vector \\Z \sim \mathcal{N}(0, R_t)\\ in a correlated
        standard normal space.

      - Mapping each component \\Z_j\\ to uniform space: \\U_j =
        \Phi(Z_j)\\.

      - Mapping to the parameter space using the empirical quantile
        functions: \\\theta^{(i)}\_{t,j} = Q\_{t,j}^{-1}(U_j)\\.

    - **Simulation and Weighting:** Simulations are run for the new
      proposals. Distances \\d^{(i)}\_t\\ are computed and the tolerance
      \\\epsilon_t\\ is set as the \\\alpha\\-quantile of the current
      wave's distances. The unnormalized weight for particle \\i\\ in
      wave \\t\\ is calculated as: \$\$ \tilde{w}^{(i)}\_t =
      \frac{P(\theta^{(i)}\_t)
      K\_{\epsilon_t}(d^{(i)}\_t)}{Q_t(\theta^{(i)}\_t)} \$\$ where
      \\P(\theta^{(i)}\_t) = \prod_j P_j(\theta^{(i)}\_{t,j})\\ is the
      prior density (assuming independence), \\K\_{\epsilon_t}\\ is the
      ABC kernel, and \\Q_t(\theta^{(i)}\_t) = \prod_j
      Q\_{t,j}(\theta^{(i)}\_{t,j})\\ is the empirical proposal density
      (also assuming independence for density calculation, consistent
      with the marginal fitting). The correlation structure is handled
      in the sampling process, and optionally in the density evaluation.

3.  **Normalization:** Weights \\w^{(i)}\_t\\ are normalized to sum to
    one. The algorithm includes a recovery mechanism: if the Effective
    Sample Size (ESS) falls below a threshold (e.g., 200), the
    acceptance rate is increased (i.e., \\\epsilon_t\\ is made larger)
    to accept more particles and improve the ESS.

4.  **Termination:** The process repeats until a maximum time is reached
    or convergence criteria based on parameter stability and credible
    interval contraction are met.

## Examples

``` r
fit = abc_adaptive(
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
#> ABC-Adaptive
#> Adaptive waves:  ■■■■■■■■                          24% | wave 2 ETA:  4s
#> Adaptive waves:  ■■■■■■■■■■■■■■■■■■■■■■■■■         81% | wave 7 ETA:  1s
#> Converged on wave: 8

summary(fit)
#> ABC adaptive fit: 8 waves - (converged)
#> Parameter estimates:
#> # A tibble: 3 × 4
#> # Groups:   param [3]
#>   param mean_sd       median_95_CrI           ESS
#>   <chr> <chr>         <chr>                 <dbl>
#> 1 mean  4.983 ± 0.033 4.986 [4.895 — 5.082]  354.
#> 2 sd1   1.984 ± 0.059 1.972 [1.834 — 2.173]  354.
#> 3 sd2   1.029 ± 0.046 1.033 [0.855 — 1.120]  354.
```
