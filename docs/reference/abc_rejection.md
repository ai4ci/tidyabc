# Perfom simple ABC rejection algorithm

This function will execute a simulation for a random selection of
parameters and identify the best matching `acceptance_rate` percent, as
defined by the summary distance metric. A large number of simulations
and a low acceptance rate are best here.

## Usage

``` r
abc_rejection(
  obsdata,
  priors_list,
  sim_fn,
  scorer_fn,
  n_sims,
  acceptance_rate,
  ...,
  converged_fn = default_termination_fn(),
  obsscores = NULL,
  distance_method = "euclidean",
  keep_simulations = FALSE,
  seed = NULL,
  parallel = FALSE,
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

- keep_simulations:

  keep the individual simulation results in the output of an ABC
  workflow. This can have large implications for the size of the result.
  It may also not be what you want and it is probably worth considering
  resampling the posteriors rather than keeping the simulations.

- seed:

  an optional random seed

- parallel:

  parallelise the simulation? If this is set to true then the simulation
  step will be parallelised using `furrr`. For this to make any
  difference it must have been set up with the following:
  `future::plan(future::multisession, workers = parallel::detectCores()-2)`

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

Parameters \\\theta^{(i)}\\ are sampled independently from the prior
distribution \\P(\theta)\\ for \\i = 1, \dots, n\_{\text{sims}}\\. For
each \\\theta^{(i)}\\, simulated data \\D_s^{(i)} = M(\theta^{(i)})\\ is
generated via the simulator function `sim_fn`, and a vector of summary
statistics \\S_s^{(i)} = \texttt{scorer\\fn}(D_s^{(i)})\\ is computed
and compared to the observed summary statistics \\S_o =
\texttt{scorer\\fn}(D_o)\\.

A distance metric \\d^{(i)} = d(S_s^{(i)}, S_o)\\ is computed. By
default, this is the Euclidean distance: \$\$ d^{(i)} = \left\\ W \circ
(S_s^{(i)} - S_o) \right\\\_2, \$\$ where \\W\\ is a vector of optional
summary statistic weights (`scoreweights`), and \\\circ\\ denotes
element-wise multiplication. Other supported metrics include Manhattan
(\\\ell_1\\) and Mahalanobis distance (using the empirical covariance
from the first wave).

The tolerance threshold \\\epsilon\\ is set to the \\\alpha =
\texttt{acceptance\\rate}\\ quantile of the distances
\\\\d^{(i)}\\\_{i=1}^{n\_{\text{sims}}}\\: \$\$ \epsilon =
\text{quantile}\big(\\d^{(i)}\\, \alpha\big). \$\$

Unnormalized ABC weights are then assigned using a kernel function
\\K\_\epsilon(\cdot)\\: \$\$ \tilde{w}^{(i)} =
K\_\epsilon\big(d^{(i)}\big), \$\$ where \\K\_\epsilon(d)\\ is one of
the kernels defined in `kernels.R` (e.g., Epanechnikov: \\K\_\epsilon(d)
= \frac{3}{4\epsilon} (1 - d^2\\/\\\epsilon^2)\\ \mathbb{I}(d \leq
\epsilon)\\). These weights are then transformed via a logistic
("expit") function and normalized to sum to one: \$\$ w^{(i)} =
\frac{\text{expit}(\log \tilde{w}^{(i)})}{ \sum_j \text{expit}(\log
\tilde{w}^{(j)})} = \frac{ \tilde{w}^{(i)} / (1 + \tilde{w}^{(i)}) }{
\sum_j \tilde{w}^{(j)} / (1 + \tilde{w}^{(j)}) }. \$\$ The resulting
weighted sample \\\\(\theta^{(i)}, w^{(i)})\\\\ approximates the ABC
posterior distribution \\P\_\epsilon(\theta \mid D_o)\\.

## Examples

``` r
fit = abc_rejection(
  example_obsdata(),
  example_priors_list(),
  example_sim_fn,
  example_scorer_fn,
  n_sims = 10000,
  acceptance_rate = 0.01
)
#> ABC rejection, 1 wave.

summary(fit)
#> ABC rejection fit: single wave
#> Parameter estimates:
#> # A tibble: 3 × 4
#> # Groups:   param [3]
#>   param mean_sd       median_95_CrI           ESS
#>   <chr> <chr>         <chr>                 <dbl>
#> 1 mean  4.974 ± 0.197 5.001 [4.559 — 5.401]  69.2
#> 2 sd1   2.382 ± 0.559 2.267 [1.548 — 4.033]  69.2
#> 3 sd2   1.209 ± 0.294 1.148 [0.717 — 2.067]  69.2
```
