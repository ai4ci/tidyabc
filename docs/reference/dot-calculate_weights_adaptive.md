# Calculate weights for particles in a new wave

The ABC weights need to be calculated for sampling from the proposal
distribution. They depend on the priors, an acceptance tolerance and the
proposal probability distribution.

## Usage

``` r
.calculate_weights_adaptive(
  sim_df,
  priors_list,
  acceptance_rate,
  proposal_list,
  kernel,
  use_proposal_correlation,
  ess_limit,
  max_recover
)
```

## Arguments

- sim_df:

  the output of a wave of simulation including a `abc_summary_distance`
  column

- priors_list:

  the list of priors as a named list of `dist_fn`s, plus one sided
  formulae as constraints and 2 sided formulae as derived values. A
  correlation matrix may also be present as an attribute `cor`.

- acceptance_rate:

  What proportion of simulations to keep in ABC rejection or hard ABC
  parts of the algorithms.

- proposal_list:

  a list of empirical probability distributions that map MVN space to
  proposal space, and are the "prior" for each adaptive wave. This is
  already used to generate the proposals and their mapping in `sim_df`

- kernel:

  one of `"epanechnikov"` (default), `"uniform"`, `"triangular"`,
  `"biweight"`, or `"gaussian"`. The kernel defines how the distance
  metric translates into the importance weight that decides whether a
  given simulation and associated parameters should be rejected or held
  for the next round. All kernels except `gaussian` have a hard cut-off
  outside of which the probability of acceptance of a particle is zero.
  Use of `gaussian` kernels can result in poor convergence.

- use_proposal_correlation:

  When calculating the weight of a particle the proposal correlation
  structure is available, to help determine how unusual or otherwise a
  particle is.

- ess_limit:

  a numeric vector of length 2 which for ABC adaptive, defines the
  limits which rate at which the algorithm will converge in terms of
  effective sample size. If for example the algorithm is converging too
  quickly and some high weight particles are dominating then the ESS
  will drop below the lower limit. In this case more particles will be
  accepted to try and offset this. On the other hand if the algorithm is
  converging too slowly low probability particles in proposal space are
  not filtered out quickly enough and this can lead to too much
  importance being given to unlikely proposals and wide bi-modal peaked
  posteriors.

- max_recover:

  if the effective sample size of SMC or adaptive algorithms drops below
  200, the algorithm will retry the wave with double the sample size to
  try and recover the shape of the distribution, up to a maximum of
  `max_recover` times.

## Value

the `sim_df` with an `abc_weight` column
