# Calculate weights for particles in a new wave

The ABC weights need to be calculated for sampling from the proposal
distribution. They depend on the priors, an acceptance tolerance and the
proposal probability distribution.

## Usage

``` r
.calculate_weights_adaptive(
  sim_df,
  priors_list,
  epsilon,
  proposal_list,
  kernel
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

- epsilon:

  epsilon is a tolerance threshold that controls how closely simulated
  summaries must match the observed ones to be considered plausible.
  This is in the unit of `abc_summary_distance`. Initially the 0.5
  quantile of distances, in subsequent waves this might be decreased. It
  is the scale parameter of the kernel function. \$K_h(\|u\|)\$

- proposal_list:

  a list of empirical probability distributions that map MVN space to
  proposal space, and are the "prior" for each adaptive wave. This is
  already used to generate the proposals and their mapping in `sim_df`

- kernel:

  one of "epanechnikov", "uniform", "triangular", "biweight", "gaussian"

## Value

the `sim_df` with an `abc_weight` column
