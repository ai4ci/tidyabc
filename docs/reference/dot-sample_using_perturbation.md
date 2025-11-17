# Generate a new set of particles from a previous wave using perturbation

This set of particles will be close to existing ones, depending on the
kernel_t parameter.

## Usage

``` r
.sample_using_perturbation(priors_list, n_sims, prev_sim_df)
```

## Arguments

- priors_list:

  the list of priors as a named list of `dist_fn`s, plus one sided
  formulae as constraints and 2 sided formulae as derived values. A
  correlation matrix may also be present as an attribute `cor`.

- n_sims:

  The number of simulations to run per wave (for SMC and Adaptive) or
  overall (for Rejection). For rejection sampling a large number is
  recommended, for the others sma

- prev_sim_df:

  the output of a previous ABC wave including a `abc_weight` column

## Value

a `sim_df`
