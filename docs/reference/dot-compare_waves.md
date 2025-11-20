# Generate comparison metrics for two sequential waves

Generate comparison metrics for two sequential waves

## Usage

``` r
.compare_waves(sim_df, prev_sim_df = NULL, priors_list, wave)
```

## Arguments

- sim_df:

  the output of a wave of simulation including a `abc_summary_distance`
  column

- prev_sim_df:

  the output of a previous ABC wave including a `abc_weight` column

- priors_list:

  the list of priors as a named list of `dist_fn`s, plus one sided
  formulae as constraints and 2 sided formulae as derived values. A
  correlation matrix may also be present as an attribute `cor`.

## Value

a nested tibble with 2 columns `summary` and `per_parameter` with stats
in each. The `summary` stats are
