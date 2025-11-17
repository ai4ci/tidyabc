# Generate a set of samples from selected posteriors

Once an ABC model fitting is complete the simulation data is generally
only one possible realisation of the parameters, which has been selected
for closeness. To properly compare the output with the observed data we
need a set of posterior re-samples, which are selected from posteriors
according to importance.

## Usage

``` r
posterior_resample(
  posteriors_df,
  sim_fn,
  n_resamples = 1,
  seed = NULL,
  parallel = FALSE,
  max_samples = 200
)
```

## Arguments

- posteriors_df:

  a dataframe of posteriors that have been selected by ABC this may
  include columns for scores, weight and/or simulation outputs (
  `abc_component_score`, `abc_summary_distance`, `abc_weight`,
  `abc_simulation` ) as well as columns matching the `priors` input
  specification.

- sim_fn:

  a user defined function that takes a set of parameters named the same
  as the list `priors`. It must return a simulated data set in the same
  format as `obsdata`, or that can be compared to `simdata` by
  `scorer_fn`. This function must not refer to global parameters, and
  will be automatically crated with `carrier`.

- n_resamples:

  the number of resamples for each parameter combination.

- seed:

  an optional random seed

- parallel:

  parallelise the simulation? If this is set to true then the simulation
  step will be parallelised using `furrr`. For this to make any
  difference it must have been set up with the following:
  `future::plan(future::multisession, workers = parallel::detectCores()-2)`

- max_samples:

  the maximum total number of resamples to pick.

## Value

a dataframe of the posteriors with an `abc_sim` list column containing
the output of `sim_fn` called with the parameters in that row.

## Examples

``` r
fit = example_adaptive_fit()

sample_df = posterior_resample(
  fit$posteriors,
  sim_fn = example_sim_fn
)

# the fitted simulations are in the `abc_sim` column
sim1 = sample_df$abc_sim[[1]]

sim1 %>% lapply(head, 10)
#> $A
#>  [1] 8.487740 3.815985 3.193516 2.381563 2.653476 4.677163 5.866941 4.230691
#>  [9] 6.327211 1.928368
#> 
#> $B
#>  [1] 4.706729 4.890665 4.435608 5.681729 5.119136 4.995636 4.043788 4.154802
#>  [9] 5.947107 5.858067
#> 
```
