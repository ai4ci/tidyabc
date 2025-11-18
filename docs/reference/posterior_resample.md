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
  as `priors_list`. It must return a simulated data set in the same
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
#>  [1] 6.219357 3.019306 5.014118 3.701336 3.322469 4.511440 6.315831 5.058636
#>  [9] 8.167976 4.697854
#> 
#> $B
#>  [1] 5.669615 6.052935 4.118884 4.131646 4.267825 5.405490 4.433625 5.404898
#>  [9] 4.138371 4.486834
#> 
```
