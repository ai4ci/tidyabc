# Spaghetti plot of resampled posterior fits

Spaghetti plot of resampled posterior fits

## Usage

``` r
plot_simulations(obsdata, fit, sim_fn, method = "auto", n = 200, ...)
```

## Arguments

- obsdata:

  The observational data. The data in this will typically be a named
  list, but could be anything, e.g. dataframe. It is the reference data
  that the simulation model is aiming to replicate.

- fit:

  A S3 `abc_fit` object as output by the `abc_XXX` functions

- sim_fn:

  a user defined function that takes a set of parameters named the same
  as the list `priors`. It must return a simulated data set in the same
  format as `obsdata`, or that can be compared to `simdata` by
  `scorer_fn`. This function must not refer to global parameters, and
  will be automatically crated with `carrier`.

- method:

  one of `"auto"`, `"count"`, or `"density"`, or a named vector one for
  each component of `obsdata`.

- n:

  the number of simulations to plot

- ...:

  Named arguments passed on to
  [`posterior_resample`](https://ai4ci.github.io/tidyabc/reference/posterior_resample.md)

  `n_resamples`

  :   the number of resamples for each parameter combination.

  `max_samples`

  :   the maximum total number of resamples to pick.

## Value

a patchwork of `ggplot`s

## Examples

``` r
plot_simulations(
  example_obsdata(),
  example_adaptive_fit(),
  example_sim_fn
)

```
