# Calculate a basket of summaries from a weighted list of posterior samples

Calculate a basket of summaries from a weighted list of posterior
samples

## Usage

``` r
posterior_summarise(
  posteriors_df,
  priors_list,
  p = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
)
```

## Arguments

- posteriors_df:

  a dataframe of posteriors that have been selected by ABC this may
  include columns for scores, weight and/or simulation outputs (
  `abc_component_score`, `abc_summary_distance`, `abc_weight`,
  `abc_simulation` ) as well as columns matching the `priors` input
  specification.

- priors_list:

  a named list of priors specified as a `abc_prior` S3 object (see
  [`priors()`](https://ai4ci.github.io/tidyabc/reference/priors.md)),
  this can include derived values as unnamed 2-sided formulae, where the
  LHS of the formula will be assigned to the value of the RHS, plus
  optionally a set of constraints as one sided formulae where the RHS of
  the formulae will resolve to a boolean value.

- p:

  a `progressr` progress bar

## Value

a dataframe indexed by parameter with useful summary metrics.

## Examples

``` r
fit = example_adaptive_fit()
summ = posterior_summarise(fit$posteriors, fit$priors)

summ %>% dplyr::glimpse()
#> Rows: 3
#> Columns: 9
#> Groups: param [3]
#> $ param     <chr> "mean", "sd1", "sd2"
#> $ mean      <dbl> 4.979629, 2.090683, 1.014513
#> $ sd        <dbl> 0.03592986, 0.06405630, 0.03351014
#> $ quantiles <list> [<tbl_df[9 x 2]>], [<tbl_df[9 x 2]>], [<tbl_df[9 x 2]>]
#> $ q.0.025   <dbl> 4.893898, 1.769493, 0.925236
#> $ q.0.5     <dbl> 4.972857, 2.094028, 1.011725
#> $ q.0.975   <dbl> 5.075677, 2.214940, 1.167002
#> $ density   <distfn[]> [mean; Median (IQR) 4.97 [4.95 — 5.01]], [sd1; Median (IQR) …
#> $ ESS       <dbl> 262.8049, 262.8049, 262.8049
```
