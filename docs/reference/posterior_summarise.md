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
#> $ mean      <dbl> 4.964376, 2.123894, 1.139277
#> $ sd        <dbl> 0.1292487, 0.2332786, 0.1733625
#> $ quantiles <list> [<tbl_df[9 x 2]>], [<tbl_df[9 x 2]>], [<tbl_df[9 x 2]>]
#> $ q.0.025   <dbl> 4.7174323, 1.5218539, 0.7531816
#> $ q.0.5     <dbl> 4.957561, 2.077933, 1.133596
#> $ q.0.975   <dbl> 5.260754, 2.904667, 1.531218
#> $ density   <distfn[]> [mean; Median (IQR) 4.96 [4.85 — 5.07]], [sd1; Median (IQR) …
#> $ ESS       <dbl> 1518.651, 1518.651, 1518.651
```
