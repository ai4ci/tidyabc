# Fit empirical distribution to posterior samples for generating more waves

This function allows "updating" of the prior with an empirical posterior
distribution which will retain the bounds of the prior, and is
effectively a spline based transform of the prior distribution, based on
the density of data in prior space. This gives a clean density when data
is close to a prior distribution limit and work better than a standard
density for

## Usage

``` r
posterior_fit_empirical(posteriors_df, priors_list, knots = NULL, bw = 0.1)
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

- knots:

  the number of knots to model the CDF with. Optional, and will be
  typically inferred from the data size. Small numbers tend to work
  better if we expect the distribution to be unimodal.

- bw:

  for Adaptive ABC data distributions are smoothed before modelling the
  CDF. Over smoothing can reduce convergence, under-smoothing may result
  in noisy posterior estimates. This is in units of the ESS and defaults
  to 0.1.

## Value

an `abc_prior` S3 object approximating the distribution of the posterior
samples, and adhering to the support of the provided priors.

## Examples

``` r
fit = example_smc_fit()
proposals = posterior_fit_empirical(fit$posteriors, fit$priors)

proposals
#> Parameters: 
#> * mean: posterior [4.97 ± 0.04]
#> * sd1: posterior [2.00 ± 0.08]
#> * sd2: posterior [0.99 ± 0.04]
#> Constraints:
#> * mean > sd2
```
