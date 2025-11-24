# Fit analytical distribution to posterior samples for generating more waves

This function allows "updating" of the prior with a posterior
distribution from the same family as the prior with updated parameters.

## Usage

``` r
posterior_fit_analytical(posteriors_df, priors_list)
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

## Value

an `abc_prior` S3 object approximating the distribution of the posterior
samples from the same family as the prior.

## Details

This takes weighted posterior samples \\\\(\theta^{(i)}, w^{(i)})\\\\
for each parameter \\\theta_j\\ and fits an analytical distribution
function \\Q_j(\theta_j)\\ that approximates the posterior marginal for
that parameter.

The resulting empirical distribution \\Q_j\\ is a `dist_fns` object that
includes the support constraints from the original prior. The set of all
fitted marginal distributions \\\\Q_j\\\\ forms the new proposal list
for the next wave.

Additionally, the weighted covariance matrix \\R\\ of the samples in the
MVN space (defined by the prior copula) is calculated: \$\$ R =
\text{Cov}\_{w}(\Phi^{-1}(P_1(\theta_1)), \dots,
\Phi^{-1}(P_d(\theta_d))) \$\$ where \\\Phi^{-1}\\ is the quantile
function of the standard normal. This covariance matrix is stored as an
attribute (`"cor"`) of the returned proposal list and is used to induce
correlation structure when sampling new proposals from the empirical
distributions in subsequent waves.

## Examples

``` r
fit = example_smc_fit()
proposals = posterior_fit_analytical(fit$posteriors, fit$priors)

proposals
#> Parameters: 
#> * mean: unif(min = 4.9, max = 5.05)
#> * sd1: unif(min = 1.81, max = 2.22)
#> * sd2: unif(min = 0.898, max = 1.09)
#> Constraints:
#> * mean > sd2
```
