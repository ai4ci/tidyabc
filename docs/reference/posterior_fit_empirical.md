# Fit empirical distribution to posterior samples for generating more waves

This function allows "updating" of the prior with an empirical posterior
distribution which will retain the bounds of the prior, and is
effectively a spline based transform of the prior distribution, based on
the density of data in prior space. This gives a clean density when data
is close to a prior distribution limit and work better than a standard
density.

## Usage

``` r
posterior_fit_empirical(
  posteriors_df,
  priors_list,
  knots = NULL,
  bw = 0.1,
  widen_by = 1
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

- knots:

  the number of knots to model the CDF with. Optional, and will be
  typically inferred from the data size. Small numbers tend to work
  better if we expect the distribution to be unimodal.

- bw:

  for Adaptive ABC data distributions are smoothed before modelling the
  CDF. Over smoothing can reduce convergence rate, under-smoothing may
  result in noisy posterior estimates, and appearance of local modes.
  This is a proportion of the ESS and defaults to 0.1.

- widen_by:

  change the dispersion of proposal distribution in ABC adaptive,
  preserving the median. This is akin to a nonlinear, heteroscedastic
  random walk in the quantile space, and can help address over-fitting
  or local modes in the ABC adaptive waves. `widen_by` is an odds ratio
  and describes how much further from the median any given part of the
  distribution is after transformation. E.g. if the median of a
  distribution is zero, and the `widen_by` is 2 then the 0.75 quantile
  will move to the position of the 0.9 quantile. The distribution will
  stay within the support of the prior. This is by default 1.05 which
  allows for some additional variability in proposals.

## Value

an `abc_prior` S3 object approximating the distribution of the posterior
samples, and adhering to the support of the provided priors.

## Details

This takes weighted posterior samples \\\\(\theta^{(i)}, w^{(i)})\\\\
for each parameter \\\theta_j\\ and constructs an empirical distribution
function \\Q_j(\theta_j)\\ that approximates the posterior marginal for
that parameter.

For each parameter \\\theta_j\\, the empirical distribution \\Q_j\\ is
fitted using the
[`empirical()`](https://ai4ci.github.io/tidyabc/reference/empirical.md)
function. The fitting is performed on the weighted samples
\\(\theta\_{j}^{(i)}, w^{(i)})\\. A key feature is the use of a link
function \\h_j(\cdot)\\ during the fitting process. This link function
is typically taken from the original prior distribution
\\P_j(\theta_j)\\ (i.e., \\h_j = P_j\\). This ensures that the fitted
empirical distribution \\Q_j\\ respects the support constraints defined
by the prior (e.g., positive values for a log-normal prior). The fitting
effectively occurs in the transformed space defined by the link:
\\h_j(\theta_j^{(i)})\\.

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
proposals = posterior_fit_empirical(fit$posteriors, fit$priors)

proposals
#> Parameters: 
#> * mean: posterior [4.98 ± 0.04]
#> * sd1: posterior [2.00 ± 0.09]
#> * sd2: posterior [1.00 ± 0.04]
#> Constraints:
#> * mean > sd2
```
