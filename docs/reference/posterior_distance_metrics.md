# Generate a set of metrics from component scores

The component scores are summary statistics output by the user supplied
`scorer_fn` as a named list. These can be variable in scale and location
and various options exist for combining them. They may need to be
weighted by scale as well as importance to get a model that works well.
Such weights can be input into the ABC algortihms using the
`scoreweights` parameter. This function helps provide diagnostics for
calibrating the `scoreweights` parameter.

## Usage

``` r
posterior_distance_metrics(posteriors_df, obsscores = NULL, keep_data = FALSE)
```

## Arguments

- posteriors_df:

  a dataframe of posteriors that have been selected by ABC this may
  include columns for scores, weight and/or simulation outputs (
  `abc_component_score`, `abc_summary_distance`, `abc_weight`,
  `abc_simulation` ) as well as columns matching the `priors` input
  specification.

- obsscores:

  Summary scores for the observational data. This will be a named list,
  and is equivalent to the output of `scorer_fn`, on the observed data.
  If not given typically it will be assumed to be all zeros.

- keep_data:

  mainly for internal use this flag gives you the component scores as a
  matrix

## Value

a list containing the following items:

- `obsscores`: the input reference scores for each component

- `means`, `sds`: the means and sds of each score component

- `cov`: the covariance matrix for the scores

- `mad`: the mean absolute differences between the `simscores` and the
  `obsscores`

- `rmsd`: the root mean squared differences between the `simscores` and
  the `obsscores`

- `scoreweights`: a the `sds` divided by the `rmsd`. This weight should
  mean that the weights of the individual summary scores have similar
  influence in the overall `abc_summary_distance` output once combined
  during SMC and adaptive waves, especially if euclidean distances are
  involved.

- `simscores`: (if `keep_data`) a matrix of all the scores from these
  input simulation posteriors

- `deltascores`: (if `keep_data`) a matrix of the differences between
  `simscores` and `obsscores`.

## Details

Given a list of component scores \\S_s^{(i)} = (s_1^{(i)}, \dots,
s_m^{(i)})\\ from \\n\\ simulations and observed scores \\S_o =
(s_1^{(o)}, \dots, s_m^{(o)})\\, this function calculates:

- The mean \\\mu_j\\ and standard deviation \\\sigma_j\\ of each score
  component \\j\\: \$\$ \mu_j = \frac{1}{n}\sum\_{i=1}^n s_j^{(i)},
  \quad \sigma_j = \sqrt{\frac{1}{n-1}\sum\_{i=1}^n (s_j^{(i)} -
  \mu_j)^2} \$\$

- The covariance matrix \\\Sigma\\ of the score components.

- The root mean squared difference (RMSD) between simulated and observed
  scores for each component: \$\$ \text{RMSD}\_j =
  \sqrt{\frac{1}{n}\sum\_{i=1}^n (s_j^{(i)} - s_j^{(o)})^2} \$\$

- A recommended vector of `scoreweights` \\w_j\\, calculated as the
  ratio of the component's standard deviation to its RMSD, normalized to
  sum to 1: \$\$ w_j = \frac{\sigma_j / \text{RMSD}\_j}{\sum\_{k=1}^m
  (\sigma_k / \text{RMSD}\_k)} \$\$ These weights help balance the
  influence of different summary statistics in the overall distance
  metric, especially when using Euclidean distance.

## Examples

``` r
fit = example_rejection_fit()
metrics = posterior_distance_metrics(fit$posteriors)

# other elements available but this is the most important and tells you what
# the relative sizes of the component scores from `scorer_fn` in this sample.
# If this is a sample from the prior then this gives us a way to judge the
# most appropriate relative weighting of each component:
metrics$scoreweights
#>         A         B 
#> 0.5567914 0.4432086 
```
