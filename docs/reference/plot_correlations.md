# A parameter posterior correlation plot

A parameter posterior correlation plot

## Usage

``` r
plot_correlations(posteriors_df, truth = NULL)
```

## Arguments

- posteriors_df:

  a dataframe of posteriors that have been selected by ABC this may
  include columns for scores, weight and/or simulation outputs (
  `abc_component_score`, `abc_summary_distance`, `abc_weight`,
  `abc_simulation` ) as well as columns matching the `priors` input
  specification.

- truth:

  a named numeric vector of known parameter values

## Value

a patchwork of ggplots including density and 2d scatters for each
combination of posteriors.

## Examples

``` r
p = plot_correlations(
  example_adaptive_fit(),
  example_truth()
)

p & ggplot2::theme(
  axis.title.y = ggplot2::element_text(angle=70,vjust=0.5),
  axis.title.x = ggplot2::element_text(angle=20,hjust=0.5)
)

```
