# Plot the evolution of the density function by wave for SMC and adaptive ABC

Plot the evolution of the density function by wave for SMC and adaptive
ABC

## Usage

``` r
plot_evolution(fit, truth = NULL, ..., what = c("posteriors", "proposals"))
```

## Arguments

- fit:

  A S3 `abc_fit` object as output by the `abc_XXX` functions

- truth:

  a named numeric vector of known parameter values

- ...:

  passed on to methods Named arguments passed on to
  [`plot.dist_fns_list`](https://ai4ci.github.io/tidyabc/reference/plot.dist_fns_list.md)

  `mapping`

  :   override default aesthetics with `name`, `id` or `group`

  `steps`

  :   resolution of the plot

  `tail`

  :   the minimum tail probability to plot

  `plot_quantiles`

  :   by default the quantiles of the distribution are plotted over the
      density sometimes this makes it hard to read

  `smooth`

  :   by default some additional smoothing is used to cover up small
      discontinuities in the PDF.

- what:

  plot posterior densities or proposal distributions?

## Value

a plot of the density functions by wave

## Examples

``` r
plot_evolution(
  example_adaptive_fit(),
  example_truth()
)

```
