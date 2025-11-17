# Plot a `dist_fns_list` S3 object

Plot a smoothed version of the PDFs of a set of `dist_fns`. These are
`ggplots` that can be facetted by `names`, `id` or `group`

## Usage

``` r
# S3 method for class 'dist_fns_list'
plot(
  x,
  ...,
  mapping = .gg_check_for_aes(...),
  steps = 200,
  tail = 0.001,
  plot_quantiles = TRUE,
  smooth = TRUE
)
```

## Arguments

- x:

  a `dist_fns_list`

- ...:

  passed to
  [`ggplot2::geom_step`](https://ggplot2.tidyverse.org/reference/geom_path.html),[`ggplot2::geom_rect`](https://ggplot2.tidyverse.org/reference/geom_tile.html)
  or
  [`ggplot2::geom_area`](https://ggplot2.tidyverse.org/reference/geom_ribbon.html),

- mapping:

  override default aesthetics with `name`, `id` or `group`

- steps:

  resolution of the plot

- tail:

  the minimum tail probability to plot

- plot_quantiles:

  by default the quantiles of the distribution are plotted over the
  density sometimes this makes it hard to read

- smooth:

  by default some additional smoothing is used to cover up small
  discontinuities in the PDF.

## Value

a ggplot
