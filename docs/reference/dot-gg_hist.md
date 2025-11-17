# Simple one liner plots from vectors with defaults

Simple one liner plots from vectors with defaults

## Usage

``` r
.gg_hist(x, binwidth = 1, ...)

.gg_dens(x, bw = NULL, ...)

.gg_bar(x, ...)

.gg_scatter(x, y, ...)

.gg_denshist(x, y, binwidth = 1)

.gg_linebar(x, y, binwidth = 1)
```

## Arguments

- x:

  a vector of data or counts

- binwidth:

  histogram binwidth

- ...:

  passed onto the geoms

- bw:

  bandwidth parameter,if given it is interpreted as a fraction of the sd
  of the data

- y:

  a vector or matrix of data or counts

## Value

a ggplot

## Functions

- `.gg_hist()`: Simple histogram

- `.gg_dens()`: A density plot

- `.gg_bar()`: A bar chart of counts

- `.gg_scatter()`: X-Y scatter plot

- `.gg_denshist()`: Histogram with densities overlaid.

- `.gg_linebar()`: Bar chart with lines overlaid
