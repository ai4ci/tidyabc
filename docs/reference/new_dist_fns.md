# Create a new `dist_fns` S3 object

Create a new `dist_fns` S3 object

## Usage

``` r
new_dist_fns(
  name = "NULL",
  pfn = pnull,
  qfn = qnull,
  rfn = NULL,
  dfn = NULL,
  ...,
  params = NULL,
  knots = NULL,
  smooth = TRUE
)
```

## Arguments

- name:

  a name for the distribution

- pfn:

  the cumulative probability function

- qfn:

  the quantile function

- rfn:

  the RNG function

- dfn:

  the density function

- ...:

  parameters for the distribution.

- params:

  parameters for the distribution as a named list.

- knots:

  in empirical distributions this holds details on the knot points

## Value

a new `dist_fns` S3 object or a `dist_fns_list` if parameters are
vectorised
