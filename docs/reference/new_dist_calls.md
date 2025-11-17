# Create a new `dist_fns` S3 object

Create a new `dist_fns` S3 object

## Usage

``` r
new_dist_calls(
  name = "NULL",
  pcall,
  qcall,
  rcall,
  dcall,
  ...,
  params = NULL,
  knots = NULL,
  smooth = TRUE
)
```

## Arguments

- name:

  a name for the distribution

- pcall:

  a call of the cumulative probability function

- qcall:

  a call of the quantile function

- rcall:

  a call of the RNG function

- dcall:

  a call of the density function

- ...:

  parameters for the distribution.

- params:

  parameters for the distribution as a named list.

- knots:

  in empirical distributions this holds details on the knot points

## Value

a new `dist_fns` S3 object or a `dist_fns_list` if parameters are
vectorised

## Unit tests


    tmp = new_dist_calls(name="norm",
      pcall=.qual("pnorm"),
      qcall=.qual("qnorm"),
      rcall=.qual("rnorm"),
      dcall=.qual("dnorm")
    )
    testthat::expect_equal(tmp$p(-2:2), c(
      0.0227501319481792,
      0.158655253931457,
      0.5,
      0.841344746068543,
      0.977249868051821
    ))
