# Create a `dist_fns` S3 object

A class wrapping a single (or set) of parametrised distributions and
provides access to the quantile, cumulative probability and random
functions of that specific distribution. Parametrisation is handled on
construction.

## Usage

``` r
# S3 method for class 'character'
as.dist_fns(x, ...)

# S3 method for class '`function`'
as.dist_fns(x, ...)

# S3 method for class 'fitdist'
as.dist_fns(x, ...)

as.dist_fns(x, ...)
```

## Arguments

- x:

  a `dist_fns` S3 object

- ...:

  passed onto methods

## Value

a `dist_fns` S3 object

## Details

`dist_fns` and `dist_fns_list` objects support `$` access for fields and
`@` access for attributes. `dist_fns_list`s can be made with the
[`c()`](https://rdrr.io/r/base/c.html) or
[`rep()`](https://rdrr.io/r/base/rep.html) functions, or with the
`purrr` style map functions, and they support subsetting. Individual
`dist_fns` members of `dist_fns_list`s can be accessed with `[[`.

## Methods (by class)

- `as.dist_fns(character)`: Construct a distribution by name

- `` as.dist_fns(`function`) ``: From a statistical function

- `as.dist_fns(fitdist)`: From a
  [`fitdistrplus::fitdist`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html)
  output

## Unit tests


    pois = as.dist_fns("pois",lambda = 8)
    n = as.dist_fns("norm",mean=4)
