# Create a `dist_fns` S3 object

A class wrapping a single (or set) of parametrised distributions and
provides access to the quantile, cumulative probability and random
functions of that specific distribution. Parametrisation is handled on
construction.

## Arguments

- x:

  a `dist_fns` S3 object

- ...:

  passed onto methods

## Details

`dist_fns` and `dist_fns_list` objects support `$` access for fields and
`@` access for attributes. `dist_fns_list`s can be made with the
[`c()`](https://rdrr.io/r/base/c.html) or
[`rep()`](https://rdrr.io/r/base/rep.html) functions, or with the
`purrr` style map functions, and they support subsetting. Individual
`dist_fns` members of `dist_fns_list`s can be accessed with `[[`.
