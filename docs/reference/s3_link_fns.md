# Create a `link_fns` S3 object

The link function class allows forwards and backwards transformation.
Link functions can be defined by name or using a statistical
distribution in which case the forward link is a logit of the cumulative
probability and the reverse is the quantile of the expit.

## Arguments

- x:

  a `link_fns` S3 object

- ...:

  passed onto methods

## Details

`link_fns` and `link_fns_list` objects support `$` access for fields and
`@` access for attributes. `link_fns_list`s can be made with the
[`c()`](https://rdrr.io/r/base/c.html) or
[`rep()`](https://rdrr.io/r/base/rep.html) functions, or with the
`purrr` style map functions, and they support subsetting. Individual
`link_fns` members of `link_fns_list`s can be accessed with `[[`.
