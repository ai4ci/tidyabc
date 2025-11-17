# Create a `link_fns` S3 object

The link function class allows forwards and backwards transformation.
Link functions can be defined by name or using a statistical
distribution in which case the forward link is a logit of the cumulative
probability and the reverse is the quantile of the expit.

## Usage

``` r
# S3 method for class 'character'
as.link_fns(x, ...)

# S3 method for class 'dist_fns'
as.link_fns(x, ...)

# S3 method for class 'family'
as.link_fns(x, ...)

# S3 method for class 'numeric'
as.link_fns(x, ..., na.rm = TRUE)

as.link_fns(x, ...)
```

## Arguments

- x:

  a range of values that for the support

- ...:

  ignored

## Value

a `link_fns` S3 object

## Details

`link_fns` and `link_fns_list` objects support `$` access for fields and
`@` access for attributes. `link_fns_list`s can be made with the
[`c()`](https://rdrr.io/r/base/c.html) or
[`rep()`](https://rdrr.io/r/base/rep.html) functions, or with the
`purrr` style map functions, and they support subsetting. Individual
`link_fns` members of `link_fns_list`s can be accessed with `[[`.

## Methods (by class)

- `as.link_fns(character)`: Link function from name

- `as.link_fns(dist_fns)`: Link function from name

- `as.link_fns(family)`: Link function from name

- `as.link_fns(numeric)`: Link function from support vector

## Unit tests


    links = c("ident", "log", "logit", "probit", "cloglog", "neginv", "inv2")
    test = seq(0.1,0.9,0.1) # within support of all links
    for (l in links) {
      lfn = as.link_fns(l)
      t = lfn$trans(test)
      i = lfn$inv(t)
      testthat::expect_equal(i,test)
    }
