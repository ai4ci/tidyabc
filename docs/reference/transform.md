# Generate a distribution from a truncation of another

Generate a distribution from a truncation of another

## Usage

``` r
transform(link, dist, ..., name = NULL)
```

## Arguments

- link:

  a link function (or name of a link function)

- dist:

  distribution(s) as a name, function or a `dist_fn` S3 object

- ...:

  parameters for the underlying distribution if `dist` is a name or
  function.

- name:

  a name for the link function

## Value

a `dist_fn` or `dist_fn_list` holding the transformed distribution(s)

## Unit tests




    t = transform("log","norm")
    ps = seq(0,1,0.1)
    qs = 0:6
    testthat::expect_equal(t$q(ps),qlnorm(ps))
    testthat::expect_equal(t$p(qs),plnorm(qs))
    testthat::expect_equal(t$d(qs),dlnorm(qs))

    t2 = transform("log","norm", 0.4, 0.1)
    testthat::expect_equal(t2$q(ps),qlnorm(ps,0.4,0.1))
