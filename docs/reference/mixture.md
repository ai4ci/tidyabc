# Construct a mixture distribution

Construct a mixture distribution

## Usage

``` r
mixture(dists, weights = 1, steps = 200, tail_p = 1e-04, ...)
```

## Arguments

- dists:

  a `dist_fn_list` of distribution functions

- weights:

  a vector of weights

- steps:

  the number of points that the mixture distribution is evaluated at to
  construct the empirical mixture

- tail_p:

  the support fo the tail of the distribution

- ...:

  Named arguments passed on to
  [`empirical_cdf`](https://ai4ci.github.io/tidyabc/reference/empirical_cdf.md)

  `smooth`

  :   fits the empirical distribution with a pair of splines for CDF and
      quantile function, creating a mostly smooth density. This
      smoothness comes at the price of potential over-fitting and will
      produce small differences between `p` and `q` functions such that
      `x=p(q(x))` is no longer exactly true. Setting this to false will
      replace this with a piecewise linear fit that is not smooth in the
      density, but is exact in forward and reverse transformation.

## Value

a `dist_fn` of the mixture distribution

## Unit tests



    dists = c(
      as.dist_fns("norm", mean=-1),
      as.dist_fns("norm", mean=1),
      as.dist_fns("gamma", shape=2)
    )
    weights = c(1,1,0.3)

    mix = mixture(dists,weights)
    testthat::expect_equal(
      format(mix),
      "mixture; Median (IQR) 0.289 [-0.886 — 1.31]"
    )
