# Quantile from weighted data

This is a moderately expensive function to call (in memory terms), as it
needs to construct the whole quantile function. if there are multiple
calls consider using `empirical_data` to build the quantile function and
using that.

## Usage

``` r
.wquantile(p, x, w = NULL, names = TRUE, ...)
```

## Arguments

- p:

  for CDF fits, a vector the same length as `x` giving `P(X <= x)`.

- x:

  either a vector of samples from a distribution `X` or cut-offs for
  cumulative probabilities when combined with `p`

- w:

  for data fits, a vector the same length as `x` giving the importance
  weight of each observation. This does not need to be normalised. There
  must be some non zero weights, and all must be finite.

- ...:

  Named arguments passed on to
  [`empirical_data`](https://ai4ci.github.io/tidyabc/reference/empirical_data.md)

  `...`

  :   Named arguments passed on to
      [`.logit_z_interpolation`](https://ai4ci.github.io/tidyabc/reference/dot-logit_z_interpolation.md)

      `bw`

      :   a bandwidth expressed in terms of the probability width, or
          proportion of observations.

## Value

a vector of quantiles

## Unit tests



    test = function(rfn,qfn, link,..., n = 100000, tol=1000/(n+10000)) {
      testthat::expect_equal(abs(
       unname(.wquantile(c(0.025, 0.5, 0.975),rfn(100000,...),link=link,names=FALSE)-
        qfn(c(0.025, 0.5, 0.975), ...))
      ), c(0,0,0), tolerance=tol)
    }

    withr::with_seed(123, {

      test(rnorm,qnorm,"identity",n = 10000)
      test(rnorm,qnorm,"identity",mean=4,n = 10000)
      test(rnorm,qnorm,"identity",sd=3, n = 100000, tol=0.05)

      test(rnorm,qnorm,"identity",n = 5000)
      test(rnorm,qnorm,"identity",n = 1000)
      test(rnorm,qnorm,"identity",n = 100)
      test(rnorm,qnorm,"identity",n = 30)

      test(rgamma,qgamma,"log", 4,n = 10000)
      test(rgamma,qgamma,"log", 4,n = 5000)
      test(rgamma,qgamma,"log", 4,n = 1000)
      test(rgamma,qgamma,"log", 4, 3,n = 100)
      test(rgamma,qgamma,"log", 4,n = 30)

       test(runif,qunif,as.link_fns(c(0,10)),0,10)

    })

    .wquantile(0.5, rnorm(1000))
