# Quantile from weighted data with link function support

This quantile function has different order of parameters from base
quantile. It takes a weight and a link function specification which
allows us to define the support of the quantile function. It is
optimised for imputing the tail of distributions and not speed.

## Usage

``` r
wquantile(p, x, w = NULL, link = "identity", names = TRUE, window = 7)
```

## Arguments

- p:

  the probabilities for which to estimate quantiles from the data

- x:

  a set of observations

- w:

  for data fits, a vector the same length as `x` giving the importance
  weight of each observation. This does not need to be normalised. There
  must be some non zero weights, and all must be finite.

- link:

  a link function. Either as name, or a `link_fns` S3 object. In the
  latter case this could be derived from a statistical distribution by
  `as.link_fns(<dist_fns>)`. This supports the use of a prior to define
  the support of the empirical function, and is designed to prevent tail
  truncation. Support for the updated quantile function will be the same
  as the provided prior.

- names:

  name the resulting quantile vector

- window:

  the number of data points to include when estimating the quantile. The
  closest `window` points are picked and used as a distance weighted
  linear interpolation of the weighted CDF in logit-link space. This
  tends to give good results for extrapolating tails.

## Value

a vector of quantiles

## Details

This is a moderately expensive function to call (in memory terms), as it
needs to construct the whole quantile function. if there are multiple
calls consider using
[`empirical()`](https://ai4ci.github.io/tidyabc/reference/empirical.md)
to build a quantile function and using that.

## Unit tests



    test = function(rfn,qfn, link,..., n = 100000, tol=1000/(n+10000)) {
      testthat::expect_equal(abs(
       unname(wquantile(c(0.025, 0.5, 0.975),rfn(100000,...),link=link,names=FALSE)-
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

## Examples

``` r
# unweighted:
wquantile(p = c(0.25,0.5,0.75), x = stats::rnorm(1000))
#>         25%         50%         75% 
#> -0.68610780 -0.04486547  0.63457453 

# weighted:
wquantile(p = c(0.25,0.5,0.75), x = seq(-2,2,0.1), w = stats::dnorm(seq(-2,2,0.1)))
#>           25%           50%           75% 
#> -6.384711e-01  2.591425e-17  6.384711e-01 
```
