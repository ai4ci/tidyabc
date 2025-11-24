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

The process involves:

1.  Link transformation: \\x\\ values are transformed using the link
    function: \\x_1 = T(x)\\.

2.  Standardization: Transformed values are standardized: \\x_2 =
    \frac{x_1 - \mu\_{w,1}}{\sigma\_{w,1}}\\, where \\\mu\_{w,1}\\ and
    \\\sigma\_{w,1}\\ are the weighted mean and standard deviation of
    \\x_1\\.

3.  Weighted CDF calculation: The empirical CDF \\y\\ is calculated from
    weights \\w\\.

4.  Logit transformation: \\y\\ is transformed: \\y_2 =
    \text{logit}(y)\\.

5.  Local interpolation: For a target probability \\p\\, \\p_2 =
    \text{logit}(p)\\ is calculated. A window of \\window\\ points is
    selected from the \\(y_2, x_2)\\ pairs around \\p_2\\. A weighted
    linear model is fitted using Gaussian kernel weights based on
    distance in \\y_2\\ space: \\K = \exp(-\frac{1}{2} u^2)\\, where
    \\u\\ is the normalized distance.

6.  Quantile estimation: The local model predicts \\q_2\\ for \\p_2\\.

7.  Back-transformation: The quantile is transformed back: \\q =
    T^{-1}(q_2 \cdot \sigma\_{w,1} + \mu\_{w,1})\\.

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

      test(stats::rnorm,stats::qnorm,"identity",n = 10000)
      test(stats::rnorm,stats::qnorm,"identity",mean=4,n = 10000)
      test(stats::rnorm,stats::qnorm,"identity",sd=3, n = 100000, tol=0.05)

      test(stats::rnorm,stats::qnorm,"identity",n = 5000)
      test(stats::rnorm,stats::qnorm,"identity",n = 1000)
      test(stats::rnorm,stats::qnorm,"identity",n = 100)
      test(stats::rnorm,stats::qnorm,"identity",n = 30)

      test(stats::rgamma,stats::qgamma,"log", 4,n = 10000)
      test(stats::rgamma,stats::qgamma,"log", 4,n = 5000)
      test(stats::rgamma,stats::qgamma,"log", 4,n = 1000)
      test(stats::rgamma,stats::qgamma,"log", 4, 3,n = 100)
      test(stats::rgamma,stats::qgamma,"log", 4,n = 30)

      test(stats::runif,stats::qunif,as.link_fns(c(0,10)),0,10)

    })

## Examples

``` r
# fit weighted data
samples = seq(0,10,0.01)
weights = dgamma2(samples, mean=5, sd=2)

wquantile(c(0.25,0.5,0.75), x = samples, w = weights, link="log")
#>      25%      50%      75% 
#> 3.523936 4.691167 6.067214 

# compared to the sampled distribution
qgamma2(c(0.25,0.5,0.75), mean=5, sd=2)
#> [1] 3.547199 4.736011 6.166153

# unweighted:
wquantile(p = c(0.25,0.5,0.75), x = stats::rnorm(1000))
#>         25%         50%         75% 
#> -0.67790098  0.03865142  0.66577384 
qnorm(p = c(0.25,0.5,0.75))
#> [1] -0.6744898  0.0000000  0.6744898

```
