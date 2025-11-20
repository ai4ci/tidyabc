# Increase the dispersion of a distribution

Increases the dispersion (spread) of a distribution by transforming its
quantile function in a standardized Q-Q space.

## Usage

``` r
widen(x, scale, knots = NULL, name = NULL)
```

## Arguments

- x:

  a distribution as a `dist_fns` S3 object

- scale:

  acts as a multiplier for the log-odds difference from the median,
  effectively acting like an odds-ratio parameter. A `scale > 1`
  increases dispersion (stretches quantiles away from the median in
  logit space), while `0 < scale < 1` decreases dispersion (compresses
  quantiles towards the median in logit space). The median value is
  preserved in the original parameter space.

- knots:

  the number of knots in the transformed distribution, if it is not
  already an empirical CDF distribution.

- name:

  a name for the widened distribution

## Value

an empirical `dist_fn` with the same median and increased SD. This
transformation will change the mean of skewed distributions.

## Details

The transformation aims to increase the standard deviation by a factor
`scale` while preserving the median. It operates on the internal Q-Q
space representation (`qx`, `qy`) of an empirical distribution generated
by `empirical_cdf`.

Applies a logit-space scaling transformation centred on the median
quantile. This transformation modifies the quantile axis (`qx`) in the
Q-Q space of an empirical distribution to change its dispersion while
preserving the median value.

Let `qx` be the original quantile coordinate in `[0, 1]`, and `qmedian`
be the quantile corresponding to the median (e.g., `qx_from_qy(0.5)`).
The transformation is defined as:

\$\$ qx_2 = \text{expit}\left( (\text{logit}(qx) -
\text{logit}(qmedian)) \times \text{scale} + \text{logit}(qmedian)
\right) \$\$

where \\\text{logit}(x) = \log(x / (1 - x))\\ and \\\text{expit}(x) = 1
/ (1 + \exp(-x))\\ are the standard logit and inverse-logit functions,
respectively.

## Examples

``` r
d1 = as.dist_fns("norm",4,2)
w1 = widen(d1, scale=1.5)


```
