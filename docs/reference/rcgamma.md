# Sampling: gamma distribution constrained to have mean \> sd

Density, distribution function, quantile function and random generation
for the Gamma distribution with parameters `shape` and `scale`.

## Usage

``` r
rcgamma(n, mean, kappa = 1/mean)
```

## Arguments

- n:

  number of observations

- mean:

  the mean value on the true scale (vectorised)

- kappa:

  a coefficient of variation. where 0 is no variability and 1 is
  maximally variability (vectorised)

## Value

`dgamma` gives the density, `pgamma` gives the distribution function,
`qgamma` gives the quantile function, and `rgamma` generates random
deviates.

Invalid arguments will result in return value `NaN`, with a warning.

The length of the result is determined by `n` for `rgamma`, and is the
maximum of the lengths of the numerical arguments for the other
functions.

The numerical arguments other than `n` are recycled to the length of the
result. Only the first elements of the logical arguments are used.

## See also

[`stats::rgamma()`](https://rdrr.io/r/stats/GammaDist.html)

## Examples

``` r
rcgamma(10, 2, 0.5)
#>  [1] 0.7464527 2.4585756 3.2335808 2.1612793 0.7633470 0.7181991 2.8403559
#>  [8] 1.2380553 6.6017349 0.2844472
```
