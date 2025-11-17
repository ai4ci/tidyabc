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
#>  [1] 1.5060685 0.5910633 1.8002079 0.6457146 2.4258825 0.8052649 1.2168476
#>  [8] 4.6360830 2.5598636 3.0773423
```
