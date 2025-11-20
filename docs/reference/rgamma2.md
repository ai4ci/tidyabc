# The Gamma Distribution

Density, distribution function, quantile function and random generation
for the Gamma distribution with parameters `shape` and `scale`.

## Usage

``` r
rgamma2(n, mean, sd = sqrt(mean))
```

## Arguments

- n:

  number of observations

- mean:

  the mean value on the true scale (vectorised)

- sd:

  the standard deviation on the true scale (vectorised)

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
rgamma2(10, 2, 1)
#>  [1] 1.361138 0.743578 2.482845 1.144367 3.898707 2.693990 2.058435 4.201515
#>  [9] 1.958853 1.356946
```
