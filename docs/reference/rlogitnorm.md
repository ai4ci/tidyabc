# Logit-normal distribution

The logit-normal distribution has a support of 0 to 1.

## Usage

``` r
rlogitnorm(n, meanlogit = 0, sdlogit = 1)
```

## Arguments

- n:

  number of observations

- meanlogit:

  the mean on the logit scale

- sdlogit:

  the sd on the logit scale

## Value

a vector of probabilities, quantiles, densities or samples.

## Examples

``` r
rlogitnorm(10, 0, 1)
#>  [1] 0.09704084 0.40588807 0.32105543 0.36672116 0.27289431 0.57001995
#>  [7] 0.33375514 0.88261186 0.33239456 0.69512937
```
