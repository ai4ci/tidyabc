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
#>  [1] 0.05615852 0.56645297 0.28313048 0.47006282 0.39084400 0.38064058
#>  [7] 0.61124899 0.46977093 0.81431682 0.36195327
```
