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
#>  [1] 0.63555887 0.38488804 0.14717321 0.46494021 0.66500435 0.39267154
#>  [7] 0.09835966 0.76532015 0.75204022 0.77632768
```
