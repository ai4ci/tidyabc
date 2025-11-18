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
#>  [1] 0.2523565 0.3592547 0.2030717 0.6602649 0.7783076 0.7992848 0.3085951
#>  [8] 0.7512198 0.7351616 0.7278198
```
