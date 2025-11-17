# Calculate the excess kurtosis of a set of data

Calculate the excess kurtosis of a set of data

## Usage

``` r
kurtosis(x, na.rm = FALSE, excess = TRUE)
```

## Arguments

- x:

  a vector of observations

- na.rm:

  remove `NA`s?

- excess:

  if false calculates raw kurtosis rather than excess

## Value

the excess kurtosis

## Examples

``` r
kurtosis(stats::rnorm(1000))
#> [1] -0.2308395
kurtosis(stats::rpois(1000, 2)) # leptokurtic > 0 (usually)
#> [1] 0.1414471
kurtosis(stats::runif(1000)) # platykurtic: < 0
#> [1] -1.191496

kurtosis(stats::rlnorm(1000))
#> [1] 32.07764
```
