# Calculate the skew of a set of data

Calculate the skew of a set of data

## Usage

``` r
skew(x, na.rm = FALSE)
```

## Arguments

- x:

  a vector of observations

- na.rm:

  remove `NA`s?

## Value

the skew

## Examples

``` r
skew(stats::rnorm(1000))
#> [1] 0.008772958
skew(stats::rbeta(1000, 1, 8)) # positively (left) skewed
#> [1] 1.38873
skew(stats::rbeta(1000, 8, 1)) # negatively (right) skewed
#> [1] -1.481538

skew(stats::rlnorm(1000))
#> [1] 4.988768
```
