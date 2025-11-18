# Weighted mean

a simple alias for base `weighted,mean`

## Usage

``` r
wmean(x, w = NULL, na.rm = TRUE)
```

## Arguments

- x:

  either a vector of samples from a distribution `X` or cut-offs for
  cumulative probabilities when combined with `p`

- w:

  for data fits, a vector the same length as `x` giving the importance
  weight of each observation. This does not need to be normalised. There
  must be some non zero weights, and all must be finite.

- na.rm:

  remove NAs (default TRUE)

## Value

a standard deviation

## Examples

``` r
#' # unweighted:
wmean(x = stats::rnorm(1000))
#> [1] -0.01197784

# weighted:
wmean(x = seq(-2,2,0.1), w = stats::dnorm(seq(-2,2,0.1)))
#> [1] 3.904271e-17
```
