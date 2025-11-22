# Weighted bandwidth selector

Weighted bandwidth selector

## Usage

``` r
wbw.nrd(x, w = NULL)
```

## Arguments

- x:

  data

- w:

  weights

## Value

a bandwidth based on weighted data

## Examples

``` r
x = runif(1000)
w = rgamma(1000,2)

wbw.nrd(x)
#> [1] 0.07698208
wbw.nrd(x,w)
#> [1] 0.07698208
```
