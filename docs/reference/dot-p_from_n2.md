# Simple sigmoid tending towards a minimum proportion of input

Simple sigmoid tending towards a minimum proportion of input

## Usage

``` r
.p_from_n2(n, n_inf, p_max = 1, n_100 = 0)
```

## Arguments

- n:

  a vector of lengths

- p_max:

  the maximum value at `n <= n_100`

- n_100:

  an offset under which the value is 1 (default is 0)

## Value

a vector of proportions for each length

## Unit tests


    testthat::expect_equal(.p_from_n2(0:5, 1), c(
      1,
      0.707106781186547,
      0.447213595499958,
      0.316227766016838,
      0.242535625036333,
      0.196116135138184
    ))
