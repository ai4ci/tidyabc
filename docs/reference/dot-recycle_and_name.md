# Recycles an input to match the length a set of names

Recycles an input to match the length a set of names

## Usage

``` r
.recycle_and_name(input, names)
```

## Arguments

- input:

  an input vector which may be named

- names:

  a set of names to match

## Value

a recycled vector with names matching the order of `names`

## Unit tests



    testthat::expect_equal(
      .recycle_and_name(c(a = 1, b = 2, c = 3), c("c", "b", "a")),
      c(c = 3, b = 2, a = 1)
    )

    testthat::expect_equal(
      .recycle_and_name(1, c("A", "B")),
      c(A = 1, B = 1)
    )
