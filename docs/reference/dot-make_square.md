# Strictly recycle the items in a list to the same length

This removes NULL entries, and zero length lists

## Usage

``` r
.make_square(lst)
```

## Arguments

- lst:

  a list input

## Value

the list with all items recycled to the same length

## Unit tests


    lst = list(a = 1,b=1:3, c= NULL, d=list(), e=character())
    testthat::expect_equal(
      .make_square(lst),
      list(a = c(1, 1, 1), b = 1:3)
    )
