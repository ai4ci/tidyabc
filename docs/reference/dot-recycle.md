# Strictly recycle function parameters

`.recycle` is called within a function and ensures the parameters in the
calling function are all the same length by repeating them using `rep`.
This function alters the environment from which it is called. It is
stricter than R recycling in that it will not repeat vectors other than
length one to match the longer ones, and it throws more informative
errors.

## Usage

``` r
.recycle(..., .min = 1, .env = rlang::caller_env())
```

## Arguments

- ...:

  the variables to recycle

- .min:

  the minimum length of the results (defaults to 1)

- .env:

  the environment to recycle within.

## Value

the length of the longest variable

## Details

NULL values are not recycled, missing values are ignored.

## Unit tests


    testfn = function(a, b, c) {
      n = .recycle(a,b,c)
      return(list(
        a=a, b=b, c=c, n=n
      ))
    }

    tmp = testfn(a=c(1,2,3), b="needs recycling", c=NULL)

    testthat::expect_equal(tmp$n, 3)
    testthat::expect_null(tmp$c)
    testthat::expect_equal(length(tmp$a), length(tmp$b))

    # no parameter
    testthat::expect_error(testfn(a=c(1,2,3), c=NULL))

    tmp = testfn(a=character(), b=integer(), c=NULL)

    testthat::expect_equal(tmp$n, 0)

    # inconsistent to have a zero length and a non zero length
    testthat::expect_error(testfn(a=c("a","b"), b=integer(), c=NULL))

    # This is currently unsupported
    #testfn = function(a, b, c, ...) {
    # n = .recycle(a,b,c, ...)
    # return(list(
    #   a=a, b=b, c=c, n=n, ...
    # ))
    #}
    #
    #tmp = testfn(a=c(1,2,3), b="needs recycling", c=NULL, d="additional")
