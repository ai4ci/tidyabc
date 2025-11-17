# Fit a weighted 1D linear model and predict output

If all the weights are NA they are ignored.

## Usage

``` r
.fit_lm_1d(y, x, ..., w = NULL, new_x = NULL)
```

## Arguments

- y:

  the y values. At least 2.

- x:

  the x values. At least 2.

- ...:

  must be empty

- w:

  weights (optional)

- new_x:

  (optional) the values to predict at (vectorised).

## Value

either a vector with intercept (`c`) and gradient (`m`) or a vector of
predictions at points `new_x`

## Unit tests


    testthat::expect_equal(
      .fit_lm_1d(c(1, 2, 3, 4), c(5, 6, 7, 8)),
      c(c = -4, m = 1)
    )

    testthat::expect_equal(
       .fit_lm_1d(c(0, 1), c(0, 1)),
       c(c = 0, m = 1)
    )

    withr::with_seed(123,{
      x = 1:100
      y = 2*x+3+stats::runif(100,-0.1,0.1)
      testthat::expect_equal(
        .fit_lm_1d(y, x),
        c(c = 3, m = 2),
        tolerance = 0.01
      )
    })

    # with prediction:
    testthat::expect_error(
      {
        .fit_lm_1d(c(6, 5, 7, 10), c(1, 2, 3, 4), c(0, 1))
      },
      structure("`...` must be empty.", names = ""),
      fixed = TRUE
    )

    testthat::expect_equal(
      .fit_lm_1d(c(0, 1), c(0, 1), new_x = -2:2),
      c(-2, -1, 0, 1, 2)
    )
