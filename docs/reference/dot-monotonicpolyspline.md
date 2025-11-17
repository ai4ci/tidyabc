# Create a `polySpline` object from a `stats::splinefun` call.

Create a `polySpline` object from a
[`stats::splinefun`](https://rdrr.io/r/stats/splinefun.html) call.

## Usage

``` r
.monotonicpolyspline(x, y)
```

## Arguments

- x, y:

  vectors giving the coordinates of the points to be interpolated.
  Alternatively a single plotting structure can be specified: see
  [`xy.coords`](https://rdrr.io/r/grDevices/xy.coords.html).

  `y` must be increasing or decreasing for `method = "hyman"`.

## Value

An object of class `polySpline`.

## Details

This function converts the output of
[`stats::splinefun`](https://rdrr.io/r/stats/splinefun.html), to a a
polynomial spline object from the `splines` package.

## Unit tests



    # strictly increasing
    spl = .monotonicpolyspline(1:10, log(1:10))
    testthat::expect_equal(
      predict(spl,4.5)$y,
      log(4.5),
      tolerance = 0.01
    )

    # strictly decreasing
    spl2 = .monotonicpolyspline(1:10, -log(1:10))
    testthat::expect_equal(
      predict(spl2,4.5)$y,
      -log(4.5),
      tolerance = 0.01
    )

    # not monotonic
    testthat::expect_error(
      {
        .monotonicpolyspline(-4:4, (-4:4)^2)
      },
      "Data is not monotonic.",
      fixed = TRUE
    )
