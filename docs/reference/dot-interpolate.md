# Interpolate from 2 coordinates and predict output

Interpolate from 2 coordinates and predict output

## Usage

``` r
.interpolate(y, x, new_x = NULL)
```

## Arguments

- y:

  the y values. length 2.

- x:

  the x values. length 2.

- new_x:

  the values to predict at (vectorised). If null will return a function
  that can predict new values.

## Value

the interpolated value or a prediction function

## Unit tests


    .interpolate(c(0,1),c(0,1),5)
