# Fit a 1D linear model and predict output

Fit a 1D linear model and predict output

## Usage

``` r
.predict_lm_1d(y, x, new_x = NULL, w = NULL)
```

## Arguments

- y:

  the y values. At least 2.

- x:

  the x values. At least 2.

- new_x:

  the values to predict at (vectorised). If null will return a function
  that can predict new values.

- w:

  weights (optional)

## Value

the predicted value of `y` at `new_x` or a prediction function

## Unit tests


    .predict_lm_1d(c(6,5,7,10), c(1,2,3,4), c(0,1))
    .predict_lm_1d(c(0,1), c(0,1), -2:2)
