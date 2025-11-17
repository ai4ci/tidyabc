# `abc_prior` S3 class

`abc_prior` S3 class

## Usage

``` r
new_abc_prior(.dists, .constraints = list(), .derived = list(), .cor = NULL)

as.abc_prior(x, ...)

is.abc_prior(x, ...)

# S3 method for class 'abc_prior'
format(x, ...)

# S3 method for class 'abc_prior'
print(x, ...)

# S3 method for class 'abc_prior'
plot(x, ...)
```

## Arguments

- .dists:

  distribution functions as a named list of S3 `dist_fns` objects

- .constraints:

  a list of one sided formulae the result each of which should evaluate
  to a boolean when compared against the names of the priors and derived
  values.

- .derived:

  a list of two sided formulae. The RHS refer to the priors, and the LHS
  as a name to derive.

- .cor:

  (optional) a correlation matrix for the priors

- x:

  an `abc_prior` S3 object

## Value

an S3 object of class `abc_prior` which contains

- a list of `dist_fns`

- a `cor` attribute describing their correlation

- a `derived` attribute describing derive values

- a `constraints` attribute listing the constraints

- a `params` attribute listing the names of the parameters

## Methods (by generic)

- `format(abc_prior)`: Format an `abc_prior`

- `print(abc_prior)`: Print an `abc_prior`

- `plot(abc_prior)`: Plot an `abc_prior`

## Functions

- `new_abc_prior()`: Create a new prior

- `as.abc_prior()`: Create a prior from a named list of `dist_fns`

- `is.abc_prior()`: Test is an `abc_prior`

## Unit tests


    p = new_abc_prior(
      .dists = list(
        mean = as.dist_fns("norm",4,2),
        sd = as.dist_fns("gamma",2)
      ),
      .derived = list(
        shape ~ mean^2 / sd^2,
        rate ~ mean / sd^2
      ),
      .constraints = list(
        ~ mean > sd
      )
    )

    testthat::expect_equal(
      format(p),
      "Parameters: \n* mean: norm(mean = 4, sd = 2)\n* sd: gamma(shape = 2, rate = 1)\nConstraints:\n* mean > sd\nDerived values:\n* shape = mean^2/sd^2\n* rate = mean/sd^2"
    )
