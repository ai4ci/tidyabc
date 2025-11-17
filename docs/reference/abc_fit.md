# `abc_fit` S3 class

A class holding the output of a single ABC model fitting, either as a
single ABC rejection round or after a set of SMC waves.

## Usage

``` r
new_abc_fit(type, iterations, converged, priors_list, wave_df, summ_df, sim_df)

# S3 method for class 'abc_fit'
format(x, ...)

# S3 method for class 'abc_fit'
summary(object, ..., truth = NULL)

tidy.abc_fit(x, ...)

# S3 method for class 'abc_fit'
print(x, ...)

# S3 method for class 'abc_fit'
plot(x, ..., truth = NULL)
```

## Arguments

- type:

  the type of ABC algorithm

- iterations:

  number of completed iterations

- converged:

  boolean - did the result meet convergence criteria

- wave_df:

  a list of dataframes of wave convergence metrics

- summ_df:

  the summary of the parameter fits after each wave.

- sim_df:

  the final wave posteriors

- x:

  a `abc_fit` object as output by the `abc_XXX` functions

- ...:

  passed on to methods

- truth:

  a named numeric vector of known parameter values

## Value

an S3 object of class `abc_fit` this contains the following:

- type: the type of ABC algorithm

- iterations: number of completed iterations

- converged: boolean - did the result meet convergence criteria

- waves: a list of dataframes of wave convergence metrics

- summary: a dataframe with the summary of the parameter fits after each
  wave.

- priors: the priors for the fit as a `abc_prior` S3 object

- posteriors: the final wave posteriors

## Methods (by generic)

- `format(abc_fit)`: S3 format method

- `summary(abc_fit)`: S3 summary method

- `print(abc_fit)`: S3 print method

- `plot(abc_fit)`: S3 plot method

## Functions

- `new_abc_fit()`: Create a `abc_fit` object

- `tidy.abc_fit()`: S3 summary method
