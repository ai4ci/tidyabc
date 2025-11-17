# Plot the evolution of the density function by wave for SMC and adaptive ABC

Plot the evolution of the density function by wave for SMC and adaptive
ABC

## Usage

``` r
plot_evolution(fit, truth = NULL, ...)
```

## Arguments

- fit:

  A S3 `abc_fit` object as output by the `abc_XXX` functions

- truth:

  a named numeric vector of known parameter values

- ...:

  passed on to methods

## Value

a plot of the density functions by wave

## Examples

``` r
plot_evolution(
  example_adaptive_fit(),
  example_truth()
)

```
