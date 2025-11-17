# Plot convergence metrics by wave for SMC and adaptive ABC

Plot convergence metrics by wave for SMC and adaptive ABC

## Usage

``` r
plot_convergence(fit)
```

## Arguments

- fit:

  A S3 `abc_fit` object as output by the `abc_XXX` functions

## Value

a patchwork plot of convergence metrics

## Examples

``` r
plot_convergence(
  example_adaptive_fit()
)

```
