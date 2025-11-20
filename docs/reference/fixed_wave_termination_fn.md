# Run the SMC or adaptive algorithm for a set number of waves

Run the SMC or adaptive algorithm for a set number of waves

## Usage

``` r
fixed_wave_termination_fn(max_wave)
```

## Arguments

- max_wave:

  the number of waves to run

## Value

A function that will test for completion

## Examples

``` r
# Declare converged after 3 waves:
converged_fn = fixed_wave_termination_fn(3)
```
