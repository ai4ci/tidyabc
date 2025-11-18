# Set up default convergence criteria for SMC and adaptive ABC

Convergence is assessed on firstly whether the central estimate of the
parameters being assessed is stable, and not changing from one wave to
the next, and secondly if the 95 percent credible interval is stable
between waves. If the parameter central estimate is stationary but the
credible intervals are still dropping then continuing simulation may get
better estimates of confidence.

## Usage

``` r
default_termination_fn(stability = 0.01, confidence = 0.1)
```

## Arguments

- stability:

  how close do sequential estimates need to be before declaring them as
  a good set of parameter estimates. This is in the units of the
  parameter. If this is given as a single number it applies to all
  parameters equally, alternatively a named vector can be used to set
  parameter specific cutoffs.

- confidence:

  how stable do the 95% confidence intervals need to be before we are
  happy with the parameter estimates. This is in the scale of the
  parameters, but represents a change in IQR from wave to wave. If this
  is given as a single number it applies to all parameters equally,
  alternatively a named vector can be used to set parameter specific
  cutoffs.

## Value

a function that specifies the convergence.

## Examples

``` r
# A more permissive definition of convergence has
# less strict stability criteria (sequential estimates varying by less than 5%)
# and confidence intervals not changing by more than 1 unit between waves.)
check = default_termination_fn(0.05, 1.0)


fit = example_smc_fit()
#> ABC-SMC
#> SMC waves:  ■■■■■■■■■                         25% | wave 2 ETA:  4s
#> SMC waves:  ■■■■■■■■■■■■■■■■■                 52% | wave 4 ETA:  2s
#> Converged on wave: 8
#> SMC waves:  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      89% | wave 7 ETA:  1s

# This is performed as an integral part of the SMC and adaptive
# fitting and is here only for example
last_wave_metrics = utils::tail(fit$waves,1)
converged = check(
  last_wave_metrics$summary[[1]],
  last_wave_metrics$per_param[[1]]
)

if (isTRUE(converged)) print("Converged (permissive definition)")
#> [1] "Converged (permissive definition)"
```
