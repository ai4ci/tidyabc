# Simulation, scoring and convergence functions

``` r
library(tidyabc)
#> 
#> Attaching package: 'tidyabc'
#> The following objects are masked from 'package:base':
#> 
#>     transform, truncate
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(ggplot2)
```

## Introduction

This vignette provides a detailed guide on writing the core user-defined
functions required for Approximate Bayesian Computation (ABC) using
`tidyabc`: the simulation function (`sim_fn`), the scoring function
(`scorer_fn`), and the convergence function (`converged_fn`).
Understanding the precise input and output structure of these functions
is crucial for correctly implementing an ABC workflow.

------------------------------------------------------------------------

## 1. Writing the Simulation Function (`sim_fn`)

The simulation function (`sim_fn`) defines your generative model. It
takes parameter values as input and returns simulated data.

### Inputs

- **Named Arguments**: The `sim_fn` must accept named arguments
  corresponding to the parameters defined in your `priors_list`. For
  example, if you have
  `priors(norm_mean ~ unif(0, 10),norm_sd ~ unif(0, 10))`, your `sim_fn`
  must accept `norm_mean` and `norm_sd` as arguments.
- **`...` (Optional)**: While not strictly required by the `abc_*`
  functions themselves, it’s good practice to allow for `...` in the
  function definition if you anticipate passing extra arguments during
  debugging or testing.

### Output

- **List**: The function must return a **list** containing the simulated
  data. The structure of this list is entirely up to your model. It
  could be a list of vectors, a list of data frames, or any other R
  object that represents the output of your simulation model.

### Example

``` r
# Example sim_fn matching the parameters from priors defined later
example_sim_fn = function(norm_mean, norm_sd, gamma_mean, gamma_sd) {
  # Simulate data based on parameters
  A = rnorm(1000, norm_mean, norm_sd)
  B = rgamma2(1000, gamma_mean, gamma_sd) # Using a custom dist from tidyabc
  C = rgamma2(1000, gamma_mean, gamma_sd)

  # Return a list containing the simulated outputs
  return(
    list(
      data1 = A + B - C,  # Example composite output
      data2 = B * C       # Example composite output
    )
  )
}
```

------------------------------------------------------------------------

## 2. Writing the Scoring Function (`scorer_fn`)

The scoring function (`scorer_fn`) compares simulated data to observed
data, producing a set of summary statistics (distances or other metrics)
that quantify the “closeness” of the simulation to the observation.

### Inputs

- **`simdata`**: This is the output list from your `sim_fn`.
- **`obsdata`**: This is the observed data provided to the `abc_*`
  function, which should be in the same format as the `simdata` list.

### Output

- **List**: The function must return a **list**. Each element of this
  list corresponds to a *component* of the overall summary statistic
  used for comparison. Commonly, these components are distances (e.g.,
  Wasserstein, RMSE) between corresponding parts of `simdata` and
  `obsdata`, but they could also be other metrics like summary moments
  (mean, variance).
- **Naming**: The names of the list elements returned by `scorer_fn` are
  crucial. They must be consistent and identifiable, as they are used
  for weighting (`scoreweights`) and internal processing.

### Combining Scorer Outputs and Comparison to Observed Scores

1.  **Execution**: The `scorer_fn` is called repeatedly within the ABC
    workflow, once for each simulation (`simdata`) against the fixed
    `obsdata`.
2.  **Component Scores**: The output of `scorer_fn(simdata, obsdata)` is
    a list of *component scores* (e.g.,
    `list(data1 = dist1, data2 = dist2)`).
3.  **Combining Components**: The ABC algorithm combines these component
    scores into a single *summary distance* for each simulation. The
    method is specified by the `distance_method` argument in `abc_*`
    functions:
    - **`"euclidean"` (default)**: Calculates the weighted Euclidean
      distance:  
      \\ d\_{summary} = \sqrt{\sum_i (w_i \cdot s_i)^2} \\ where \\s_i\\
      is the \\i\\-th component score (e.g., `dist1`, `dist2`) and
      \\w_i\\ is the corresponding weight from the `scoreweights`
      vector.
    - **`"manhattan"`**: Calculates the weighted Manhattan distance:  
      \\ d\_{summary} = \sum_i w_i \cdot \|s_i\| \\
    - **`"mahalanobis"`**: Calculates the Mahalanobis distance using the
      covariance matrix of the component scores from the first wave,
      incorporating weights.
4.  **`obsscores`**: This argument allows you to provide pre-calculated
    component scores for the `obsdata` itself (e.g.,
    `list(data1 = 0, data2 = 0)` if comparing against `obsdata` using
    `scorer_fn(obsdata, obsdata)` yields zeros). This is optional and
    often inferred internally if `scorer_fn` is run on `obsdata` itself
    initially.
5.  **Tolerance and Kernels**: The final \\d\_{summary}\\ is compared
    against a tolerance threshold (`epsilon`) using a kernel function
    (defined by the `kernel` argument) to calculate the ABC weights for
    each simulation.

### Example

``` r
# Example scorer_fn matching the simdata structure from example_sim_fn
example_scorer_fn = function(simdata, obsdata) {
  # Compare corresponding elements of simdata and obsdata
  # Using calculate_wasserstein from tidyabc as an example distance
  return(list(
    data1 = calculate_wasserstein(simdata$data1, obsdata$data1),
    data2 = calculate_wasserstein(simdata$data2, obsdata$data2)
  ))
  # Output: list(data1 = <dist_val1>, data2 = <dist_val2>)
}
```

### Score Weights (`scoreweights`)

The `scoreweights` argument allows you to assign relative importance to
each component score returned by `scorer_fn`.

- **Type**: A named numeric vector, where names correspond to the names
  of the list elements returned by `scorer_fn`.
- **Purpose**: If `data1` and `data2` have different scales or
  importance, you can use weights to balance their contribution to the
  summary distance. For instance, `c(data1 = 1, data2 = 2)` means
  `data2`’s distance contributes twice as much as `data1`’s distance to
  the final \\d\_{summary}\\ (before applying the square root for
  Euclidean).
- **Automatic Calculation**: The
  [`posterior_distance_metrics()`](https://ai4ci.github.io/tidyabc/reference/posterior_distance_metrics.md)
  function can analyze results from a trial run (e.g., `abc_rejection`)
  and suggest suitable `scoreweights` based on the scale of the summary
  statistics.

## 3. Debugging and Parallelisation

It is a good idea to run your simulation with some test parameters
first. Similarly testing the scoring function is correctly working is a
good idea. An end to end test harness is supplied in
[`test_simulation()`](https://ai4ci.github.io/tidyabc/reference/test_simulation.md)
that will run your simulation and scorer for a fixed set of parameters.
It can be configured to drop to an interactive debugging browser if it
encounters an error.

Although testing against one set of parameters will catch major problems
edge cases often arise during the main ABC. Each of the ABC functions
can also drop the user to a debugging browser if an issue is detected.
This is not compatible with parallelisation though so best to run a
small set of simulations using
e.g. [`abc_rejection()`](https://ai4ci.github.io/tidyabc/reference/abc_rejection.md).

For parallel execution the simulation and scorer functions must be
completely self contained, and make no reference to global environment
variables, or unqualified references to package functions. Global
references will in general be automatically detected, and you will be
alerted and your functions will be rewritten with common package names
qualified so that they are self contained.

If you need to use globals
[`carrier::crate`](https://rdrr.io/pkg/carrier/man/crate.html)ing your
simulation function will achieve the same thing.

------------------------------------------------------------------------

## 4. The Convergence Function (`converged_fn`)

The convergence function (`converged_fn`) determines when the iterative
ABC-SMC or ABC-Adaptive algorithm should stop. It evaluates the change
in the posterior approximation between waves.

### Inputs to `converged_fn`

The `converged_fn` is called internally by `abc_smc` and `abc_adaptive`
at the end of each wave. It receives three arguments:

- **`wave`**: The wave number as an integer
- **`summary`**: A **single-row data frame** containing summary metrics
  aggregated across *all parameters* for the *current wave*. Common
  columns include:
  - `abs_distance`: The absolute tolerance threshold (`epsilon`) for the
    current wave.
  - `abs_distance_redn`: The absolute reduction in `epsilon` compared to
    the previous wave (`epsilon_{prev} - epsilon_{current}`).
  - `rel_distance_redn`: The relative reduction in `epsilon` compared to
    the previous wave
    (`(epsilon_{prev} - epsilon_{current}) / epsilon_{prev}`).
  - `ESS`: The Effective Sample Size of the current wave’s particle set.
  - Other aggregated metrics might be present depending on the
    algorithm’s implementation.
- **`per_param`**: A **single-row data frame** containing metrics
  calculated *for each individual parameter* in the current wave,
  compared to the previous wave. Common columns include:
  - `param`: The name of the parameter (e.g., “norm_mean”).
  - `IQR_95_redn`: The reduction in the 95% credible interval (IQR)
    range (`q0.975 - q0.025`) for this parameter compared to the
    previous wave (`IQR_{prev} - IQR_{current}`).
  - `abs_variance_redn`: The absolute reduction in the posterior
    variance for this parameter compared to the previous wave
    (`var_{prev} - var_{current}`).
  - `rel_mean_change`: The relative change in the posterior mean for
    this parameter compared to the previous wave
    (`|mean_{current} - mean_{prev}| / |mean_{prev}|`).
  - This data frame has one row for each parameter being estimated.

### Output of `converged_fn`

- **Logical**: The function must return a single `TRUE` or `FALSE`.
  - `TRUE`: Signals that the algorithm has converged, and the iterative
    process stops.
  - `FALSE`: Signals that the algorithm should continue to the next
    wave.

### The `default_termination_fn`

`tidyabc` provides a
[`default_termination_fn()`](https://ai4ci.github.io/tidyabc/reference/default_termination_fn.md)
which implements a common convergence strategy:

- It checks if the `rel_mean_change` for *all* parameters is below a
  specified `stability` threshold.
- *AND* it checks if the `IQR_95_redn` for *all* parameters is below a
  specified `confidence` threshold.
- If both conditions are met, it returns `TRUE`, indicating convergence.

``` r
# Example of using default_termination_fn with custom thresholds
default_converged_fn = default_termination_fn(stability = 0.01, confidence = 0.1)
# This means: converge if mean changes are < 1% and 95% CI shrinks by < 0.1 units for all params.

# If you have 2 parameters `norm_mean` and `norm_sd` you can specify different
# convergence criteria by nameing the stability and confidence inputs:

default_converged_fn2 = default_termination_fn(
  stability = c(norm_mean = 0.01, norm_sd = 0.1),
  confidence = c(norm_mean = 0.01, norm_sd = 0.1)
)

# this would accept more variation in the `norm_sd` estimate compared to that
# of the `norm_mean`
```

------------------------------------------------------------------------

## 5. Writing a Custom Convergence Function

Writing a custom convergence function allows you to tailor the stopping
criteria to your specific needs.

### Structure

A custom convergence function must adhere to the input signature:
`function(wave, summary, per_param)`. It processes the data frames
provided and returns a single `TRUE` or `FALSE`.

### Example: Convergence Based on ESS and Distance Reduction

This example defines convergence as achieving a high Effective Sample
Size (ESS) *and* a minimal reduction in the absolute distance threshold.

``` r
# Define a custom convergence function
custom_converged_fn = function(wave, summary, per_param) {
  # Extract values from the summary data frame (it's a single row)
  current_ess = summary$ESS
  abs_dist_redn = summary$abs_distance_redn # Absolute reduction in epsilon

  # Define criteria
  min_ess_threshold = 500
  max_dist_reduction_threshold = 0.001 # Stop if epsilon isn't reducing much

  # Check conditions
  ess_condition = current_ess >= min_ess_threshold
  dist_condition = abs_dist_redn <= max_dist_reduction_threshold

  # Return TRUE only if BOTH conditions are met
  return(ess_condition && dist_condition)
}

# Example usage (conceptual - requires obsdata, priors, etc. defined)
# smc_result = abc_smc(
#   obsdata = test_obsdata,
#   priors_list = test_priors,
#   sim_fn = example_sim_fn,
#   scorer_fn = example_scorer_fn,
#   n_sims = 1000,
#   acceptance_rate = 0.5,
#   converged_fn = custom_converged_fn # Use the custom function
# )
```

### Example: Convergence Based on Parameter Variance

This example focuses on the stability of parameter uncertainty,
converging when the absolute variance reduction for *all* parameters is
small.

``` r
# Define a custom convergence function based on variance reduction
custom_converged_fn_variance = function(wave, summary, per_param) {
  # Extract the 'abs_variance_redn' column from the per_param data frame
  abs_variance_reductions = per_param$abs_variance_redn

  # Define a threshold for variance reduction
  variance_reduction_threshold = 0.0001

  # Check if ALL absolute variance reductions are below the threshold
  all_stable <- all(abs_variance_reductions < variance_reduction_threshold)

  return(all_stable)
}

# Example usage (conceptual)
# smc_result = abc_smc(
#   obsdata = test_obsdata,
#   priors_list = test_priors,
#   sim_fn = example_sim_fn,
#   scorer_fn = example_scorer_fn,
#   n_sims = 1000,
#   acceptance_rate = 0.5,
#   converged_fn = custom_converged_fn_variance
# )
```

### Key Considerations for Custom Functions

- **Data Frame Access**: Remember that `summary` and `per_param` are
  data frames. Access columns using `$` (e.g., `summary$ESS`,
  `per_param$rel_mean_change`).
- **Aggregation**: The `per_param` data frame has one row per parameter.
  Use functions like [`all()`](https://rdrr.io/r/base/all.html),
  [`any()`](https://rdrr.io/r/base/any.html), or
  [`max()`](https://rdrr.io/r/base/Extremes.html) on the relevant column
  to combine information across parameters for a single logical output.
- **Robustness**: Consider edge cases, like what happens in the very
  first wave if previous values are not available (though the algorithm
  typically handles initialization).
- **Monitoring**: Convergence criteria are crucial for computational
  efficiency and result quality. Choose metrics that reflect the
  stability you require for your specific inference problem.
- **[`browser()`](https://rdrr.io/r/base/browser.html)**: You can get a
  better feel for the statistics available by running an interactive
  browser from within a convergence function while testing. —

## Conclusion

Understanding the precise structure and purpose of `sim_fn`,
`scorer_fn`, and `converged_fn` is important for effectively using
`tidyabc`. The simulation function generates data, the scoring function
quantifies similarity, and the convergence function determines when the
iterative process halts. By carefully specifying these functions, you
can tailor the ABC workflow to your model and inference goals.
