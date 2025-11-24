# Run the simulation for one set of parameters

Run the simulation for one set of parameters

## Usage

``` r
test_simulation(
  sim_fn,
  scorer_fn,
  ...,
  params = NULL,
  obsdata = NULL,
  seed = NULL,
  debug = FALSE,
  debug_errors = !debug
)
```

## Arguments

- sim_fn:

  a user defined function that takes a set of parameters named the same
  as `priors_list`. It must return a simulated data set in the same
  format as `obsdata`, or that can be compared to `simdata` by
  `scorer_fn`. This function must not refer to global parameters, and
  will be automatically crated with `carrier`.

- scorer_fn:

  a user supplied function that matches the following signature
  `scorer_fn(simdata, obsdata, ....)`, i.e. it takes data in the format
  of `simdata` paired with the original `obsdata` and returns a named
  list of component scores per simulation. This function can make use of
  the `calculate_*()` set of functions to compare components of the
  simulation to the original data. This function must not refer to
  global parameters, and will be automatically crated with `carrier`. If
  this is a purrr style function then `.x` will refer to simulation
  output and `.y` to original observation data.

- ...:

  simulation parameters, must be named

- params:

  a named list of simulation parameters to test (as an alternative to
  including them in `...`)

- obsdata:

  The observational data. The data in this will typically be a named
  list, but could be anything, e.g. dataframe. It is the reference data
  that the simulation model is aiming to replicate.

- seed:

  an optional random seed

- debug:

  start the simulation function in debug mode. This will step through
  both the `sim_fn` and the `scorer_fn` line by line to check that the
  behaviour is as intended.

- debug_errors:

  Errors that crop up in `sim_fn` during a simulation due to anomolous
  value combinations are hard to debug. If this flag is set, whenever a
  `sim_fn` or `scorer_fn` throws an error an interactive debugging
  session is started with the failing parameter combinations. This is
  not compatible with running in parallel.

## Value

a list containing the parameters as `truth` and an instance of the
simulation as `obsdata`, `obsscores` is the result of comparing the
`obsdata` with itself and is usually going to result in zeros. Both
`sim_fn` and `scorer_fn` are rewritten to make sure that all package
names are fully qualified. The rewritten versions are refurned in the
result list (as `sim_fn` and `scorer_fn` respectively)

## Examples

``` r
test = test_simulation(
  example_sim_fn,
  example_scorer_fn,
  # Model parameters to test with:
  mean = 4, sd1 = 3, sd2 = 2,
  obsdata = example_obsdata()
)

# the rewritten function:
cat(test$sim_fn, sep="\n")
#> function (mean, sd1, sd2) 
#> {
#>     return(list(A = stats::rnorm(1000, mean, sd1), B = tidyabc::rgamma2(1000, 
#>         mean, sd2)))
#> }

# The scores resulting from this one simulation, when compared to the
# reference `obsdata`.
test$obsscores
#> $A
#> [1] 0.3578723
#> 
#> $B
#> [1] 0.5690237
#> 
```
