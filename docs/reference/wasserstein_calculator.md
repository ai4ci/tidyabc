# Generate a function to calculate a wasserstein distance

**\[experimental\]**

## Usage

``` r
wasserstein_calculator(obs, debias = FALSE, bootstraps = 1)
```

## Arguments

- obs:

  A vector of observations

- debias:

  Should the simulations be shifted to match the mean of the observed
  data

- bootstraps:

  Randomly resample from the simulated data points to match the observed
  size this many times and combine the output by averaging. The
  alternative, when this is 1 (the default) matches the sizes by
  selecting and/or repeating the simulated data points in order
  (deterministically)

## Value

a function that takes parameter `sim` and returns a length normalised
wasserstein distance. This is the average distance an individual data
point must be shifted to match the reference data normalised by the
average distance of the reference data from the mean. The function also
takes a `p` parameter which can be a `progressr` progress bar which must
be named.

## Details

This function takes reference data in the forms of individual
observations in the form of for example event times and returns a crated
function to calculate the wasserstein distance from simulated data to
the observed data. For large simulations this will be quicker than
`calculate_wasserstein` but must be used with a crated `scorer_fn`

In the comparison unequal lengths of the data can be accommodated. The
simulated data is recycled, and sampled, until the same length as the
observed data before the comparison.

## Unit tests


    tmp = wasserstein_calculator(0:10)

    # zero if no distance
    testthat::expect_equal(tmp(0:10), 0)
    testthat::expect_equal(tmp(10:0), 0)

## Examples

``` r
# example case counts from an exponential growth process
sim = rexpgrowth(1000, 0.05, 40, 0)
obs = rexpgrowth(1000, 0.075, 40, 0)
obs2 = rexpgrowth(1000, 0.05, 40, 0)

calc = wasserstein_calculator(sim)

# obs is a different distribution to sim (larger growth)
calc(obs)
#> [1] 0.2760168

# obs2 is from the same distribution as sim so the RMSE should be lower:
calc(obs2)
#> [1] 0.1103527
```
