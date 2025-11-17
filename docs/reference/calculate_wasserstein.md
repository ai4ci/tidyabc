# Calculate a wasserstein distance

This function takes simulation and observed data and calculates a
normalised wasserstein distance.

## Usage

``` r
calculate_wasserstein(sim, obs, debias = FALSE, bootstraps = 1)
```

## Arguments

- sim:

  A vector of simulated data points

- obs:

  A vector of observed data points

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

a length normalised wasserstein distance. This is the average distance
an individual simulated data point must be shifted to match the observed
data normalised by the average distance of the observed data from the
mean.

## Details

In the comparison unequal lengths of the data can be accommodated. The
simulated data is recycled, and sampled, until the same length as the
observed data before the comparison.

## Unit tests



    testthat::expect_equal(calculate_wasserstein(0:10, 10:0), 0)

    # zero if no distance
    testthat::expect_equal(calculate_wasserstein(0:10, 0:10), 0)
    testthat::expect_equal(calculate_wasserstein(10:0, 0:10), 0)

    # normalised so that all mass at mean = 1
    testthat::expect_equal(calculate_wasserstein(rep(5, 11), 0:10), 1)

    # smaller sample recycled and normalises to same value
    testthat::expect_equal(calculate_wasserstein(rep(5, 5), 0:10), 1)

    # should be ((0+1+0+1+0+0+0+1+0+1+0) / 11) / ((5+4+3+2+1+0+1+2+3+4+5) / 11) = 0.1333...
    testthat::expect_equal(
      calculate_wasserstein(
        c(0, 0, 2, 2, 4, 5, 6, 8, 8, 10, 10),
        0:10
     ),
      0.133333333333333
    )

    withr::with_seed(100, {
       ref = rnorm(1000)
       cmp1 = rnorm(1000)
       cmp2 = rnorm(1000, sd=2)
       cmp3 = rnorm(100)
    })

    testthat::expect_equal(calculate_wasserstein(cmp1,ref), 0.0576417503974498)
    testthat::expect_equal(calculate_wasserstein(cmp2,ref), 1.04950621760385)
    testthat::expect_equal(calculate_wasserstein(cmp3,ref), 0.212977231452674)

    calculate_wasserstein(cmp1,ref,bootstraps=10)
