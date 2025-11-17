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

    # normalised so that all mass at mean = 1
    testthat::expect_equal(tmp(rep(5, 11)), 1)

    # smaller sample recycled and normalises to same value
    testthat::expect_equal(tmp(rep(5, 5)), 1)

    # should be ((0+1+0+1+0+0+0+1+0+1+0) / 11) / ((5+4+3+2+1+0+1+2+3+4+5) / 11) = 0.1333...
    testthat::expect_equal(
      tmp(c(0, 0, 2, 2, 4, 5, 6, 8, 8, 10, 10)),
      0.133333333333333
    )

    withr::with_seed(100, {
       ref = rnorm(1000)
       cmp1 = rnorm(1000)
       cmp2 = rnorm(1000, sd=2)
       cmp3 = rnorm(100)
    })

    tmp2 = wasserstein_calculator(ref)

    testthat::expect_equal(tmp2(cmp1), 0.0576417503974498)
    testthat::expect_equal(tmp2(cmp2), 1.04950621760385)
    testthat::expect_equal(tmp2(cmp3), 0.212977231452674)
