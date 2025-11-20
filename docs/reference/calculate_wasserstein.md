# Calculate a Wasserstein distance

This function takes simulation and observed data and calculates a
normalised Wasserstein distance.

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
simulated data is sorted and linearly interpolated to the same length as
the observed data before the comparison.

\$\$ W = \frac{1}{N\_{obs} \cdot \sigma\_{obs}} \sum\_{i=1}^{N\_{obs}}
\| \hat{x}\_i - y_i \| \$\$

where \\y_i\\ are the ordered observed data points, \\\hat{x}\_i\\ are
the simulated data points after matching size (potentially via
interpolation), debiasing (optional), and sorting, \\N\_{obs}\\ is the
number of observed points, and \\\sigma\_{obs}\\ is the standard
deviation of the observed data (used for normalisation).

Size matching via linear interpolation:

Let \\x\\ be the sorted simulated data of length \\N\_{sim}\\, and
\\N\_{obs}\\ be the target length.

Indices \\idx\\ are generated as \\idx = \frac{(i-1) \cdot (N\_{sim} -
1)}{N\_{obs} - 1} + 1\\ for \\i = 1, ..., N\_{obs}\\. Then
\\\hat{x}\_i\\ is calculated as:

\$\$ \hat{x}\_i = (1 - p_i) \cdot x\_{\lfloor idx_i \rfloor} + p_i \cdot
x\_{\lceil idx_i \rceil} \$\$

where \\p_i = idx_i - \lfloor idx_i \rfloor\\. If \\\lfloor idx_i
\rfloor = \lceil idx_i \rceil\\, then \\\hat{x}\_i = x\_{\lfloor idx_i
\rfloor}\\.

For bootstrapping (`bootstraps > 1`), the process is repeated
`bootstraps` times with random sampling, and the final distance is the
average of the results.

## Unit tests



    testthat::expect_equal(calculate_wasserstein(0:10, 10:0), 0)

    # zero if no distance
    testthat::expect_equal(calculate_wasserstein(0:10, 0:10), 0)
    testthat::expect_equal(calculate_wasserstein(10:0, 0:10), 0)

    # normalised so that all mass at mean = 1
    ref = mean(abs(0:10-5))/sd(0:10)
    testthat::expect_equal(calculate_wasserstein(rep(5, 11), 0:10), ref)

    # smaller sample recycled and normalises to same value
    testthat::expect_equal(calculate_wasserstein(rep(5, 5), 0:10), ref)

    # should be ((0+1+0+1+0+0+0+1+0+1+0) / 11) / ((5+4+3+2+1+0+1+2+3+4+5) / 11) = 0.1333...
    testthat::expect_equal(
      calculate_wasserstein(c(0, 0, 2, 2, 4, 5, 6, 8, 8, 10, 10), 0:10),
      0.109640488937369
    )

    withr::with_seed(100, {
       ref = rnorm(1000)
       cmp1 = rnorm(1000)
       cmp2 = rnorm(1000, sd=2)
       cmp3 = rnorm(100)
    })

    testthat::expect_equal(
      calculate_wasserstein(cmp1, ref),
      0.0458673541615467
    )
    testthat::expect_equal(
      calculate_wasserstein(cmp2, ref),
      0.835125114099775
    )
    testthat::expect_equal(
      calculate_wasserstein(cmp3, ref),
      0.180171816429487
    )

    tmp = withr::with_seed(100, {
     calculate_wasserstein(cmp1,ref,bootstraps=10)
    })
    testthat::expect_equal(tmp, 0.0678606864134826)


    gen = function(n, mean=0, sd=1) {
      sample(pnorm(seq(-1,1,length.out = n - n
    }

    # there should be approximately zero
    testthat::expect_equal(calculate_wasserstein(gen(1000), gen(1000)), 0)
    testthat::expect_equal(calculate_wasserstein(gen(1000), gen(100)), 0)
    testthat::expect_equal(
      calculate_wasserstein(gen(100), gen(1000)),
      0,tolerance = 0.01
    )
    testthat::expect_equal(
      calculate_wasserstein(gen(200), gen(1000)),
      0,tolerance = 0.01
    )

    # these should be approximately equal:
    tmp2 = max(abs(diff(c(
    calculate_wasserstein(gen(100,0.1), gen(1000)),
    calculate_wasserstein(gen(200,0.1), gen(1000)),
    calculate_wasserstein(gen(1000,0.1), gen(1000)),
    calculate_wasserstein(gen(1000, 0.1), gen(200)),
    calculate_wasserstein(gen(1000, 0.1), gen(100))
    ))))

    testthat::expect_equal(tmp2, 0, tolerance = 0.01)
