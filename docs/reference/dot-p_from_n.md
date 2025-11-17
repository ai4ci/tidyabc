# Proportions from a vector of lengths

Typically useful in finding alpha values for plotting, or bandwidths for
smoothing data this is a decreasing sigmoid with tunable parameters.

## Usage

``` r
.p_from_n(n, n_50, n_01, p_max = 1, n_100 = 0)
```

## Arguments

- n:

  a vector of lengths

- n_50:

  the size at 50 percent

- n_01:

  the size at 1 percent, and the limit of the returned `p` value such
  that \\p \times n \ge n\_{01}\\

- p_max:

  the maximum value at `n <= n_100`

- n_100:

  an offset under which the value is 1 (default is 0)

## Value

a vector of proportions for each length

## Details

This is using a composite of a [Hill
function](https://en.wikipedia.org/wiki/Hill_equation_(biochemistry))
and a reciprocal:

\$\$ f(n) = \frac{1}{1 + \left( \dfrac{n - n\_{100}}{n\_{50}} \right)^b}
\times p\_{max}\\ b = \frac{\ln(99)}{\ln(n\_{01}/n\_{50})} \$\$

If \\n-n\_{100} \le 0\\ then the maximum value is returned

## Unit tests



    testthat::expect_equal(.p_from_n(c(0, 5, 100), 5, 100), c(1, 0.5, 0.01))

    testthat::expect_error(
      {
        .p_from_n(0, 50, 1)
      },
      "inconsistent parameters, n_01 must be larger than n_50.",
      fixed = TRUE
    )

    testthat::expect_equal(.p_from_n(c(0:10), 5, 10, 0.1), c(
      0.1,
      0.0999976758268449,
      0.0997704297474442,
      0.0967278219151741,
      0.081446655140766,
      0.05,
      0.0229935645768342,
      0.0097036542590398,
      0.00424593238101868,
      0.00199056073397113,
      0.001
    ))

    # does not go below
    testthat::expect_equal(
      .p_from_n(c(1000, 10000, 100000), 5, 50),
      c(5e-04, 5e-05, 5e-06)
    )
