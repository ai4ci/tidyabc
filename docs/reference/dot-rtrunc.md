# Sample from an interval of a distribution

Sample from an interval of a distribution

## Usage

``` r
.rtrunc(n, x_left, x_right, distr, ..., simplify = FALSE)
```

## Arguments

- n:

  The number of samples or a vector the same length as the number of
  samples

- x_left:

  The lower end of the interval or NA for open

- x_right:

  The upper end of the interval or NA for open

- distr:

  a name e.g. (`"norm"`) or statistical function (e.g. `rnorm`). `p` and
  `q` equivalents must exist

- ...:

  parameters for the underlying distribution.

- simplify:

  instead of a list return a matrix with `n` columns and
  `length(x_left)` rows, or a plain vector of length `n`.

## Value

a list of samples the same length as `x_left` and `x_right`

## Unit tests



    testthat::expect_equal(
      withr::with_seed(123, .rtrunc(5, 0, NA, rnorm, simplify = TRUE)),
      c(
        0.36860458186587,
        1.24891855422508,
        0.537354048150005,
        1.56756538386253,
        1.88423854662693
      )
    )

    testthat::expect_equal(
      withr::with_seed(
        123,
        .rtrunc(5, 100, 101, "norm", mean = 0, sd = 2, simplify = TRUE)
      ),
      c(
        100.287577520125,
        100.788305135444,
        100.408976921812,
        100.883017404005,
        100.940467284294
      )
    )

    testthat::expect_equal(
      withr::with_seed(
        123,
        .rtrunc(5, 0, 0, "norm", mean = 0, sd = 2, simplify = TRUE)
      ),
      c(0, 0, 0, 0, 0)
    )

    testthat::expect_equal(
      withr::with_seed(123, .rtrunc(5, 0.2, 0.3, stats::qexp, r = 2, simplify = TRUE)),
      c(
        0.226768410134133,
        0.277097702724268,
        0.238513388505999,
        0.28721473064001,
        0.29345270111699
      )
    )

    # Some part of range is possible => works and only produces positives.
    testthat::expect_equal(
      withr::with_seed(
        123,
        .rtrunc(5, -1, 1, "gamma", shape = 1, simplify = TRUE)
      ),
      c(
        0.200628506352649,
        0.689760686697004,
        0.29911075895734,
        0.81683790576972,
        0.902606552073028
      )
    )

    # No part of range is possible => only produces NAs.
    testthat::expect_equal(
      withr::with_seed(
        123,
        .rtrunc(5, -2, -1, "gamma", shape = 1, simplify = TRUE)
      ),
      c(NA, NA, NA, NA, NA)
    )

    testthat::expect_error(
      {
        .rtrunc(5, 2, 1, "gamma", shape = 1, simplify = TRUE)
      },
      "all `x_left` must be smaller than `x_right`",
      fixed = TRUE
    )
