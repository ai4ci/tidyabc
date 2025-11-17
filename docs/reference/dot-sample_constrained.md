# Apply derived values and constraints to samples

Repeatedly samples and removes invalid until the desired number of valid
samples is reached. This also calculates derived values

## Usage

``` r
.sample_constrained(proposal_list, n_sims, sampler_fn = .sample_priors, ...)
```

## Arguments

- proposal_list:

  a list of empirical probability distributions that map MVN space to
  proposal space, and are the "prior" for each adaptive wave. This is
  already used to generate the proposals and their mapping in `sim_df`

- n_sims:

  The number of simulations to run per wave (for SMC and Adaptive) or
  overall (for Rejection). For rejection sampling a large number is
  recommended, for the others sma

- sampler_fn:

  a function that creates random samples

- ...:

  passed on to `sampler_fn`

## Value

a data frame of samples in MVN (prefixed `abc_mvn_`) and parameter
space.

## Unit tests


    p = new_abc_prior(
      .dists = list(
        mean = as.dist_fns("norm",4,2),
        sd = as.dist_fns("gamma",2)
      ),
      .derived = list(
        shape ~ mean^2 / sd^2,
        rate ~ mean / sd^2
      ),
      .constraints = list(
        ~ mean > sd
      )
    )

    s = .sample_constrained(p,1000)

    testthat::expect_equal(
      colnames(s),
      c("abc_mvn_mean", "mean", "abc_mvn_sd", "sd", "shape", "rate")
    )

    testthat::expect_equal(all(s$mean > s$sd), TRUE)
