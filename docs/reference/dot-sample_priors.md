# Sample for the prior distribution

This uses a multivariate normal and a copula to generate correlated
structure. The correlation is held as an attribute in the proposals.
N.B. this used to construct derived values and apply constraints but
this was complex and has been deferred. Constraints change the
probability distribution also and will affect validity of proposal
distribution.

## Usage

``` r
.sample_priors(proposal_list, n_sims)
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
      )
    )

    s = .sample_priors(p,10)
    testthat::expect_equal(
      colnames(s),
      c("abc_mvn_mean", "mean", "abc_mvn_sd", "sd")
    )
