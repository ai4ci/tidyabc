# Calculate the probability of new proposals based on previous

Calculate the probability of new proposals based on previous

## Usage

``` r
.log_q_proposal_vectorized(theta_new, theta_prev, w_prev, Sigma_pert)
```

## Arguments

- theta_new:

  a set of new particles as a matrix

- theta_prev:

  a set of previous particles as a matrix

- w_prev:

  the weights of the previous particles assumed calculated using the
  same kernel

- Sigma_pert:

  the covariance of the kernel.

## Value

a vector of probabilities for each row of theta_new

## Unit tests


    # Previous particles: 100 samples from N(0, 1)
    N_prev <- 100
    theta_prev <- matrix(rnorm(N_prev, 0, 1), ncol = 1)
    w_prev <- rep(1/N_prev, N_prev)  # uniform weights
    Sigma_pert <- matrix(0.5^2, nrow = 1, ncol = 1)  # scalar variance

    # New point
    theta_new <- matrix(c(0.5), ncol = 1)

    # Compute via vectorized function
    q_val_fast <- .log_q_proposal_vectorized(theta_new, theta_prev, w_prev, Sigma_pert)

    # Compute via slow loop (ground truth)
    q_val_slow <- 0
    for (j in 1:N_prev) {
      q_val_slow <- q_val_slow + w_prev[j] * dnorm(theta_new, mean = theta_prev[j], sd = sqrt(Sigma_pert))
    }
    expect_equal(exp(q_val_fast), as.vector(q_val_slow))
