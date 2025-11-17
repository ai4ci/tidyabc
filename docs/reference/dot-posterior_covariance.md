# Covariance matrix from the posteriors

Covariance matrix from the posteriors

## Usage

``` r
.posterior_covariance(theta, weights, kernel_t = 1)
```

## Arguments

- theta:

  a set of particles as a matrix

- weights:

  the weights of the previous particles assumed calculated using the
  same kernel

- kernel_t:

  A kernel bandwidth parameter for proposals. This controls the amount
  of noise that particles are perturbed by (and hence the spread of the
  PDF of the proposal distribution), and 1 is approximately 34% of the
  proposal distribution at the centre. Smaller values (default is 0.2)
  give a smaller step size in generating new proposals, and proposals
  will be closer to currently accepted particles.

## Value

a covariance matrix, optionally scaled
