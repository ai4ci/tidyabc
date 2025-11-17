# Generate a non normalised log of kernel for ABC weight calculation

\$\$ \alpha\_{i,j} = \frac{1}{\beta+1}\\ \beta = \alpha +1 \$\$

## Usage

``` r
.log_kernel(
  u,
  epsilon,
  kernel = c(c("epanechnikov", "uniform", "triangular", "biweight", "gaussian"))
)
```

## Arguments

- u:

  the distance (dissimilarilty) between

- epsilon:

  epsilon is a tolerance threshold that controls how closely simulated
  summaries must match the observed ones to be considered plausible.
  This is in the unit of `abc_summary_distance`. Initially the 0.5
  quantile of distances, in subsequent waves this might be decreased. It
  is the scale parameter of the kernel function. \$K_h(\|u\|)\$

- kernel:

  one of "epanechnikov", "uniform", "triangular", "biweight", "gaussian"

## Value

the log of the kernel function
