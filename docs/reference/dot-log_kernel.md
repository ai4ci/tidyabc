# Generate a log of a kernel function for ABC weight calculation

\$\$ \text{Uniform: } K(d, \epsilon) = \frac{1}{2\epsilon} \mathbb{I}(d
\leq \epsilon) \$\$ \$\$ \text{Triangular: } K(d, \epsilon) =
\frac{2}{\epsilon}\left(1 - \frac{d}{\epsilon}\right) \mathbb{I}(d \leq
\epsilon) \$\$ \$\$ \text{Epanechnikov: } K(d, \epsilon) =
\frac{3}{4\epsilon}\left(1 - \frac{d^2}{\epsilon^2}\right) \mathbb{I}(d
\leq \epsilon) \$\$ \$\$ \text{Biweight: } K(d, \epsilon) =
\frac{15}{16\epsilon}\left(1 - \frac{d^2}{\epsilon^2}\right)^2
\mathbb{I}(d \leq \epsilon) \$\$ \$\$ \text{Gaussian: } K(d, h) =
\frac{1}{\sqrt{2\pi h^2}} \exp\left(-\frac{d^2}{2h^2}\right) \$\$ where
\\d\\ is the distance between observed and simulated summary statistics,
\\\epsilon\\ is the tolerance level, \\h\\ is the bandwidth parameter
for the Gaussian kernel, and \\\mathbb{I}\\ is the indicator function
taking the value 1 if \\d \leq \epsilon\\ and 0 otherwise (this applies
to all kernels except Gaussian).

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
