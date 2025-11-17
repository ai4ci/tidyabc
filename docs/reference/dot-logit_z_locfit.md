# Weighted distribution function interpolation in a logit z space

Weighted cumulative probabilities are mapped to a logit space, data is
transformed to Z space (assumes support is -Inf..Inf). This does not use
a link function and the resulting interpolation functions are not
vectorised. Importance weighting is done during CDF construction.
Prediction is done using a weighted linear interpolation of nearby
points. Weighting for interpolation is a distance based gaussian kernel
from data points to interpolation point. OOB interpolation is supported.

## Usage

``` r
.logit_z_locfit(x, w = NULL, bw = NULL)
```

## Arguments

- x:

  either a vector of samples from a distribution `X` or cut-offs for
  cumulative probabilities when combined with `p`

- w:

  for data fits, a vector the same length as `x` giving the importance
  weight of each observation. This does not need to be normalised. There
  must be some non zero weights, and all must be finite.

- bw:

  a bandwidth expressed in terms of the probability width, or proportion
  of observations.

## Value

a function that will predict a quantile assuming infinite support
