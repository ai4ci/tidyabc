# NEWS

# tidyabc (development version)

* ABC
  + auto calibrate kernel epsilon (? and sample size) on ESS in adaptive fits
* Migrate 
  + scoring utilities from `ggoutbreak`
  + use of empirical CDF in `ggoutbreak`
* Empirical distributions 
  + display splines of empirical fits 
  + update splines when using `empirical` as a link function
  + see: (https://en.wikipedia.org/wiki/Copula_(statistics)#Gaussian_copula)
  + and: (https://en.wikipedia.org/wiki/Sigmoid_function)
* `dist_fns` extensions:
  + convolution
  + skew
  

# tidyabc 0.0.1

* `ai4ci` r-universe release.
* `ABC` estimators with rejection, SMC or adaptive sampling
  + distance functions for Wasserstein and RMSE. 
  + modern support for parallel simulation execution using `furrr`
  + posterior visualisation and resampling
  + convergence statistics and visualisation
* statistic functions wrappers `dist_fns` S3 class
  + support for distributions as dataframe column
  + empirical, transformed and mixture distributions in unified framework
  + distribution plotting
* Empirical distributions
  + from weighted data or from cumulative probabilities
  + link function support for constraints on support
  + monotonic spline based fits in Q-Q (cdf) or logit-Z (data) space with 
    support for `q`, `p`, `d`, `r` functions.
* standalone distribution functions (migrated from `ggrrr`)
  + e.g. re-parametrisations of standard distributions.
  + constrained gamma distribution
  + logit normal distribution
* S3 class for priors
  + formula based constructor
  + print methods
  + search library paths for distributions.
  + constraints
