# Changelog

## tidyabc (development version)

- S3 class for priors
  - formula based constructor
  - print methods
  - search library paths for distributions.
  - constraints
- Migrate scoring utilities from `ggoutbreak`
  - migrate use of empirical CDF in `ggoutbreak`
- Convert empirical to use 0..1 range rather than `-Inf..Inf`
  - consider a new quantile link class that maps `0..1` (copula? / PIT?)
    or convert existing
  - store splines in accessible form.
  - update splines when using `empirical` as a link function
  - see:
    (<https://en.wikipedia.org/wiki/Copula_(statistics)#Gaussian_copula>)
  - and: (<https://en.wikipedia.org/wiki/Sigmoid_function>)
  - ?return environment rather than list with splines built in.
- `dist_fns` extensions:
  - convert function wrappers to store call rather than function.
  - convolution
  - skew

## tidyabc 0.0.1

- `ai4ci` r-universe release.
- `ABC` estimators with rejection, SMC or adaptive sampling
  - distance functions for Wasserstein and RMSE.
  - modern support for parallel simulation execution using `furrr`
  - posterior visualisation and resampling
  - convergence statistics and visualisation
- statistic functions wrappers `dist_fns` S3 class
  - support for use as dataframe column
  - empirical and mixture distributions
- standalone distribution functions (migrated from `ggrrr`)
  - e.g. re-parametrisations of standard distributions.
