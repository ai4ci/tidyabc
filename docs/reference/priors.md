# Construct a set of priors

`abc_prior` S3 objects are used to hold the specification of prior and
intermediate proposal distributions. They are inputs to the main
`abc_...()` workflow functions.

## Usage

``` r
priors(...)
```

## Arguments

- ...:

  a list of formulae. Two sided will be interpreted as distribution or
  derived value specifications. One sided as constraints between
  parameters. A distribution is specified as the name of the family of
  statistical distributions and their parameters: e.g.:
  `x ~ norm(mean=3,sd=2)`. The name will be matched to the first hit on
  the search path.

## Value

an S3 object of class `abc_prior` which contains

- a list of `dist_fns`

- a `cor` attribute describing their correlation

- a `derived` attribute describing derive values

- a `constraints` attribute listing the constraints

- a `params` attribute listing the names of the parameters

## Examples

``` r
p = priors(
  mean ~ tidyabc::rgamma2(4,2),
  sd ~ gamma2(2,1),
  shape ~ mean^2 / sd^2,
  rate ~ mean / sd^2,
  ~ mean > sd
)

print(p)
#> Parameters: 
#> * mean: gamma2(mean = 4, sd = 2)
#> * sd: gamma2(mean = 2, sd = 1)
#> Constraints:
#> * mean > sd
#> Derived values:
#> * shape = mean^2/sd^2
#> * rate = mean/sd^2

# Plot methods are also provided:
if (interactive()) plot(p)

# constraints:
p@constraints
#> [[1]]
#> ~mean > sd
#> <environment: 0x566a0bb11538>
#> 
```
