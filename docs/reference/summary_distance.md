# Combine simulation scores by euclidean distance.

calculate the overall score for an individual simulation, based on the
components output by `score_list`. This function is applied to a list of
simulation scores and combines them into a single metric by comparing
them to the scores of the original observed data

## Usage

``` r
summary_distance(
  component_scores,
  obsscores = NULL,
  distance_method = c("euclidean", "manhattan", "mahalanobis"),
  wave1_cov = NULL
)
```

## Arguments

- component_scores:

  - a list column of scores. each entry is itself a list

- obsscores:

  Summary scores for the observational data. This will be a named list,
  and is equivalent to the output of `scorer_fn`, on the observed data.
  If not given typically it will be assumed to be all zeros.

- distance_method:

  what metric is used to combine `simscores` and `obsscores` and is one
  of `"euclidean"`, `"manhattan"`, or `"mahalanobis"`.

## Value

a vector of distances (with a `cov` attribute of the wave 1 covariance
matrix)

## Details

TODO: needs refactoring to work with whole sets of scores rather than
individual ones, so we can do a mahalabois distance.
