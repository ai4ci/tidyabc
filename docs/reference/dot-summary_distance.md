# Combine simulation scores to a single distance metric.

This function calculates the overall score for an individual simulation,
based on combining the components output by `score_list` comparing them
to the target component scores of the original observed data. The
Mahalanobis distance is calculated relative to the distribution of the
first wave (i.e. when the priors were directly sampled).

## Usage

``` r
.summary_distance(
  component_scores,
  obsscores = NULL,
  distance_method = c("euclidean", "manhattan", "mahalanobis"),
  wave1_cov = NULL,
  scoreweights = NULL
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

- scoreweights:

  A named vector with names matching output of `scorer_fn` that defines
  the importance of this component of the scoring in the overall
  distance and weighting of any given simulation. This can be used to
  assign more weight on certain parts of the model output.

## Value

a vector of distances
