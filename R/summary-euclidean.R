#' Combine simulation scores by euclidean distance.
#'
#'
#' calculate the overall score for an individual
#' simulation, based on the components output by `score_list`. This function
#' is applied to a list of simulation scores and combines them into a
#' single metric by comparing them to the scores of the original observed data
#'
#' TODO: needs refactoring to work with whole sets of scores rather than
#' individual ones, so we can do a mahalabois distance.
#'
#' @param component_scores - a list column of scores. each entry is itself a list
#' @inheritParams tidyabc_common
#'
#' @returns a vector of distances (with a `cov` attribute of the wave 1
#'   covariance matrix)
#'
#' @keywords internal
summary_distance = function(
  component_scores,
  obsscores = NULL,
  distance_method = c("euclidean", "manhattan", "mahalanobis"),
  wave1_cov = NULL
) {
  method = match.arg(distance_method)

  # unnest scores from list of lists into a matrix:
  # this relies on the naming order being consistent, which it should be as output
  # from purrr...
  simscores = matrix(
    unname(unlist(component_scores)),
    nrow = length(component_scores),
    byrow = TRUE
  )
  # simscores = as.matrix(as.data.frame(Reduce(
  #   f = function(init, x2) {
  #     setNames(
  #       lapply(names(init), function(nm) c(init[[nm]], x2[[nm]])),
  #       names(init)
  #     )
  #   },
  #   x = component_scores
  # )))

  if (is.null(obsscores)) {
    obsscores = rep(0, ncol(simscores))
  }
  obsscores = matrix(
    unlist(obsscores),
    nrow = nrow(simscores),
    ncol = length(obsscores),
    byrow = TRUE
  )

  if (is.null(wave1_cov)) {
    wave1_cov = cov(simscores)
  }

  # if (ncol(simscores) != ncol(obsscores)) {
  #   browser()
  #
  # }

  if (method == "euclidean") {
    delta = simscores - obsscores
    return(
      structure(
        sqrt(rowSums(delta^2)),
        cov = wave1_cov
      )
    )
  }

  if (method == "manhattan") {
    delta = unlist(simscores[names(obsscores)]) - unlist(obsscores)
    return(
      structure(
        rowSums(abs(delta)),
        cov = wave1_cov
      )
    )
  }

  if (method == "mahalanobis") {
    return(
      structure(
        stats::mahalanobis(simscores, obsscores, wave1_cov),
        cov = wave1_cov
      )
    )
  }
}

# x = list(list(a=1,b=2,c=3),list(a=1,b=2,c=3),list(a=1,b=2,c=3))
