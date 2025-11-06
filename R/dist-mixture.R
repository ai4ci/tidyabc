#' Construct a mixture distribution
#'
#' @param dists a `dist_fn_list` of distribution functions
#' @param weights a vector of weights
#'
#' @returns a `dist_fn` of the mixture distribution
#' @export
#' @concept empirical
#'
#' @unit
#'
#' dists = c(
#'   as.dist_fns("norm", mean=-1),
#'   as.dist_fns("norm", mean=1),
#'   as.dist_fns("gamma", shape=2)
#' )
#' weights = c(1,1,0.3)
#'
#' mix = mixture(dists,weights)
#' testthat::expect_equal(
#'   format(mix),
#'   "mixture; Median (IQR) 0.289 [-0.886 — 1.31]"
#' )
#'
mixture = function(dists, weights = 1, steps = 200, tail_p = 0.00001) {
  stopifnot(
    is.dist_fns_list(dists)
  )

  .recycle(dists, weights)
  weights = weights / sum(weights)

  steps = steps + (steps + 1) %% 2
  qeval = c(tail_p, 1 - tail_p)

  knots_list = Filter(Negate(is.null), dists@knots)
  knots_df = dplyr::bind_rows(
    lapply(seq_along(dists), function(i) {
      if (is.null(dists[[i]]@knots)) {
        dplyr::tibble(
          p = qeval,
          x = dists[[i]]$q(qeval),
          index = i
        )
      } else {
        dplyr::mutate(dists[[i]]@knots, index = i)
      }
    })
  )

  xeval = unique(knots_df$x)
  xeval = seq(min(xeval), max(xeval), length.out = steps)

  combined = dplyr::bind_rows(lapply(seq_along(dists), function(i) {
    dplyr::tibble(
      p = dists[[i]]$p(xeval),
      x = xeval,
      wt = weights[[i]]
    )
  }))

  combined = combined %>%
    dplyr::group_by(x) %>%
    dplyr::summarise(p = sum(p * wt)) %>%
    dplyr::arrange(x)
  out = empirical(combined$x, combined$p, tail_p = tail_p, name = "mixture")
  return(out)
}
