#' Construct a mixture distribution
#'
#' Constructs a mixture distribution from a list of component distributions
#' \eqn{\text{dist}_i} with corresponding weights \eqn{w_i}.
#' The CDF \eqn{F_{\text{mix}}} of the mixture is a weighted sum of the component CDFs:
#' \deqn{
#' F_{\text{mix}}(x) = \sum_{i=1}^{k} w_i \cdot F_i(x)
#' }
#' where \eqn{F_i} is the CDF of the \eqn{i}-th component distribution and
#' \eqn{\sum w_i = 1}. The implementation first evaluates the weighted CDF on a grid
#' of \eqn{x} values (including tail points defined by `tail_p` and potentially
#' knot points from empirical components). The resulting \eqn{(x, F_{\text{mix}}(x))}
#' pairs are then used as input to `empirical_cdf` to create the final smooth or
#' piecewise linear `dist_fns` object representing the mixture distribution.
#'
#' @param dists a `dist_fn_list` of distribution functions
#' @param weights a vector of weights
#' @param steps the number of points that the mixture distribution is
#'   evaluated at to construct the empirical mixture
#' @param tail_p the support fo the tail of the distribution
#' @param name a name for the mixture
#' @inheritDotParams empirical_cdf smooth
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
mixture = function(
  dists,
  weights = 1,
  steps = 200,
  tail_p = 0.0001,
  ...,
  name = "mixture"
) {
  stopifnot(
    is.dist_fns_list(dists)
  )

  .recycle(dists, weights)
  weights = weights / sum(weights)

  steps = steps + (steps + 1) %% 2
  qeval = c(tail_p, 1 - tail_p)

  # knots_list = Filter(Negate(is.null), dists@knots)
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
  xeval = xeval[is.finite(xeval)]
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
  out = empirical_cdf(
    combined$x,
    combined$p,
    name = name,
    ...
  )
  return(out)
}
