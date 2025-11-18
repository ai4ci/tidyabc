#' Generate a distribution from a truncation of another
#'
#' Generates a truncated distribution from an existing distribution `dist`.
#' The new distribution is restricted to the interval \eqn{[x_{left}, x_{right}]}
#' (or \eqn{(x_{left}, x_{right})} if open bounds are intended, though the
#' quantile function will map 0 to \eqn{x_{left}} and 1 to \eqn{x_{right}}).
#' The probability density function (PDF) of the truncated distribution \eqn{f_T(x)}
#' is related to the original PDF \eqn{f(x)} by a normalization constant:
#' \deqn{
#' f_T(x) = \frac{f(x)}{F(x_{right}) - F(x_{left})} \mathbb{I}_{[x_{left}, x_{right}]}(x)
#' }
#' where \eqn{F} is the original CDF and \eqn{\mathbb{I}} is the indicator function.
#' Consequently, the CDF \eqn{F_T} and quantile function \eqn{Q_T} are:
#' \deqn{
#' F_T(x) = \frac{F(x) - F(x_{left})}{F(x_{right}) - F(x_{left})}
#' }
#' \deqn{
#' Q_T(p) = Q(F(x_{left}) + p \cdot (F(x_{right}) - F(x_{left})))
#' }
#' where \eqn{Q} is the quantile function of the original distribution.
#' This function implements these transformations for the `p`, `q`, and `r` functions
#' of the resulting `dist_fns` object.
#'
#' @param dist a distribution as a name, function or a `dist_fn` S3 object
#' @param x_left The lower end of the interval or NA for open
#' @param x_right The upper end of the interval or NA for open
#' @param ... parameters for the underlying distribution if `dist` is a name or
#'  function.
#'
#' @returns a `dist_fn` or `dist_fn_list` holding the truncated distribution(s)
#' @export
#' @concept empirical
#'
#' @unit
#'
#' shape2_gamma = as.dist_fns(pgamma, shape=2)
#' g2 = truncate(shape2_gamma, 0.5, 4)
#'
#' testthat::expect_equal(
#'   format(g2),
#'   "trunc(gamma(shape = 2, rate = 1), 0.50 - 4.00); Median (IQR) 1.68 [1.08 — 2.46]"
#' )
#'
#' testthat::expect_equal(
#'   g2$p(0:5),
#'   c(0, 0.212702667019627, 0.615716430100344, 0.868531239886668, 1, 1)
#' )
#'
#' testthat::expect_equal(g2$q(seq(0, 1, 0.2)), c(
#'   0.5,
#'   0.971743593113941,
#'   1.42711359317907,
#'   1.95304152540128,
#'   2.6642447272259,
#'   4
#' ))
#'
#'
#' g3 = truncate(shape2_gamma, NA, 4)
#' testthat::expect_equal(
#'   format(g3),
#'   "trunc(gamma(shape = 2, rate = 1), 0.00 - 4.00); Median (IQR) 1.54 [0.899 — 2.35]"
#' )
#'
#' g4 = truncate(shape2_gamma, 2, NA)
#' testthat::expect_equal(
#'   format(g4),
#'   "trunc(gamma(shape = 2, rate = 1), 2.00 - Inf); Median (IQR) 2.97 [2.42 — 3.87]"
#' )
#'
truncate = function(dist, x_left, x_right, ...) {
  if (!is.dist_fns_list(dist)) {
    if (!is.dist_fns(dist)) {
      dist = as.dist_fns(dist, ...)
    }
    dist = as.dist_fns_list(dist)
  }

  if (is.na(x_left) && is.na(x_right)) {
    return(dist)
  }

  stopifnot(
    length(x_left) == 1,
    length(x_right) == 1,
    all(x_left < x_right, na.rm = TRUE)
  )

  p_low = ifelse(is.na(x_left), 0, dist$p(x_left))
  p_high = ifelse(is.na(x_right), 1, dist$p(x_right))
  p_delta = p_high - p_low

  x_low = if (is.na(x_left)) dist$q(0) else x_left
  x_high = if (is.na(x_right)) dist$q(1) else x_right
  # x_delta = x_high - x_low

  qfn = NULL
  out = .super_crate(
    x_low = x_low,
    x_high = x_high,
    p_low = p_low,
    p_delta = p_delta,
    dist = dist,
    .fns = list(
      pfn = function(q) {
        base_p = dist$p(q)
        return(
          ifelse(
            q > x_high,
            1,
            ifelse(
              q < x_low,
              0,
              (base_p - p_low) / p_delta
            )
          )
        )
      },
      qfn = function(p) {
        base_p = (p * p_delta) + p_low
        return(dist$q(base_p))
      },
      rfn = function(n) {
        qfn(stats::runif(n))
      }
    )
  )

  name = sprintf("trunc(%s, %1.2f - %1.2f)", dist$name, x_low, x_high)

  return(
    new_dist_fns(
      name = name,
      pfn = out$pfn,
      qfn = out$qfn,
      rfn = out$rfn,
      ...
    )
  )
}
