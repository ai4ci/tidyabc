#' Generate a distribution from a truncation of another
#'
#' @param dist distribution(s) as a name, function or a `dist_fn` S3 object
#' @param link a link function (or name of a link function)
#' @param ... parameters for the underlying distribution if `dist` is a name or
#'  function.
#' @param name a name for the link function
#'
#' @returns a `dist_fn` or `dist_fn_list` holding the transformed distribution(s)
#' @export
#' @concept empirical
#'
#' @unit
#'
#'
#' t = transform("log","norm")
#' ps = seq(0,1,0.1)
#' qs = 0:6
#' testthat::expect_equal(t$q(ps),qlnorm(ps))
#' testthat::expect_equal(t$p(qs),plnorm(qs))
#' testthat::expect_equal(t$d(qs),dlnorm(qs))
#'
#' t2 = transform("log","norm", 0.4, 0.1)
#' testthat::expect_equal(t2$q(ps),qlnorm(ps,0.4,0.1))
transform = function(link, dist, ..., name = NULL) {
  if (!is.dist_fns_list(dist)) {
    if (!is.dist_fns(dist)) {
      dist = as.dist_fns(dist, ...)
    }
    dist = as.dist_fns_list(dist)
  }

  link = as.link_fns(link)

  if (is.null(name)) {
    name = link$name
  }

  pfn = carrier::crate(
    function(q) {
      q2 = link$trans(q)
      p2 = dist$p(q2)
      return(p2)
    },
    dist = dist,
    link = link
  )
  qfn = carrier::crate(
    function(p) {
      q2 = dist$q(p)
      q = link$inv(q2)
      return(q)
    },
    dist = dist,
    link = link
  )
  rfn = function(n) {
    qfn(stats::runif(n))
  }
  dfn = carrier::crate(
    function(x) {
      tmp1 = dist$d(link$trans(x))
      tmp2 = link$ddxtrans(x)
      ifelse(
        (tmp1 == 0 & !is.finite(tmp2) | !is.finite(tmp1) & tmp2 == 0),
        0,
        tmp1 * tmp2
      )
    },
    dist = dist,
    link = link
  )
  name = sprintf("%s(%s)", name, dist$name)

  return(
    new_dist_fns(
      name = name,
      pfn = pfn,
      qfn = qfn,
      rfn = rfn,
      dfn = dfn
    )
  )
}
