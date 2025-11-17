# Unified interface ----

#' Fit a piecewise logit transformed linear model to cumulative data
#'
#' This creates statistical distribution functions fitting to data in a
#' transformed space. Inputs can either be weighted data observations (`x`,`w` pairs)
#' of a CDF (`x`,`p` pairs). X axis transformation is specified in the `link`
#' parameter and is either something like "log", "logit", etc or can also be
#' specified as a statistical distribution, or even a length 2 numeric vector
#' defining support.
#'
#' If the empirical distribution if from data it is fully transformed to a
#' `logit-z` space where the cumulative weighted data is interpolated with a
#' kernel weighted linear model.
#'
#' In the case of CDF data the empirical distribution is fitted with a piecewise
#' linear or monotonically increasing spline fit to data in `Q-Q` space. The
#' spline is bounded by 0 and 1 in both dimensions.
#'
#' This function imputes tails of distributions within the constraints of the
#' link functions. Given perfect data as samples or as quantiles it should
#' well approximate the tail.
#'
#' @param x either a vector of samples from a distribution `X` or cut-offs for
#'   cumulative probabilities when combined with `p`
#' @param p for CDF fits, a vector the same length as `x` giving
#'   `P(X <= x)`.
#' @param w for data fits, a vector the same length as `x` giving the
#'   importance weight of each observation. This does not need to be normalised.
#'   There must be some non zero weights, and all must be finite.
#' @param link a link function. Either as name, or a `link_fns` S3 object. In
#'   the latter case this could be derived from a statistical distribution by
#'   `as.link_fns(<dist_fns>)`. This supports the use of a prior to define the
#'   support of the empirical function, and is designed to prevent tail
#'   truncation. Support for the updated quantile function will be the same as
#'   the provided prior.
#' @param fit_spline for distributions from data should we fit a spline to
#'   reduce memory usage and speed up sampling? If `knots` is given this is
#'   assumed to be true, set this to TRUE to allow knots to be determined by the
#'   data.
#' @param knots for spline fitting from data how many points do we use to model
#'   the cdf? By default this will be estimated from the data. I recommend an
#'   uneven number, without a lot of data this will tend to overfit 9 knots for
#'   1000 samples seems OK, max 7 for 250, 5 for 100. Less is usually more.
#' @param name a label for this distribution. If not given one will be generated
#'   from the distribution and parameters. This can be used as part of plotting.
#' @inheritDotParams .logit_z_interpolation
#' @inheritDotParams empirical_cdf
#'
#'
#' @returns a `dist_fns` S3 object containing functions `p()` for CDF, `q()`
#'   for quantile, and `r()` for a RNG, and `d()` for density. The density
#'   function may be discontinuous.
#' @concept empirical
#' @export
#'
#' @unit
#'
#' # from samples:
#' withr::with_seed(123,{
#'  e2 = empirical(rnorm(10000))
#'  testthat::expect_equal(e2$p(-5:5), pnorm(-5:5), tolerance=0.01)
#'  testthat::expect_equal(e2$q(seq(0,1,0.1)), qnorm(seq(0,1,0.1)), tolerance=0.05)
#' })
#'
#' p2 = seq(0,1,0.1)
#' testthat::expect_equal( e2$p(e2$q(p2)), p2, tolerance = 0.001)
#'
#' # updating a prior, with a horribly skewed gamma distribution
#' # not a great fit but not great data
#' withr::with_seed(124,{
#'   data = rgamma(500,1)
#'   e4 = empirical(data, link=as.dist_fns("unif",0,10))
#'   if (interactive()) plot(e4)+ggplot2::geom_function(fun = ~ dgamma(.x, 1))
#'   testthat::expect_equal(mean(e4$r(10000)), 1, tolerance=0.1)
#'
#'   e5 = empirical(data,link="log")
#'   testthat::expect_equal(mean(e5$r(10000)), 1, tolerance=0.1)
#'   if (interactive()) plot(e5)+ggplot2::geom_function(fun = ~ dgamma(.x, 1))
#' })
#'
#' withr::with_seed(123,{
#'   data = c(rnorm(200,4),rnorm(200,7))
#'   weights = c(rep(0.1,200), rep(0.3,200))
#'   e6 = empirical(x=data,w = weights)
#'   plot(e6)
#'   testthat::expect_equal(
#'     format(e6),
#'     "empirical; Median (IQR) 6.56 [5.23 — 7.42]"
#'   )
#' })
#'
#' # Construct a normal using a sequence and density as weight.
#' e7 = empirical(x=seq(-10,10,length.out=1000),w=dnorm(seq(-10,10,length.out=1000)),knots = 20)
#' testthat::expect_equal(e7$p(-5:5), pnorm(-5:5), tolerance=0.01)
#'
empirical = function(
  x,
  ...,
  w = NULL,
  p = NULL,
  link = "ident",
  fit_spline = !is.null(knots),
  knots = NULL,
  name = NULL
) {
  if (!is.null(p) && !is.null(w)) {
    stop("Cannot provide both `p` and `w` parameters.")
  }

  if (!is.null(p)) {
    if (length(x) != length(p) || length(unique(x)) < 2) {
      stop(
        "`x` and `p` must be the same length and there must be at least 2 different values for `x`"
      )
    }
    # empirical_cdf
    return(
      empirical_cdf(
        x = x,
        p = p,
        link = link,
        name = name,
        ...
      )
    )
  } else {
    if (!is.null(w) && length(x) != length(w)) {
      stop("`w` must be the same length as `x` if given.")
    }
    # empirical_data
    fit = empirical_data(
      x = x,
      w = w,
      link = link,
      name = name,
      ...
    )

    if (isTRUE(fit_spline)) {
      # empirical_cdf from empirical data

      # must be fewer knots than data, and must be odd number.
      if (is.null(knots)) {
        knots = .knots_from_data(x)
      }
      knots = min(c(knots, length(x) - 1))
      knots = knots %/% 2 * 2 + 1

      tail_p = 0.0005
      knotp = seq(tail_p, 1 - tail_p, length.out = knots)
      knotx = fit$q(knotp)

      return(empirical_cdf(
        x = knotx,
        p = knotp,
        link = link,
        name = name,
        ...
      ))
    } else {
      return(fit)
    }
  }
}

.knots_from_data = function(x) {
  round(log(length(x)) / 2) * 2 + 1
}


### Old fitting from data logic
#   # fitting from data
#
#   valid = weights > 0 & !is.na(x) & !is.na(weights)
#   x = x[valid]
#   weights = weights[valid]
#   weights = weights / sum(weights)
#
#   tail_p = 0.005
#
#   if (method == "sgolay") {
#     knotp = seq(tail_p, 1 - tail_p, length.out = knots)
#     # knotp = seq(tail_p, 1 - tail_p, length.out = knots * 10 +1)
#     knotx = wquantile(
#       p = knotp,
#       x = x,
#       w = weights,
#       link = link,
#       bw = 2 / (2 + knots)
#     )
#     # # Smoothing in quantiles.
#     # if (smooth) {
#     #   knotx = signal::sgolayfilt(knotx, n = 5, p = 1)
#     # }
#     # knotp = knotp[seq_along(knotp) %% 10 == 1]
#     # knotx = knotx[seq_along(knotx) %% 10 == 1]
#   } else {
#     # locfit approach:
#     rangex = wquantile(
#       p = c(tail_p, 1 - tail_p),
#       x = x,
#       w = weights,
#       link = link,
#       bw = 1 / knots
#     )
#
#     # try to stop a newsplit: out of vertex space error
#     minx = rangex[1]
#     maxx = rangex[2]
#     knotx = seq(minx, maxx, length.out = knots)
#
#     nn = 1 / knots * 1.25
#     # h = (maxx - minx) / (knots / 2) * 1.25
#
#     # manually calculate weights
#     ord = order(x)
#     x = x[ord]
#     weights = weights[ord]
#     p = cumsum(weights)
#     p = p * length(x) / (length(x) + 1) + 0.5 / (length(x) + 1)
#
#     # browser()
#     # denominator = weights * length(x)
#     # numerator = rank(x) / (length(x) + 1) * denominator
#     # browser()
#
#     #browse_on_error({
#     cdf_fit = locfit::locfit.raw(
#       x = locfit::lp(x, nn = nn, deg = 1),
#       y = p,
#       family = "binomial"
#       # maxk = 500
#     )
#     #})
#     knotp = stats::predict(cdf_fit, knotx)
#     # get rid of any points which are decreasing:
#     decreasing = (knotp < cummax(knotp))
#     knotx = knotx[!decreasing]
#     knotp = knotp[!decreasing]
#     # if (any(decreasing)) {
#     #   warning("Smoothed empirical CDF was not strictly monotonic.")
#     # }
#   }
#
#   return(
#     empirical_cdf(
#       x = knotx,
#       p = knotp,
#       link = link,
#       smooth = smooth,
#       name = name,
#       ...
#     )
#   )
# }

# Empirical CDF fitting ----

#' Fit a piecewise logit transformed linear model to a CDF
#'
#' This fits a CDF and quantile function to cumulative probabilities in a transformed space. X value
#' transformation is specified in the `link` parameter and is either something
#' like "log", "logit", etc or can also be specified as the logit transformed cdf and quantile function from a statistical distribution.
#'
#' The empirical distribution fitted is a piecewise linear or monotonically
#' increasing spline fit to CDF in transformed X and logit Y space. The end
#' points are linearly interpolated in this space to the `tail_p`th quantile.
#' The function can fit data provided as `x, P(X<=x)` pairs.
#'
#' This function imputes tails of distributions. Given perfect data as samples
#' or as quantiles it should recover the tail
#'
#' @inheritParams empirical
#' @param smooth fits the empirical distribution with a pair of splines for CDF
#'   and quantile function, creating a mostly smooth density. This smoothness
#'   comes at the price of potential over-fitting and will produce small
#'   differences between `p` and `q` functions such that `x=p(q(x))` is no
#'   longer exactly true. Setting this to false will replace this with a
#'   piecewise linear fit that is not smooth in the density, but is exact in
#'   forward and reverse transformation.
#' @param ... not used
#'
#' @returns a `dist_fns` S3 object containing functions `p()` for CDF, `q()`
#'   for quantile, and `r()` for a RNG, and `d()` for density. The density
#'   function may be discontinuous.
#' @concept empirical
#' @export
#'
#' @unit
#'
#' #from cdf:
#' xs = c(2,3,6,9)
#' ps = c(0.1,0.4,0.6,0.95)
#' e = empirical_cdf(xs, ps, link="log")
#'
#' testthat::expect_equal(e$p(xs), ps)
#' testthat::expect_equal(e$q(ps), xs)
#'
#' # quantiles:
#' p = c(0.025,0.05,0.10,0.25,0.5,0.75,0.9,0.95,0.975)
#' q = stats::qgamma(p, shape=2)
#' shape2_gamma = as.dist_fns(pgamma, shape=2)
#' gemp = empirical_cdf(q,p,link = shape2_gamma)
#' withr::with_seed(123, {
#'   testthat::expect_equal(mean(gemp$r(100000)),2, tolerance=0.01)
#'   testthat::expect_equal(sd(gemp$r(100000)), sqrt(2), tolerance=0.01)
#' })
#'
#' # With perfect input can recover the underlying distribution including tails:
#' tmp = empirical_cdf(x=seq(0.01,0.99,0.01),p=seq(0.01,0.99,0.01),link = as.dist_fns(punif,0, 1), knots = 100)
#' testthat::expect_equal(
#'   tmp$q(c(0.01, 0.1, 0.25, 0.75, 0.9, 0.99)),
#'   c(0.01, 0.1, 0.25, 0.75, 0.9, 0.99),
#'   tolerance = 0.002
#' )
#'
#' # bimodal with log link and end defined
#' tmp3 = empirical_cdf(x = 1:7, c(0.1,0.3,0.4,0.4,0.5,0.9,1),link="log")
#' testthat::expect_equal(
#'   tmp3$p(-1:8),
#'   c(0, 0, 0.1, 0.3, 0.4, 0.4, 0.5, 0.9, 1, 1),
#'   tolerance = 0.01
#' )
#'
#' testthat::expect_equal(
#'   tmp3$q(seq(0, 1, 0.2)),
#'   c(0, 1.63476839034733, 3, 5.05332769159444, 5.51891960240613, 7)
#' )
#'
#'
empirical_cdf = function(
  x,
  p,
  link = "ident",
  smooth = TRUE,
  name = NULL,
  ...
) {
  if (is.dist_fns(link)) {
    qx_from_x = link$p
    # TODO: dqx_from_x = link$d
    x_from_qx = link$q
    support = c(link$q(0), link$q(1))
  } else {
    link = as.link_fns(link)
    support = link$support
    if (identical(link$support, c(0, 1))) {
      # in Q-Q space already
      qx_from_x = function(x) x
      # TODO: dqx_from_x = ??? 1
      x_from_qx = function(qx) qx
    } else {
      # In arbitrary links we need to shift the centre and range of the
      # transformed distribution to prevent underflow when putting into logit
      # space, and allowing a more natural fit.

      # This is done with a piecewise linear model on the quantiles to look for
      # median and 35%-65% quantiles (approx 1 SD)

      # fit a normal to quantiles:
      qmu = stats::approx(p, link$trans(x), 0.5, ties = "mean")$y
      qsd = diff(
        stats::approx(
          p,
          link$trans(x),
          stats::pnorm(c(-0.5, 0.5)),
          ties = "mean"
        )$y
      )

      f = g = NULL
      qx_from_x = carrier::crate(
        function(x) {
          # TODO: this does no checking which may be a bad idea.
          x2 = suppressWarnings(f(x))
          return(stats::pnorm(x2, !!qmu, !!qsd))
        },
        f = link$trans
      )
      # TODO: dqx_from_x = dnorm(f(x)) * f_dash(x)
      # f_dash - link$ddxtrans
      x_from_qx = carrier::crate(
        function(q) {
          x2 = stats::qnorm(q, !!qmu, !!qsd)
          return(g(x2))
        },
        g = link$inv
      )
    }
  }

  # fitting from CDF in a Q-Q space

  qy = p
  qx = qx_from_x(x)

  order_by = order(qx, qy)
  qx = qx[order_by]
  qy = qy[order_by]

  # Set some absolute limits in Q-Q space

  # get rid of duplicates
  # this keeps the highest other value of a duplicated pair
  # so c(0.1,0.2,0.2,0.4) and c(0.1,0.2,0.3,0.4) will give
  # c(0.1,0.2,0.4) & c(0.1,0.3,0.4)
  # This is in line with the idea of the input being a CDF and the P(X<=x)

  # inputs have been sorted at this point, so qx is increasing or flat.
  # qx must be unique (otherwise we have a infinite PDF at the point)
  while (anyDuplicated(qx)) {
    qx[duplicated(qx)] = qx[duplicated(qx)] + sqrt(sqrt(.Machine$double.eps))
  }
  # because we need to invert it is also a bad idea to have completely
  # flat qy
  while (anyDuplicated(qy)) {
    qy[duplicated(qy)] = qy[duplicated(qy)] + sqrt(sqrt(.Machine$double.eps))
  }

  # qx was sorted so is increasing
  # in a CDF qy must be strictly flat or increasing.
  invalid = qy < cummax(qy)
  if (any(invalid)) {
    warning("CDF was not stricly increasing. Ignoring invalid points.")
  }
  qx = qx[!invalid]
  qy = qy[!invalid]

  # linear prediction of the tails is not directly possible in Q-Q unlike
  # logit logit but we can estimate where a point might be if current head / tail
  # rate of change of QQ is maintained and shifting this more centrally.
  # however we know the fixed end points of 0,0 and 1,1 if they have not been
  # defined in the CDF as hard truncation points:

  if (min(qx) > 0 & min(qy) > 0) {
    qy0 = .interpolate(utils::head(qy, 2), utils::head(qx, 2), 0)
    qx0 = .interpolate(utils::head(qx, 2), utils::head(qy, 2), 0)
    # shift centrally
    if (qx0 <= 0) {
      qx0 = qx[1] / 5
    }
    if (qy0 <= 0) {
      qy0 = qy[1] / 5
    }

    qx = c(0, qx0, qx)
    qy = c(0, qy0, qy)
  }

  if (max(qx) < 1 & max(qy) < 1) {
    qy1 = .interpolate(utils::tail(qy, 2), utils::tail(qx, 2), 1)
    qx1 = .interpolate(utils::tail(qx, 2), utils::tail(qy, 2), 1)
    # shift centrally
    if (qx1 >= 1) {
      qx1 = 1 - (1 - utils::tail(qx, 1)) / 5
    }
    if (qy1 >= 1) {
      qy1 = 1 - (1 - utils::tail(qy, 1)) / 5
    }

    qx = c(qx, qx1, 1)
    qy = c(qy, qy1, 1)
  }

  if (isFALSE(smooth)) {
    qy_from_qx = carrier::crate(function(qx2) {
      stats::approx(
        x = !!qx,
        y = !!qy,
        xout = qx2,
        yleft = 0,
        yright = 1
      )$y
    })
    qx_from_qy = carrier::crate(function(qy2) {
      stats::approx(
        x = !!qy,
        y = !!qx,
        xout = qy2,
        yleft = 0,
        yright = 1
      )$y
    })
  } else {
    qy_from_qx_spl = .monotonicpolyspline(x = qx, y = qy) # cdf zy from zx
    # qx_from_qy_spl = splines::backSpline(qy_from_qx_spl) # quantile fn zx from zy
    qx_from_qy_spl = .monotonicpolyspline(x = qy, y = qx) # quantile fn zx from zy

    qy_from_qx = carrier::crate(function(qx2) {
      tmp = rep(NA, length(qx2))
      non.na = which(!is.na(qx2))
      tmp[non.na] = stats::predict(!!qy_from_qx_spl, qx2[non.na])$y
      tmp[tmp < 0] = 0
      tmp[tmp > 1] = 1
      return(tmp)
    })
    qx_from_qy = carrier::crate(function(qy2) {
      tmp = rep(NA, length(qy2))
      non.na = which(!is.na(qy2))
      tmp[non.na] = stats::predict(!!qx_from_qy_spl, qy2[non.na])$y
      tmp[tmp < 0] = 0
      tmp[tmp > 1] = 1
      return(tmp)
    })

    # TODO: in theory the differential is possible to find using predict(deriv=1)
    # however we would need to think through chain rule of expit(predict(trans(x)))
    # which will probably require differential of link trans.
    # dqydqx_from_qx =
    # need to incorporate this into the crate also.
  }

  qfn = minx = maxx = NULL
  out = .super_crate(
    # support
    minx = support[1],
    maxx = support[2],
    #transforms
    qy_from_qx = qy_from_qx,
    qx_from_x = qx_from_x,
    qx_from_qy = qx_from_qy,
    x_from_qx = x_from_qx,
    #functions
    .fns = list(
      pfn = function(q) {
        ifelse(
          q <= minx,
          0,
          ifelse(q >= maxx, 1, qy_from_qx(qx_from_x(q)))
        )
      },
      qfn = function(p) {
        ifelse(p < 0 | p > 1, NaN, x_from_qx(qx_from_qy(p)))
      },
      rfn = function(n) {
        qfn(stats::runif(n))
      }
    )
  )

  return(
    new_dist_fns(
      name = if (is.null(name)) "empirical" else name,
      pfn = out$pfn,
      qfn = out$qfn,
      rfn = out$rfn,
      knots = dplyr::tibble(x = x_from_qx(qx), p = qy),
      smooth = smooth
    )
  )
}


#' Interpolate from 2 coordinates and predict output
#'
#' @param y the y values. length 2.
#' @param x the x values. length 2.
#' @param new_x the values to predict at (vectorised). If null will return a
#'   function that can predict new values.
#'
#' @returns the interpolated value or a prediction function
#' @keywords internal
#'
#' @unit
#' .interpolate(c(0,1),c(0,1),5)
.interpolate = function(y, x, new_x = NULL) {
  nas = is.na(x) | is.na(y)
  x = x[!nas]
  y = y[!nas]
  if (length(x) != length(y) || length(unique(x)) != 2) {
    stop("coordinates must be length 2")
  }
  fn = function(x2) {
    (y[2] - y[1]) / (x[2] - x[1]) * (x2 - x[1]) + y[1]
  }
  if (is.null(new_x)) {
    return(fn)
  }
  return(fn(new_x))
}


#' Create a `polySpline` object from a `stats::splinefun` call.
#'
#' @inheritParams stats::splinefun
#'
#' @details This function converts the output of `stats::splinefun`, to a
#'   a polynomial spline object from the `splines` package.
#'
#' @return An object of class `polySpline`.
#' @keywords internal
#' @unit
#'
#' # strictly increasing
#' spl = .monotonicpolyspline(1:10, log(1:10))
#' testthat::expect_equal(
#'   predict(spl,4.5)$y,
#'   log(4.5),
#'   tolerance = 0.01
#' )
#'
#' # strictly decreasing
#' spl2 = .monotonicpolyspline(1:10, -log(1:10))
#' testthat::expect_equal(
#'   predict(spl2,4.5)$y,
#'   -log(4.5),
#'   tolerance = 0.01
#' )
#'
#' # not monotonic
#' testthat::expect_error(
#'   {
#'     .monotonicpolyspline(-4:4, (-4:4)^2)
#'   },
#'   "Data is not monotonic.",
#'   fixed = TRUE
#' )
#'
.monotonicpolyspline = function(
  x,
  y
) {
  stopifnot(length(x) == length(y))
  if (any(is.na(y))) {
    browser()
  }

  o = order(x, y)
  x = x[o]
  y = y[o]

  mono = !any(x < cummax(x)) &&
    (!any(y < cummax(y)) || !any(y > cummin(y)))

  if (!mono) {
    stop("Data is not monotonic.", call. = FALSE)
  }

  # Y ties are OK:
  # e.g. plot(splinefun(1:4, c(1,2,2,4),method = "monoH.FC"),from = 1,to=4)
  # Can't have X ties (infinite gradient)
  # e.g. plot(splinefun(c(1,2,2+sqrt(.Machine$double.eps),4),1:4,method = "monoH.FC"),from = 1,to=4)
  while (anyDuplicated(x)) {
    x[duplicated(x)] = x[duplicated(x)] + sqrt(.Machine$double.eps)
  }

  sf = stats::splinefun(x, y, method = "monoH.FC")

  m = environment(sf)$m
  delta = environment(sf)$dx
  x = environment(sf)$x0
  y = environment(sf)$y0

  n <- length(y)
  y_i <- y[-n]
  y_ip1 <- y[-1]
  delta2 <- delta^2
  delta3 <- delta^3

  m_i <- m[-n] * delta
  m_ip1 <- m[-1] * delta

  ct_0 <- y_i
  ct_1 <- m_i
  ct_2 <- (-3 * y_i - 2 * m_i + 3 * y_ip1 - m_ip1)
  ct_3 <- (2 * y_i + m_i - 2 * y_ip1 + m_ip1)

  coefficients = rbind(
    cbind(ct_0, ct_1 / delta, ct_2 / delta2, ct_3 / delta3),
    matrix(c(utils::tail(y, 1), utils::tail(m, 1), 0, 0), nrow = 1)
  )
  colnames(coefficients) = c("constant", "linear", "quadratic", "cubic")

  spl <- structure(
    list(
      knots = x,
      coefficients = coefficients
    ),
    class = c("npolySpline", "polySpline", "spline")
  )
  return(spl)
}


# Empirical weighted data fitting ----

#' Quantile from weighted data
#'
#' @param p the quantile probabilities
#' @param names should the result be named
#' @param na.rm remove NAs from observations (and weights)
#'
#' @param x the observations
#' @param w importance weights. Must be some non zero weights, and all must be
#'   finite. They don't need to sum to one.
#' @param link a link function mapping observations to infinite support
#' @keywords internal
#' @name wt_quant
NULL


#' Fit a piecewise logit transformed linear model to weighted data
#'
#' This fits a CDF and quantile function to ranked data in a
#' transformed space. X value transformation is specified in the `link`
#' parameter and is either something like "log", "logit", etc.
#'
#' The empirical distribution fitted is a piecewise linear in z transformed X
#' and logit Y space. The evaluation points are linearly interpolated in this
#' space given a bandwidth for interpolation.
#'
#' This function imputes tails of distributions. Given perfect data as samples
#' or as quantiles it should recover the tail
#'
#' @inheritParams empirical
#' @param bw a bandwidth expressed in terms of the probability width, or proportion
#'   of observations.
#' @concept empirical
#' @export
#'
#' @returns a `dist_fns` S3 object that function that contains statistical
#'   distribution functions for this data.
#'
#' @unit
#'
#' # from samples:
#' withr::with_seed(123,{
#'  e2 = empirical_data(rnorm(10000), bw=0.1)
#'  testthat::expect_equal(e2$p(-5:5), pnorm(-5:5), tolerance=0.01)
#'  testthat::expect_equal(e2$d(-5:5), dnorm(-5:5), tolerance=0.05)
#'  testthat::expect_equal(e2$q(seq(0,1,0.1)), qnorm(seq(0,1,0.1)), tolerance=0.025)
#' })
#'
#'
#' # Construct a normal using a sequence and density as weight.
#' e7 = empirical_data(
#'   x=seq(-10,10,length.out=1000),
#'   w=dnorm(seq(-10,10,length.out=1000))
#' )
#' testthat::expect_equal(e7$p(-5:5), pnorm(-5:5), tolerance=0.01)
#'
empirical_data = function(
  x,
  w = NULL,
  link = "identity",
  ...,
  name = NULL,
  bw = NULL
) {
  link = as.link_fns(link)

  #TODO: review this:
  bw_min = .p_from_n(length(x), 40, 1000)
  if (is.null(bw) || bw < bw_min) {
    bw = bw_min
  }

  # dist = .logit_z_interpolation(link$trans(x), w, ..., bw=bw)
  dist = .logit_z_locfit(link$trans(x), w, ..., bw = bw)

  if (is.null(name)) {
    name = if (link$name == "I") {
      "empirical"
    } else {
      sprintf("empirical (%s link))", link$name)
    }
  }

  qfn = NULL
  out = .super_crate(
    dist = dist,
    link = link,
    .fns = list(
      pfn = function(q) {
        q2 = link$trans(q)
        p2 = dist$p(q2)
        return(p2)
      },
      qfn = function(p) {
        q2 = dist$q(p)
        q = link$inv(q2)
        return(q)
      },
      rfn = function(n) {
        qfn(stats::runif(n))
      },
      dfn = function(x) {
        x2 = link$trans(x)
        tmp1 = dist$d(x2)
        tmp2 = link$ddxtrans(x)
        ifelse(
          (tmp1 == 0 & !is.finite(tmp2) | !is.finite(tmp1) & tmp2 == 0),
          0,
          # chain rule
          tmp1 * tmp2
        )
      }
    )
  )

  return(new_dist_fns(
    name = name,
    pfn = out$pfn,
    qfn = out$qfn,
    rfn = out$rfn,
    dfn = out$dfn
  ))
}


#' Weighted distribution function interpolation in a logit z space
#'
#' Weighted cumulative probabilities are mapped to a logit space, data is transformed
#' to Z space (assumes support is -Inf..Inf). This does not use a link function
#' and the resulting interpolation functions are not vectorised. Importance weighting
#' is done during CDF construction. Prediction
#' is done using a weighted linear interpolation of nearby points. Weighting
#' for interpolation is a distance based gaussian kernel from data points to
#' interpolation point. OOB interpolation is supported.
#'
#' @inheritParams empirical
#' @param bw a bandwidth expressed in terms of the probability width, or proportion
#'   of observations.
#'
#' @returns a function that will predict a quantile assuming infinite support
#' @keywords internal
.logit_z_interpolation = function(x, w = NULL, bw = NULL) {
  if (is.null(w)) {
    w = rep(1, length(x))
  }

  if (length(x) != length(w)) {
    stop("Observations and weights must be the same length.")
  }

  if ((anyNA(x) || anyNA(w))) {
    cc = !is.na(x) && !is.na(w)
    x = x[cc]
    w = w[cc]
  }

  if (any(!is.finite(w))) {
    stop("Infinite weight detected.")
  }

  zerow = (w == 0)
  x = x[!zerow]
  w = w[!zerow]

  o = order(x)
  x = x[o]
  w = w[o]

  # map data into Z space
  x1mu = wmean(x, w)
  x1sigma = wsd(x, w)

  y = cumsum(w)
  y = (y - 0.5 * w) / sum(w)

  dup = duplicated(x, fromLast = TRUE) | duplicated(y)
  oob = y <= 0 | y >= 1
  x = x[!dup & !oob]
  y = y[!dup & !oob]

  if (length(x) < 9) {
    stop("Insufficient data for interpolation.")
  }

  x2 = (x - x1mu) / x1sigma
  y2 = stats::qlogis(y)

  n = length(x)

  return(
    # This is not a dist_fn because it is not vectorised and not link capable
    .super_crate(
      n = n,
      bw = bw,
      y = y,
      x = x,
      y2 = y2,
      x2 = x2,
      x1sigma = x1sigma,
      x1mu = x1mu,
      .fit_lm_1d = .fit_lm_1d,
      .fns = list(
        q = function(pv) {
          q2v = sapply(pv, function(p) {
            p2 = stats::qlogis(p)
            hbw = bw / 2
            i = which.min(abs(y - p + hbw))
            j = which.min(abs(y - p - hbw))
            while (j - i < 5) {
              hbw = hbw^0.5
              i = which.min(abs(y - p + hbw))
              j = which.min(abs(y - p - hbw))
            }

            # data weighting is included in CDF. kernel weighting for distance
            # gaussian normalised to bandwidth
            u = abs(y2[i:j] - p2)
            u = 3 * u / max(u)
            kernel_w = exp(-0.5 * u^2)

            q2 = .fit_lm_1d(
              y = x2[i:j],
              x = y2[i:j],
              new_x = p2,
              w = kernel_w
            )
            return(q2)
          })

          # map quantiles out of Z space:
          q1 = q2v * x1sigma + x1mu

          return(unname(q1))
        },
        p = function(qv) {
          hwin = floor(bw / 2 * n)

          p2v = sapply(qv, function(q) {
            q2 = (q - x1mu) / x1sigma

            m = which.min(abs(x - q))
            i = max(c(m - hwin, 1))
            j = min(c(m + hwin, n))
            while (j - i < 5) {
              i = max(c(i - 1, 1))
              j = min(c(j + 1, n))
            }
            # data weighting is included in CDF. kernel weighting for distance
            # kernel width:
            x2win = x2[i:j]
            dx2win = diff(range(x2win))
            # gaussian normalised to bandwidth
            u = abs((x2win - q2))
            u = 3 * u / max(u)
            kernel_w = exp(-0.5 * u^2)
            p2 = .fit_lm_1d(y = y2[i:j], x = x2win, new_x = q2, w = kernel_w)
          })
          # map probabilities out of logit space:
          p1 = 1 / (1 + exp(-p2v))

          return(unname(p1))
        },
        d = function(qv) {
          # double window for density
          hwin = floor(bw * n)

          d1v = sapply(qv, function(q) {
            m = which.min(abs(x - q))
            i = max(c(m - hwin, 1))
            j = min(c(m + hwin, n))
            while (j - i < 9) {
              i = max(c(i - 1, 1))
              j = min(c(j + 1, n))
            }

            q2 = (q - x1mu) / x1sigma

            # data weighting is included in CDF. kernel weighting for distance
            # kernel width:
            x2win = x2[i:j]
            dx2win = diff(range(x2win))
            # gaussian normalised to bandwidth
            u = abs((x2win - q2))
            u = 3 * u / max(u)
            kernel_w = exp(-0.5 * u^2)
            beta = .fit_lm_1d(y = y2[i:j], x = x2win, w = kernel_w)

            d2 = beta["m"]
            p2 = beta["m"] * q2 + beta["c"]

            # map probabilities out of logit space (chain rule):
            # expit'(F(q2)) * f(q2)
            d1 = exp(-p2) / (1 + exp(-p2))^2 * d2
          })

          return(unname(d1v) / x1sigma)
        }
      )
    )
  )
}


# Locfit log-z interpolation ----

#' Weighted distribution function interpolation in a logit z space
#'
#' Weighted cumulative probabilities are mapped to a logit space, data is transformed
#' to Z space (assumes support is -Inf..Inf). This does not use a link function
#' and the resulting interpolation functions are not vectorised. Importance weighting
#' is done during CDF construction. Prediction
#' is done using a weighted linear interpolation of nearby points. Weighting
#' for interpolation is a distance based gaussian kernel from data points to
#' interpolation point. OOB interpolation is supported.
#'
#' @inheritParams empirical
#' @param bw a bandwidth expressed in terms of the probability width, or proportion
#'   of observations.
#'
#' @returns a function that will predict a quantile assuming infinite support
#' @keywords internal
.logit_z_locfit = function(x, w = NULL, bw = NULL) {
  if (is.null(w)) {
    w = rep(1, length(x))
  }

  if (length(x) != length(w)) {
    stop("Observations and weights must be the same length.")
  }

  if ((anyNA(x) || anyNA(w))) {
    cc = !is.na(x) && !is.na(w)
    x = x[cc]
    w = w[cc]
  }

  if (any(!is.finite(w))) {
    stop("Infinite weight detected.")
  }

  zerow = (w == 0)
  x = x[!zerow]
  w = w[!zerow]

  o = order(x)
  x = x[o]
  w = w[o]

  # map data into Z space
  x1mu = wmean(x, w)
  x1sigma = wsd(x, w)

  y = cumsum(w)
  y = (y - 0.5 * w) / sum(w)

  dup = duplicated(x, fromLast = TRUE) | duplicated(y)
  oob = y <= 0 | y >= 1
  x = x[!dup & !oob]
  y = y[!dup & !oob]

  if (length(x) < 9) {
    stop("Insufficient data for interpolation.")
  }

  x2 = (x - x1mu) / x1sigma
  y2 = stats::qlogis(y)

  # n = length(x)

  # fit locfit models forward and backwards:
  cdf_fit = locfit::locfit.raw(
    x = locfit::lp(x2, nn = bw, deg = 1),
    y = y2
  )

  pdf_fit = locfit::locfit.raw(
    x = locfit::lp(x2, nn = bw, deg = 1),
    y = y2,
    deriv = 1
  )

  qf_fit = locfit::locfit.raw(
    x = locfit::lp(y2, nn = bw, deg = 1),
    y = x2
  )

  return(
    # This is not a dist_fn because it is not vectorised and not link capable
    .super_crate(
      qf_fit = qf_fit,
      cdf_fit = cdf_fit,
      pdf_fit = pdf_fit,
      x1sigma = x1sigma,
      x1mu = x1mu,
      .fns = list(
        q = function(p) {
          p2 = stats::qlogis(p)
          valid = is.finite(p2)
          q2 = p2 # copy infinites and NAs
          q2[valid] = stats::predict(qf_fit, p2[valid])
          # map quantiles out of Z space:
          q1 = q2 * x1sigma + x1mu
          return(unname(q1))
        },
        p = function(q) {
          q2 = (q - x1mu) / x1sigma
          valid = is.finite(q2)
          p2 = q2 # copy infinites and NAs
          p2[valid] = stats::predict(cdf_fit, q2[valid])
          # map probabilities out of logit space:
          p1 = 1 / (1 + exp(-p2))
          return(unname(p1))
        },
        d = function(q) {
          q2 = (q - x1mu) / x1sigma
          valid = is.finite(q2)
          p2 = q2 # copy infinites and NAs
          d2 = rep(0, length(q2))
          p2[valid] = stats::predict(cdf_fit, q2[valid])
          d2[valid] = stats::predict(pdf_fit, q2[valid])
          # map probabilities out of logit space (chain rule):
          # expit'(F(q2)) * f(q2)
          d1 = ifelse(
            p2 == -Inf,
            0,
            exp(-p2) / (1 + exp(-p2))^2 * d2
          )
          return(unname(d1) / x1sigma)
        }
      )
    )
  )
}

# General data distribution functions ----

#' Calculate the skew of a set of data
#'
#' @param x a vector of observations
#' @param na.rm remove `NA`s?
#'
#' @returns the skew
#' @export
#' @concept empirical
#'
#' @examples
#' skew(stats::rnorm(1000))
#' skew(stats::rbeta(1000, 1, 8)) # positively (left) skewed
#' skew(stats::rbeta(1000, 8, 1)) # negatively (right) skewed
#'
#' skew(stats::rlnorm(1000))
skew = function(x, na.rm = FALSE) {
  mu = mean(x, na.rm = na.rm)
  sigma = stats::sd(x, na.rm = na.rm)
  if (na.rm) {
    x = x[!is.na(x)]
  }
  mean((x - mu)^3) / (mean((x - mu)^2))^1.5
}

#' Calculate the excess kurtosis of a set of data
#'
#' @param x a vector of observations
#' @param na.rm remove `NA`s?
#' @param excess if false calculates raw kurtosis rather than excess
#'
#' @returns the excess kurtosis
#' @export
#' @concept empirical
#'
#' @examples
#' kurtosis(stats::rnorm(1000))
#' kurtosis(stats::rpois(1000, 2)) # leptokurtic > 0 (usually)
#' kurtosis(stats::runif(1000)) # platykurtic: < 0
#'
#' kurtosis(stats::rlnorm(1000))
kurtosis = function(x, na.rm = FALSE, excess = TRUE) {
  mu = mean(x, na.rm = na.rm)
  sigma = stats::sd(x, na.rm = na.rm)
  if (na.rm) {
    x = x[!is.na(x)]
  }
  kappa = mean((x - mu)^4) / (mean((x - mu)^2))^2
  if (excess) {
    kappa = kappa - 3
  }
  return(kappa)
}


#' Weighted standard deviation
#'
#' @inheritParams empirical
#' @param na.rm remove NAs (default TRUE)
#'
#' @returns a standard deviation
#' @export
#' @concept empirical
#' @examples
#'
#' # unweighted:
#' wsd(x = stats::rnorm(1000))
#'
#' # weighted:
#' wsd(x = seq(-2,2,0.1), w = stats::dnorm(seq(-2,2,0.1)))
#'
wsd = function(x, w = NULL, na.rm = TRUE) {
  if (is.null(w)) {
    return(stats::sd(x, na.rm = na.rm))
  }
  if (na.rm) {
    nas = is.na(x) | is.na(w)
    x = x[!nas]
    w = w[!nas]
  }
  sqrt(stats::cov.wt(matrix(x, ncol = 1), wt = w)$cov[1, 1])
}

#' Weighted mean
#'
#' a simple alias for base `weighted,mean`
#'
#' @inheritParams empirical
#' @param na.rm remove NAs (default TRUE)
#'
#' @returns a standard deviation
#' @export
#' @concept empirical
#' @examples
#' #' # unweighted:
#' wmean(x = stats::rnorm(1000))
#'
#' # weighted:
#' wmean(x = seq(-2,2,0.1), w = stats::dnorm(seq(-2,2,0.1)))
#'
wmean = function(x, w = NULL, na.rm = TRUE) {
  if (is.null(w)) {
    return(mean(x, na.rm = na.rm))
  }
  if (na.rm) {
    nas = is.na(x) | is.na(w)
    x = x[!nas]
    w = w[!nas]
  }
  stats::weighted.mean(x, w)
}

#' Quantile from weighted data with link function support
#'
#' This quantile function has different order of parameters from
#' base quantile. It takes a weight and a link function specification which
#' allows us to define the support of the quantile function. It is
#' optimised for imputing the tail of distributions and not speed.
#'
#' This is a moderately expensive function to call (in memory terms), as it
#' needs to construct the whole quantile function. if there are multiple calls
#' consider using `empirical()` to build a quantile function and using that.
#'
#' @param p the probabilities for which to estimate quantiles from the data
#' @param x a set of observations
#' @inheritParams empirical
#' @param names name the resulting quantile vector
#' @param window the number of data points to include when estimating the
#'   quantile. The closest `window` points are picked and used as a distance
#'   weighted linear interpolation of the weighted CDF in logit-link space. This
#'   tends to give good results for extrapolating tails.
#' @concept empirical
#'
#' @returns a vector of quantiles
#' @export
#' @examples
#'
#' # unweighted:
#' wquantile(p = c(0.25,0.5,0.75), x = stats::rnorm(1000))
#'
#' # weighted:
#' wquantile(p = c(0.25,0.5,0.75), x = seq(-2,2,0.1), w = stats::dnorm(seq(-2,2,0.1)))
#'
#' @unit
#'
#' test = function(rfn,qfn, link,..., n = 100000, tol=1000/(n+10000)) {
#'   testthat::expect_equal(abs(
#'    unname(wquantile(c(0.025, 0.5, 0.975),rfn(100000,...),link=link,names=FALSE)-
#'     qfn(c(0.025, 0.5, 0.975), ...))
#'   ), c(0,0,0), tolerance=tol)
#' }
#'
#' withr::with_seed(123, {
#'
#'   test(stats::rnorm,stats::qnorm,"identity",n = 10000)
#'   test(stats::rnorm,stats::qnorm,"identity",mean=4,n = 10000)
#'   test(stats::rnorm,stats::qnorm,"identity",sd=3, n = 100000, tol=0.05)
#'
#'   test(stats::rnorm,stats::qnorm,"identity",n = 5000)
#'   test(stats::rnorm,stats::qnorm,"identity",n = 1000)
#'   test(stats::rnorm,stats::qnorm,"identity",n = 100)
#'   test(stats::rnorm,stats::qnorm,"identity",n = 30)
#'
#'   test(stats::rgamma,stats::qgamma,"log", 4,n = 10000)
#'   test(stats::rgamma,stats::qgamma,"log", 4,n = 5000)
#'   test(stats::rgamma,stats::qgamma,"log", 4,n = 1000)
#'   test(stats::rgamma,stats::qgamma,"log", 4, 3,n = 100)
#'   test(stats::rgamma,stats::qgamma,"log", 4,n = 30)
#'
#'   test(stats::runif,stats::qunif,as.link_fns(c(0,10)),0,10)
#'
#' })
#'
#'
wquantile = function(
  p,
  x,
  w = NULL,
  link = "identity",
  names = TRUE,
  window = 7
) {
  link = as.link_fns(link)

  if (is.null(w)) {
    w = rep(1, length(x))
  }

  if ((anyNA(x) || anyNA(w))) {
    cc = !is.na(x) && !is.na(w)
    x = x[cc]
    w = w[cc]
  }

  o = order(x)
  x = x[o]
  w = w[o]

  x1 = link$trans(x)

  x1mu = wmean(x1, w)
  x1sigma = wsd(x1, w)

  y = cumsum(w)
  y = (y - 0.5 * w) / sum(w)

  dup = duplicated(x, fromLast = TRUE) | duplicated(y)
  oob = y <= 0 | y >= 1
  x = x[!dup & !oob]
  y = y[!dup & !oob]

  if (length(x) < window * 2) {
    stop("Insufficient data for interpolation.")
  }

  x2 = (x1 - x1mu) / x1sigma
  y2 = stats::qlogis(y)

  n = length(x1)
  p2v = stats::qlogis(p)

  q2v = sapply(p2v, function(p2) {
    hwin = window %/% 2
    m = which.min(abs(y2 - p2))
    i = max(c(m - hwin, 1))
    j = min(c(m + hwin, n))
    while (j - i < window) {
      i = max(c(i - hwin, 1))
      j = min(c(j + hwin, n))
    }

    # data weighting is included in CDF. kernel weighting for distance
    # gaussian normalised to bandwidth
    u = abs(y2[i:j] - p2)
    u = 3 * u / max(u)
    kernel_w = exp(-0.5 * u^2)

    q2 = .fit_lm_1d(
      y = x2[i:j],
      x = y2[i:j],
      new_x = p2,
      w = kernel_w
    )
    return(q2)
  })

  q1 = q2v * x1sigma + x1mu
  q = link$inv(q1)

  if (names) {
    names(q) = sprintf("%1.3g%%", p * 100)
  }
  return(q)
}

# Utility ----

#' Fit a weighted 1D linear model and predict output
#'
#' If all the weights are NA they are ignored.
#'
#' @param y the y values. At least 2.
#' @param x the x values. At least 2.
#' @param ... must be empty
#' @param w weights (optional)
#' @param new_x (optional) the values to predict at (vectorised).
#'
#' @returns either a vector with intercept (`c`) and gradient (`m`) or
#'   a vector of predictions at points `new_x`
#' @keywords internal
#'
#' @unit
#' testthat::expect_equal(
#'   .fit_lm_1d(c(1, 2, 3, 4), c(5, 6, 7, 8)),
#'   c(c = -4, m = 1)
#' )
#'
#' testthat::expect_equal(
#'    .fit_lm_1d(c(0, 1), c(0, 1)),
#'    c(c = 0, m = 1)
#' )
#'
#' withr::with_seed(123,{
#'   x = 1:100
#'   y = 2*x+3+stats::runif(100,-0.1,0.1)
#'   testthat::expect_equal(
#'     .fit_lm_1d(y, x),
#'     c(c = 3, m = 2),
#'     tolerance = 0.01
#'   )
#' })
#'
#' # with prediction:
#' testthat::expect_error(
#'   {
#'     .fit_lm_1d(c(6, 5, 7, 10), c(1, 2, 3, 4), c(0, 1))
#'   },
#'   structure("`...` must be empty.", names = ""),
#'   fixed = TRUE
#' )
#'
#' testthat::expect_equal(
#'   .fit_lm_1d(c(0, 1), c(0, 1), new_x = -2:2),
#'   c(-2, -1, 0, 1, 2)
#' )
#'
.fit_lm_1d = function(y, x, ..., w = NULL, new_x = NULL) {
  rlang::check_dots_empty()
  if (length(x) != length(y)) {
    stop("unequal length LM coordinates")
  }

  if (is.null(w) || length(unique(w)) == 1 || all(is.na(w))) {
    w = NULL
  }

  missing = is.na(x) | is.na(y)

  if (!is.null(w)) {
    if (length(w) != length(x)) {
      stop("weights length does not match data length")
    }
    if (sum(w) == 0) {
      stop("sum of weights is zero.")
    }
    missing = missing | is.na(w)
    w = w[!missing]
  }

  # Handle missing values
  x = x[!missing]
  y = y[!missing]

  # Check for sufficient data
  if (length(unique(x)) < 2 || length(unique(y)) < 2) {
    stop("not enough data in LM")
  }

  # Prepare design matrix
  X = matrix(c(rep(1, length(x)), x), ncol = 2)
  Y = matrix(y, ncol = 1)

  if (is.null(w)) {
    # Ordinary least squares
    beta = solve(t(X) %*% X) %*% t(X) %*% Y
  } else {
    # Weighted least squares
    W = diag(w) # Diagonal weight matrix
    beta = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y
  }

  if (is.null(new_x)) {
    return(c(c = beta[1, 1], m = beta[2, 1]))
  } else {
    return(beta[2, 1] * new_x + beta[1, 1])
  }
}
