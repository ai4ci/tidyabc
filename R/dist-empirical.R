#' Fit a piecewise logit transformed linear model to cumulative data
#'
#' This fits a CDF and quantile function to data in a transformed space. X value
#' transformation is specified in the `link` parameter and is either something
#' like "log", "logit", etc or can also be specified as the logit transformed cdf and quantile function from a statistical distribution.
#'
#' The empirical distribution fitted is a piecewise linear or monotonically
#' increasing spline fit to data in transformed X and logit Y space. The end
#' points are linearly interpolated in this space to the `tail_p`th quantile.
#' The function can fit data provided either as samples in a vector of `x`, or
#' as quantiles specified as `x, P(X<=x)` pairs.
#'
#' This function imputes tails of distributions. Given perfect data as samples
#' or as quantiles it should recover the tail
#'
#' @param x either a vector of samples from a distribution `X` or cut-offs for
#'   cumulative probabilities when combined with `probs`
#' @param probs (optional) If present a vector the same length as x giving
#'   `P(X <= x)`.
#' @param link a link function. Either as name, or a `link_fns` S3 object. In
#'   the latter case this could be derived from a statistical distribution by
#'   `as.link_fns(<dist_fns>)`. This supports the use of a prior to define the
#'   support of the empirical function, and is designed to prevent tail
#'   truncation. Support for the updated quantile function will be the same as
#'   the provided prior.
#' @param knots for distributions from data how many points do we use to model
#'   the cdf? I recommend an uneven number, without a lot of data this will tend
#'   to overfit 9 knots for 1000 samples seems OK, max 7 for 250, 5 for 100.
#'   Less is usually more.
#' @param tail_p what is the minimum tail probability modelled.
#' @param weights if the distribution is from data then an importance weighting
#'   for each point can be supplied here. This will be normalised, so it is a
#'   relative importance with respect to the mean of the weights.
#' @param smooth fits the empirical distribution with a spline, and generates a
#'   back-spline for the inverse, creating a mostly smooth density. This
#'   smoothness comes at the price of potential over-fitting and will produce
#'   small differences between `p` and `q` functions such that `x=p(q(x))` is no
#'   longer exactly true. Setting this to false will replace this with a
#'   piecewise linear fit that is not smooth in the density, but is exact in
#'   forward and reverse transformation.
#' @param ... not used
#'
#' @returns a `dist_fns` S3 object containing 3 functions `p()` for CDF, `q()`
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
#' e = empirical(xs, ps, link="log")
#'
#' testthat::expect_equal(e$p(xs), ps)
#' testthat::expect_equal(e$q(ps), xs)
#'
#' # from samples:
#' withr::with_seed(123,{
#'  e2 = empirical(rnorm(10000),knots = 20)
#'  testthat::expect_equal(e2$p(-5:5), pnorm(-5:5), tolerance=0.01)
#'  testthat::expect_equal(e2$q(seq(0,1,0.1)), qnorm(seq(0,1,0.1)), tolerance=0.025)
#' })
#'
#' p2 = seq(0,1,0.1)
#' testthat::expect_equal( e2$p(e2$q(p2)), p2, tolerance = 0.001)
#'
#' # quantiles:
#' p = c(0.025,0.05,0.10,0.25,0.5,0.75,0.9,0.95,0.975)
#' q = stats::qgamma(p, shape=2)
#' shape2_gamma = as.dist_fns(pgamma, shape=2)
#' gemp = empirical(q,p,link = shape2_gamma)
#' withr::with_seed(123, {
#'   testthat::expect_equal(mean(gemp$r(100000)),2, tolerance=0.01)
#'   testthat::expect_equal(sd(gemp$r(100000)), sqrt(2), tolerance=0.01)
#' })
#'
#' # With perfect input can recover the underlying distribution including tails:
#' tmp = empirical(x=seq(0.01,0.99,0.01),link = as.dist_fns(punif,0, 1), knots = 100)
#' testthat::expect_equal(
#'   tmp$q(c(0.01, 0.1, 0.25, 0.75, 0.9, 0.99)),
#'   c(0.01, 0.1, 0.25, 0.75, 0.9, 0.99),
#'   tolerance = 0.002
#' )
#'
empirical = function(
  x,
  probs = NULL,
  link = "ident",
  knots = NULL,
  tail_p = 0.00001,
  weights = NULL,
  smooth = TRUE,
  name = NULL,
  ...
) {
  link = as.link_fns(link)
  # if (is.null(limits)) {
  #   limits = link$domain
  # }

  trans_fn = link$trans
  inv_fn = link$inv

  p_trans_fn = .logit
  p_inv_fn = .expit

  # p_trans_fn = qnorm
  # p_inv_fn = pnorm

  if (is.null(knots)) {
    knots = .knots_from_data(x)
  }
  knots = min(c(knots, length(x)))
  knots = knots %/% 2 * 2 + 1

  if (is.null(weights)) {
    weights = rep(1, length(x))
  }

  if (is.null(probs)) {
    # fitting from data

    q = 0.01
    minx = stats::quantile(x, q)
    maxx = stats::quantile(x, 1 - q)
    #knotx = qnorm(seq(q, 1 - q, length.out = knots), mean(x), sd(x))

    knotx = seq(minx, maxx, length.out = knots)
    h = (maxx - minx) / (knots / 2) * 1.25
    cdf_fit = locfit::locfit.raw(
      x = locfit::lp(x, h = h, deg = 1),
      # x = locfit::lp(x, nn = 1 / (knots * 5)),
      y = rank(x) / (length(x) + 1) * weights, # numerator
      family = "binomial",
      weights = weights #* length(x) # denominator
    )
    tmp = predict(cdf_fit, knotx)
    # browser()
    if (any(tmp != cummax(tmp))) {
      tmp = cummax(tmp)
      knotx = knotx[!duplicated(tmp)]
      tmp = tmp[!duplicated(tmp)]
      warning("Smoothed empirical CDF was not strictly monotonic.")
    }

    # TODO: can i extract the spline from the locfit without

    # dens = locfit::density.lf(
    #   x,
    #   weights = weights,
    #   width = (maxx - minx) / knots,
    #   n = knots * 10 + 1,
    #   from = minx,
    #   to = maxx,
    #   ...
    # )
    # x = dens$x
    # tmp = cumsum(dens$y)
    # tmp = tmp - utils::head(tmp, 1)
    y = tmp #/
    #   utils::tail(tmp, 1) *
    #   (1 - 2 * q) +
    #   q

    x2 = trans_fn(knotx)
    # x2 = trans_fn(x)
    y2 = p_trans_fn(y)

    predict_with = 2
  } else {
    # fitting from CDF

    y = probs[order(x)]
    x = x[order(x)]

    x2 = trans_fn(x)
    y2 = p_trans_fn(y)

    predict_with = 2
  }

  y2_eps_hi = p_trans_fn(1 - tail_p)
  y2_eps_lo = p_trans_fn(tail_p)

  x2_eps_lo = .predict_lm_1d(
    utils::head(x2, predict_with),
    utils::head(y2, predict_with),
    y2_eps_lo
  )
  # x2_eps_lo = max(c(x2_eps_lo), min(trans_fn(limits)))

  x2_eps_hi = .predict_lm_1d(
    utils::tail(x2, predict_with),
    utils::tail(y2, predict_with),
    y2_eps_hi
  )
  # x2_eps_hi = min(c(x2_eps_hi), max(trans_fn(limits)))

  # filter down size of data until it is knots long
  # if (is.null(probs)) {
  #   x2 = x2[seq_along(x2) %% 10 == 1]
  #   y2 = y2[seq_along(y2) %% 10 == 1]
  # }

  # browser()
  if (min(y2) > y2_eps_lo) {
    x2 = c(x2_eps_lo, x2)
    y2 = c(y2_eps_lo, y2)
  }

  if (max(y2) < y2_eps_hi) {
    x2 = c(x2, x2_eps_hi)
    y2 = c(y2, y2_eps_hi)
  }

  dup =
    abs(x2 - c(x2[-1], Inf)) < sqrt(.Machine$double.eps) |
    abs(y2 - c(y2[-1], Inf)) < sqrt(.Machine$double.eps)

  x2 = x2[!dup]
  y2 = y2[!dup]

  if (isFALSE(smooth)) {
    fwdfn = stats::approxfun(
      x2,
      y2,
      yleft = p_trans_fn(0),
      yright = p_trans_fn(1)
    )
    backfn = stats::approxfun(
      y2,
      x2,
      yleft = p_trans_fn(0),
      yright = p_trans_fn(1)
    )
    samplerfn = stats::approxfun(y2, x2, yleft = min(y2), yright = max(y2))
  } else {
    # bspl = .monotonicpolyspline(y2, x2) # quantile function
    # ispl = splines::backSpline(bspl)

    ispl = .monotonicpolyspline(x2, y2) # cdf
    bspl = splines::backSpline(ispl)

    fwdfn = carrier::crate(function(vx2) {
      ifelse(
        vx2 < !!min(x2),
        !!p_trans_fn(0),
        ifelse(vx2 > !!max(x2), !!p_trans_fn(1), stats::predict(!!ispl, vx2)$y)
      )
    })
    backfn = carrier::crate(function(vy2) {
      ifelse(
        vy2 < !!min(y2),
        !!p_trans_fn(0),
        ifelse(vy2 > !!max(y2), !!p_trans_fn(1), stats::predict(!!bspl, vy2)$y)
      )
    })
    samplerfn = carrier::crate(function(vy2) {
      ifelse(
        vy2 < !!min(y2),
        !!min(y2),
        ifelse(vy2 > !!max(y2), !!max(y2), stats::predict(!!bspl, vy2)$y)
      )
    })
    # TODO: in theory the differential is possible to find using predict(deriv=1)
    # however we would need to think through chain rule of expit(predict(trans(x)))
    # which will probably require differential of trans, which we don't have
    # analytically anyway.
  }

  x = y = y2_eps = window = predict_with = NULL
  return(
    # TODO: look at supporting log option.
    new_dist_fns(
      name = if (is.null(name)) "empirical" else name,
      pfn = carrier::crate(
        function(q) {
          p_inv_fn(fwdfn(trans_fn(q)))
        },
        fwdfn = fwdfn,
        trans_fn = trans_fn,
        p_inv_fn = p_inv_fn
      ),
      qfn = carrier::crate(
        function(p) {
          inv_fn(backfn(p_trans_fn(p)))
        },
        backfn = backfn,
        inv_fn = inv_fn,
        p_trans_fn = p_trans_fn
      ),
      rfn = carrier::crate(
        function(n) {
          inv_fn(samplerfn(p_trans_fn(stats::runif(n))))
        },
        samplerfn = samplerfn,
        inv_fn = inv_fn,
        p_trans_fn = p_trans_fn
      ),
      knots = dplyr::tibble(x = inv_fn(x2), p = p_inv_fn(y2))
    )
  )
}


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
#' skew(rnorm(1000))
#' skew(rbeta(1000, 1, 8)) # positively (left) skewed
#' skew(rbeta(1000, 8, 1)) # negatively (right) skewed
#'
#' skew(rlnorm(1000))
skew = function(x, na.rm = FALSE) {
  mu = mean(x, na.rm = na.rm)
  sigma = sd(x, na.rm = na.rm)
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
#' kurtosis(rnorm(1000))
#' kurtosis(rpois(1000, 2)) # leptokurtic > 0 (usually)
#' kurtosis(runif(1000)) # platykurtic: < 0
#'
#' kurtosis(rlnorm(1000))
kurtosis = function(x, na.rm = FALSE, excess = TRUE) {
  mu = mean(x, na.rm = na.rm)
  sigma = sd(x, na.rm = na.rm)
  if (na.rm) {
    x = x[!is.na(x)]
  }
  kappa = mean((x - mu)^4) / (mean((x - mu)^2))^2
  if (excess) {
    kappa = kappa - 3
  }
  return(kappa)
}

# Empirical function utilities ----

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
  if (length(unique(x)) != length(unique(y)) || length(unique(y)) != 2) {
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

#' Fit a 1D linear model and predict output
#'
#' @param y the y values. At least 2.
#' @param x the x values. At least 2.
#' @param new_x the values to predict at (vectorised). If null will return a
#'   function that can predict new values.
#'
#' @returns the predicted value or a prediction function
#' @keywords internal
#'
#' @unit
#' .predict_lm_1d(c(6,5,7,10), c(1,2,3,4), c(0,1))
#' .predict_lm_1d(c(0,1), c(0,1), -2:2)
.predict_lm_1d = function(y, x, new_x = NULL) {
  if (length(y) == 2 && length(x) == 2) {
    return(.interpolate(y, x, new_x))
  }
  if (length(x) != length(y)) {
    stop("unequal length LM coordinates")
  }
  if (length(unique(x)) < 2 || length(unique(y)) < 2) {
    stop("not enough data in LM")
  }
  missing = is.na(x) | is.na(y)
  x = x[!missing]
  y = y[!missing]
  x = matrix(c(rep(1, length(x)), x), ncol = 2)
  y = matrix(y, ncol = 1)
  beta = solve(t(x) %*% x) %*% t(x) %*% y
  fn = function(x2) return(beta[1, 1] + beta[2, 1] * x2)
  if (is.null(new_x)) {
    return(fn)
  }
  return(fn(new_x))
}


# fast
.logit = stats::qlogis

# slow (but faster that plogis)
.expit = function(x) {
  return(1 / (1 + exp(-x)))
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
  y,
  ties = mean
) {
  stopifnot(length(x) == length(y))
  if (any(is.na(y))) {
    browser()
  }

  o = order(x)
  x = x[o]
  y = y[o]

  dup =
    abs(x - c(x[-1], Inf)) < sqrt(.Machine$double.eps) |
    abs(y - c(y[-1], Inf)) < sqrt(.Machine$double.eps)

  x = x[!dup]
  y = y[!dup]

  if (
    !(all(utils::head(y, -1) > y[-1]) ||
      all(utils::head(y, -1) < y[-1]))
  ) {
    stop("Data is not monotonic.", call. = FALSE)
  }

  sf = stats::splinefun(x, y, method = "monoH.FC", ties)
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


.knots_from_data = function(x) {
  round(log(length(x)) / 2) * 2 + 1
}
