#' `dist_fns` S3 Class operations
#'
#' A function wrapping a single (or set) of parametrised distributions and
#' provides access to the quantile, cumulative probability and random functions
#' of that specific distribution. Parametrisation is handled on construction.
#'
#' @param x a `dist_fns` S3 object
#' @param ... passed onto methods
#' @name s3_dist_fns
NULL

#' @describeIn s3_dist_fns Create a new `dist_fns` S3 object
#' @param name a name for the distribution
#' @param pfn the cumulative probability function
#' @param qfn the quantile function
#' @param rfn the RNG function
#' @param knots in empirical distributions this holds details on the knot points
#' @returns a new `dist_fns` S3 object
#' @keywords internal
#' @concept dist_fns_s3
new_dist_fns = function(name, pfn, qfn, rfn, ..., knots = NULL) {
  dots = rlang::list2(...)
  dots = dots[!names(dots) %in% c("lower.tail", "log.p", "p", "q", "n")]
  return(
    structure(
      list(
        name = name,
        p = function(q, lower.tail = TRUE, log.p = FALSE) {
          tmp = do.call(
            pfn,
            args = c(list(q = q), dots)
          )
          if (!lower.tail) {
            tmp = 1 - tmp
          }
          if (log.p) {
            tmp = log(tmp)
          }
          return(tmp)
        },
        q = function(p, lower.tail = TRUE, log.p = FALSE) {
          if (log.p) {
            p = exp(p)
          }
          if (!lower.tail) {
            p = 1 - p
          }
          tmp = do.call(
            qfn,
            args = c(list(p = p), dots)
          )
          return(tmp)
        },
        r = function(n) {
          do.call(rfn, args = c(list(n = n), dots))
        }
      ),
      knots = knots,
      class = c("dist_fns")
    )
  )
}

#' Format a `dist_fns` S3 object
#'
#' @param name description
#'
#' @export
#' @returns a character value
#' @concept dist_fns_s3
format.dist_fns = function(x, ..., digits = 3) {
  return(sprintf(
    "Dist: %s; Median (IQR) %s [%s \u2014 %s]",
    x$name,
    format(x$q(0.5), ..., digits = digits),
    format(x$q(0.25), ..., digits = digits),
    format(x$q(0.75), ..., digits = digits)
  ))
}

#' @describeIn s3_dist_fns Format a `dist_fns` S3 object
#' @export
#' @param steps resolution of the plot
#' @keywords internal
#' @returns a ggplot
#' @concept dist_fns_s3
plot.dist_fns = function(x, ..., steps = 200) {
  ys = seq(0, 1, length.out = steps + 1)
  xs = x$q(ys)
  dxs = utils::tail(xs, -1) - utils::head(xs, -1)
  dys = (1 / steps) / dxs
  ggplot2::ggplot(
    data.frame(x = xs, y = c(0, dys)),
    ggplot2::aes(x = x, y = y)
  ) +
    ggplot2::geom_step(direction = "vh")
}

#' @export
#' @describeIn as.dist_fns From a name
#' @concept dist_fns_s3
#' @param x the name of a distribution
#' @param ... the distribution parameters
#' @unit
#' x = as.dist_fns("norm")
#' testthat::expect_equal(x$p(-2:2), pnorm(-2:2))
as.dist_fns.character = function(x, ...) {
  distr = x
  rawfns = lapply(c("p", "q", "r"), function(type) {
    distname <- paste(type, distr, sep = "")
    pkg = rlang::caller_env()

    if (!exists(distname, mode = "function", envir = pkg)) {
      stop(
        "The ",
        distname,
        " function must be defined in ",
        environmentName(pkg)
      )
    }
    return(get(distname, mode = "function", envir = pkg))
  })

  return(
    new_dist_fns(distr, rawfns[[1]], rawfns[[2]], rawfns[[3]], ...)
  )
}

#' @export
#' @describeIn as.dist_fns From a statistical function
#' @concept dist_fns_s3
#' @param x a function from a statistical family, e.g. `stats::rnorm`
#' @param ... parameters passed to the function family defined in `x`
#' @unit
#' x = as.dist_fns(stats::rnorm)
#' testthat::expect_equal(x$p(-2:2), pnorm(-2:2))
as.dist_fns.function = function(x, ...) {
  distr = deparse(substitute(x)) #rlang::enexpr(x)
  pkg = environment(x)
  # distr = rlang::as_label(distr)
  name = gsub("^[^:]+::[dpqr]", "", distr, perl = TRUE)
  rawfns = lapply(c("p", "q", "r"), function(type) {
    if (grepl("::", distr)) {
      distname = gsub("^[^:]+::[dpqr]", type, distr, perl = TRUE)
    } else {
      distname = gsub("^[dpqr]", type, distr, perl = TRUE)
    }
    if (!exists(distname, mode = "function", envir = pkg)) {
      stop(
        "The ",
        distname,
        " function must be defined in ",
        environmentName(pkg)
      )
    }
    return(get(distname, mode = "function", envir = pkg))
  })

  return(
    new_dist_fns(name, rawfns[[1]], rawfns[[2]], rawfns[[3]], ...)
  )
}

#' @param x a `dist_fns` S3 class
#' @returns a numeric or other orderable item
#' @noRd
.order_dist_fns = function(x) {
  stop("ordering not supported")
}

#TODO: implement additional object specific exported functions

# S3 dist_fns class ----

#' @describeIn s3_dist_fns_list Create an empty `dist_fns_list`
#' @export
#' @returns an empty `dist_fns_list`
#' @concept dist_fns_s3
dist_fns = function() {
  return(as.dist_fns_list(NULL))
}

# removes this class from the object - internal use
.unclass_dist_fns = function(x) {
  class(x) <- setdiff(class(x), "dist_fns")
  return(x)
}


#' Create a `dist_fns` object
#'
#' @export
#' @returns a `dist_fns` S3 object
#' @concept dist_fns_s3
as.dist_fns = function(x, ...) {
  UseMethod("as.dist_fns", x)
}

#' @export
#' @concept dist_fns_s3
as.dist_fns.default = function(x, ...) {
  if (is.dist_fns(x)) {
    return(x)
  }
  stop("Don't know how to convert a `", class(x)[[1]], "` to a `dist_fns`.")
}

#' @describeIn s3_dist_fns Cast a `dist_fns` S3 object to a character vector
#' @export
#' @keywords internal
#' @returns a character vector
#' @concept dist_fns_s3
as.character.dist_fns = function(x, ...) {
  return(format(x, ...))
}

#' @describeIn s3_dist_fns Print a `dist_fns` S3 object
#' @export
#' @keywords internal
#' @returns nothing
#' @concept dist_fns_s3
print.dist_fns = function(x, ...) {
  cat(suppressWarnings(format(x, ...)))
  invisible(NULL)
}

#' @describeIn s3_dist_fns Print a `dist_fns` S3 object in an Rd document
#' @exportS3Method knitr::knit_print dist_fns
#' @keywords internal
#' @returns an `as-is` knitr chunk
#' @concept dist_fns_s3
knit_print.dist_fns = function(x, ...) {
  structure(format(x), class = "knit_asis")
}

#' @describeIn s3_dist_fns `dist_fns` S3 objects length is always 1
#' @export
#' @keywords internal
#' @returns 1 always
#' @concept dist_fns_s3
length.dist_fns = function(x, ...) {
  return(1)
}

#' @exportS3Method pillar::type_sum dist_fns
type_sum.dist_fns = function(x, ...) {
  abbreviate("dist_fns", 3)
}

#' @describeIn s3_dist_fns Check is a `dist_fns` S3 object
#' @export
#' @returns TRUE or FALSE
#' @concept dist_fns_s3
is.dist_fns = function(x) {
  return(inherits(x, "dist_fns"))
}

#' @describeIn s3_dist_fns Extract named attribute from a `dist_fns`
#' @param y item to retrieve
#' @returns an attribute value for `x`
#' @export
`@.dist_fns` = function(x, y) {
  if (is.character(y)) {
    ylab = y
  } else {
    ylab = deparse(substitute(y))
  }
  return(attr(x, ylab))
}

#' @describeIn s3_dist_fns Support for auto suggests on `dist_fns_list`s
#' @param pattern a regular expression
#' @returns the names of the attributes
#' @export
.AtNames.dist_fns = function(x, pattern) {
  return(.DollarNames(attributes(x), pattern))
}

#' @describeIn s3_dist_fns Concatenate a `dist_fns` S3 object or `dist_fns_list`s
#' @param ... some of `dist_fns`, `dist_fns_list`s or lists of `dist_fns`
#' @returns a `dist_fns_list`
#' @export
c.dist_fns = function(...) {
  dots = rlang::list2(...)
  if (length(dots) == 0) {
    return(as.dist_fns_list(NULL))
  }
  dots[[1]] = as.dist_fns_list(dots[[1]])
  return(do.call(c.dist_fns_list, dots))
}

#' @describeIn s3_dist_fns Repeat an `dist_fns` S3 object
#' @returns a `dist_fns_list`
#' @export
rep.dist_fns = function(x, ...) {
  rep(as.dist_fns_list(x), ...)
}

#' @describeIn s3_dist_fns Convert a `dist_fns` S3 object into a plain list
#' @export
#' @returns the internal structure of the object as a plain list
#' @concept dist_fns_s3
as.list.dist_fns = function(x, ...) {
  .unclass_dist_fns(x)
}

# S3 dist_fns_list class ----

#' Manipulate `dist_fns` S3 object lists
#'
#' These functions allow generic list behaviour.
#'
#' @param x a `dist_fns_list` S3 object
#' @param ... passed onto methods
#'
#' @concept dist_fns_s3
#' @name s3_dist_fns_list
NULL

#' @describeIn s3_dist_fns_list Unformat the `dist_fns_list`
#' @export
#' @keywords internal
#' @concept dist_fns_s3
#' a plain list of `dist_fns` S3 objects
as.list.dist_fns_list = function(x, ...) {
  unclass(x)
}

#' @describeIn s3_dist_fns_list Length of a `dist_fns_list`
#' @export
#' @keywords internal
#' @returns the length
#' @concept dist_fns_s3
length.dist_fns_list = function(x, ...) {
  return(length(unclass(x)))
}

#' @describeIn s3_dist_fns_list Check is a `dist_fns_list`
#' @export
#' @concept dist_fns_s3
is.dist_fns_list = function(x) {
  inherits(x, "dist_fns_list")
}

# x is a list but not a dist_fns. This means it is a plain list
# or a dist_fns_list
.is_list_excl_dist_fns = function(x) {
  if (!is.list(x)) {
    return(FALSE)
  }
  if (is.dist_fns(x)) {
    return(FALSE)
  }
  return(TRUE)
}

#' Cast to a list of `dist_fns` S3 objects
#'
#' This function wraps `dist_fns` and unwraps plain lists such that
#' the result is a flat `dist_fns_list` containing `dist_fns` objects only
#'
#' @param x a list
#' @return a `dist_fns_list` S3 object
#' @export
#' @concept dist_fns_s3
as.dist_fns_list = function(x) {
  if (length(unlist(x)) == 0) {
    # .class may be asserted when creating zero size / NULL `dist_fns`s
    return(structure(list(), class = c("dist_fns_list", "list")))
  }
  if (is.dist_fns(x)) {
    return(structure(list(x), class = c("dist_fns_list", "list")))
  }
  if (is.list(x)) {
    # x is a list or dist_fns_list (but cannot be a single dist_fns at this point)
    while (any(sapply(x, .is_list_excl_dist_fns))) {
      # if there are any nested dist_fns_lists or plain lists we will collapse
      # them . A dist_fns_list must be a list of dist_fnss without hierarchy.
      # We also have to make sure that plain dist_fnss are wrapped
      x = lapply(x, as.dist_fns_list)
      x = unlist(x, recursive = FALSE)
    }
    return(structure(x, class = c("dist_fns_list", "list")))
  }

  stop(
    "Not convertible to a `dist_fns_list` x is not `dist_fns_list`, a `dist_fns` or a uniform list of `dist_fns`s",
    call. = FALSE
  )
}

#' @describeIn s3_dist_fns_list Format a `dist_fns_list`
#' @export
#' @returns a character vector for the list
#' @keywords internal
#' @concept dist_fns_s3
format.dist_fns_list = function(x, ...) {
  unlist(lapply(x, format))
}

#' @describeIn s3_dist_fns_list Print a `dist_fns_list`
#' @export
#' @returns nothing
#' @keywords internal
#' @concept dist_fns_s3
print.dist_fns_list = function(x, ...) {
  cat(sprintf("dist_fns(%s)\n", length(x)))
  cat(suppressWarnings(format.dist_fns_list(x, ...)), "\n")
  invisible(NULL)
}

#' @describeIn s3_dist_fns_list Convert a `dist_fns_list` to character
#' @export
#' @returns a character vector
#' @keywords internal
#' @concept dist_fns_s3
as.character.dist_fns_list = function(x, ...) {
  format.dist_fns_list(x, ...)
}

#' @exportS3Method pillar::type_sum dist_fns_list
#' @keywords internal
#' @concept dist_fns_s3
type_sum.dist_fns_list = function(x, ...) {
  I(sprintf("<%s[]>", abbreviate("dist_fns", 3, named = FALSE)))
}

#' @exportS3Method pillar::pillar_shaft
#' @keywords internal
#' @concept dist_fns_s3
pillar_shaft.dist_fns_list <- function(x, ...) {
  out <- format.dist_fns_list(x, ...)
  pillar::new_pillar_shaft_simple(out, align = "right")
}

#' @describeIn s3_dist_fns_list Concatenate or construct a `dist_fns_list`
#' @param ... some of `dist_fns_list` and `dist_fns` or list of `dist_fns`s
#' @returns a fully flattened `dist_fns_list` S3 object
#' @export
c.dist_fns_list = function(...) {
  dots = rlang::list2(...)
  if (is.dist_fns_list(dots)) {
    return(dots)
  }
  if (length(dots) == 1) {
    return(as.dist_fns_list(dots))
  }
  # remove empty items
  dots = dots[sapply(dots, length) > 0]
  # make sure all list entries are a dist_fns list (dots is list of dist_fns_list)
  tmp = lapply(dots, as.dist_fns_list)
  # convert to plain list of lists
  tmp = lapply(tmp, as.list.dist_fns_list)
  # collapse one level
  tmp = unlist(tmp, recursive = FALSE)
  # convert to dist_fns_list. this should throw an error if types are mixed.
  return(as.dist_fns_list(tmp))
}

#' @describeIn s3_dist_fns_list Repeat a `dist_fns_list`
#' @keywords internal
#' @export
rep.dist_fns_list = function(x, ...) {
  tmp = NextMethod()
  return(as.dist_fns_list(tmp))
}


#' @describeIn s3_dist_fns_list Repeat a `dist_fns_list`
#' @param decreasing reverse the sort order
#' @keywords internal
#' @export
sort.dist_fns_list = function(x, decreasing = FALSE, ...) {
  indx = order(unlist(sapply(x, .order_dist_fns)))
  if (decreasing) {
    indx = rev(indx)
  }
  return(x[indx])
}

## dist_fns_list Subsetting functions ----

#' @describeIn s3_dist_fns_list Extract named item(s) from a `dist_fns_list`
#' @param y item to retrieve
#' @keywords internal
#' @export
`$.dist_fns_list` = function(x, y) {
  if (is.character(y)) {
    ylab = y
  } else {
    ylab = deparse(substitute(y))
  }
  if (length(x) == 1) {
    return(x[[1]][[ylab]])
  }
  return(sapply(seq_along(x), function(i) x[[i]][[ylab]], USE.NAMES = FALSE))
}

#' Support for auto suggests on `dist_fns_list`s
#' @keywords internal
#' @param x a `dist_fns_list`
#' @keywords internal
#' @returns the names of the children
#' @export
.DollarNames.dist_fns_list = function(x, pattern) {
  if (length(x) == 0) {
    return(character())
  }
  return(.DollarNames(x[[1]], pattern))
}

#' @describeIn s3_dist_fns_list Extract named item(s) from a `dist_fns_list`
#' @param y attribute to retrieve
#' @keywords internal
#' @returns a vector or list of the underlying `dist_fns` attribute values
#' @export
`@.dist_fns_list` = function(x, y) {
  if (is.character(y)) {
    ylab = y
  } else {
    ylab = deparse(substitute(y))
  }
  if (length(x) == 1) {
    return(attr(x[[1]], ylab))
  }
  return(unname(sapply(
    seq_along(x),
    function(i) attr(x[[i]], ylab),
    USE.NAMES = FALSE
  )))
}

#' Support for auto suggests on `dist_fns_list`s
#' @keywords internal
#' @param x a `dist_fns_list`
#' @keywords internal
#' @returns the names of the attributes
#' @export
.AtNames.dist_fns_list = function(x, pattern) {
  if (length(x) == 0) {
    return(character())
  }
  return(.DollarNames(attributes(x[[1]]), pattern))
}

#' @describeIn s3_dist_fns_list Subset a `dist_fns_list`
#' @returns a `dist_fns_list` S3 object
#' @keywords internal
#' @export
`[.dist_fns_list` = function(x, ...) {
  y = `[`(unclass(x), ...)
  return(as.dist_fns_list(y))
}

#' @describeIn s3_dist_fns_list Assign a subset to a `dist_fns_list`
#' @param value the value as `dist_fns_list` or `dist_fns` S3 objects
#' @keywords internal
#' @returns the updated `dist_fns_list` S3 object
#' @export
`[<-.dist_fns_list` = function(x, ..., value) {
  if (!is.dist_fns_list(value) && !is.dist_fns(value)) {
    stop(
      "cannot add a `",
      class(value)[[1]],
      "` to a `dist_fns_list`"
    )
  }
  y = `[<-`(unclass(x), ..., value)
  return(as.dist_fns_list(y))
}

#' @describeIn s3_dist_fns_list get a value from a `dist_fns_list`
#' @returns a `dist_fns` S3 object
#' @keywords internal
#' @export
`[[.dist_fns_list` = function(x, ...) {
  y = `[[`(unclass(x), ...)
  return(y)
}

#' @describeIn s3_dist_fns_list set a single value in a `dist_fns_list`
#' @param value the value
#' @keywords internal
#' @returns the updated `dist_fns_list` S3 object
#' @export
`[[<-.dist_fns_list` = function(x, ..., value) {
  if (is.dist_fns_list(value) && length(value) == 1) {
    value = value[[1]]
  }
  if (!is.dist_fns(value)) {
    stop(
      "cannot add a `",
      class(value)[[1]],
      "` to a `dist_fns_list`"
    )
  }
  y = `[[<-`(unclass(x), ..., value)
  return(as.dist_fns_list(y))
}

#' Apply a function to each element of a vector returning a `dist_fns_list`
#'
#' Analogous to `purrr::map_dbl()`
#'
#' @inheritParams purrr::map
#' @param .f a function to apply that returns a `dist_fns` S3 object (usually an
#'   ``as.dist_fns()`` call)
#' @seealso [purrr::map()]
#'
#' @returns a `dist_fns_list`
#' @export
#' @concept dist_fns_s3
map_dist_fns = function(.x, .f, ..., .progress = FALSE) {
  # This will flatten any nested dist_fns_lists. This is good as .f may return a
  # single dist_fns or more likely a 1 element dist_fns_list.
  return(purrr::map(.x, .f, ..., .progress = .progress) %>% as.dist_fns_list())
}


#' Map over two inputs returning a `dist_fns_list`
#'
#' Analogous to `purrr::map2_dbl()`
#'
#' @inheritParams purrr::map2
#' @param .f a function to apply to each `.x`, `.y` pair that returns a `dist_fns`
#'   S3 object (usually an ``as.dist_fns()`` call)
#' @seealso [purrr::map2()]
#'
#' @returns a `dist_fns_list`
#' @export
#' @concept dist_fns_s3
map2_dist_fns = function(.x, .y, .f, ..., .progress = FALSE) {
  return(
    purrr::map2(.x, .y, .f, ..., .progress = .progress) %>% as.dist_fns_list()
  )
}

#' Map over multiple inputs returning a `dist_fns_list`
#'
#' Analogous to `purrr::pmap_dbl()`
#'
#' @inheritParams purrr::pmap
#' @param .f a function to apply to each `.l` item (usually an ``as.dist_fns()`` call)
#' @seealso [purrr::map()]
#' @returns a `dist_fns_list`
#' @export
#' @concept dist_fns_s3
pmap_dist_fns = function(.l, .f, ..., .progress = FALSE) {
  return(
    purrr::pmap(
      .l,
      .f,
      ...,
      .progress = .progress
    ) %>%
      as.dist_fns_list()
  )
}
