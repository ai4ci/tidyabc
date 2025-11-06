#' Crate a function that may not have all the namespaces qualified.
#'
#' Regular crating requires everything be namespaced and which is easy to
#' forget, particularly for `stats::` and `dplyr::`.
#'
#' @param expr an expression defining a function
#' @param ... named list of global variables
#'
#' @returns a crated function
#' @keywords internal
#'
#' @unit
#'
#' z = 13
#' tmp = .autocrate(function(x) {
#'   y <- runif(x)
#'   stats::rnorm(x,y)
#'   x[1] <- 12
#'   x[2] <- !!z # inlining happens during expression parsing...?
#' })
#'
#' # unqualified reference is qualified
#' testthat::expect_equal(format(body(tmp)), c(
#'   "{",
#'   "    y <- stats::runif(x)",
#'   "    stats::rnorm(x, y)",
#'   "    x[1] <- 12",
#'   "    x[2] <- 13",
#'   "}"
#' ))
#'
#' tmpfn = function(x) {
#'   "temp"
#' }
#' tmpvar = 10
#'
#' testthat::expect_error(
#'   {
#'     .autocrate(function(y) {
#'       if (tmpvar == 10) tmpfn(y)
#'     })
#'   },
#'   "References to undefined globals in function: tmpvar,tmpfn",
#'   fixed = TRUE
#' )
#'
#' with_ref = .autocrate(
#'   function(y) {
#'     if (tmpvar == 10) tmpfn(y)
#'   },
#'   tmpvar = tmpvar,
#'   tmpfn=tmpfn
#' )
#'
#' testthat::expect_equal(
#'   ls(rlang::fn_env(with_ref)),
#'   c("tmpfn", "tmpvar")
#' )
#'
#' library(dplyr)
#' # Global definitions in dplyr must be explicitly assigned like in .
#' a = .autocrate(function(x) {
#'   Species = NULL
#'   x %>% pull(Species) %>% head(5)
#' })
#'
#' testthat::expect_equal(a(iris), structure(
#'   c(1L, 1L, 1L, 1L, 1L),
#'   levels = c("setosa", "versicolor", "virginica"),
#'   class = "factor"
#' ))
#'
#' # Autocrating is recursive on functions supplied as globals:
#' not_crated_fn = function(x) {rnorm(x); return("success")}
#'
#' crated_fn = .autocrate(function(x) {
#'    generator(x)
#' }, generator = not_crated_fn)
#'
#' tmp = format(body(rlang::fn_env(crated_fn)$generator))
#'
#' testthat::expect_equal(
#'   tmp,
#'   c("{", "    stats::rnorm(x)", "    return(\"success\")", "}")
#' )
#'
#'
.autocrate = function(expr, ...) {
  # splicing happens before the expression is passed to this function
  # so we don;t need to know about it.
  expr = rlang::enexpr(expr)

  # evaluate the function expression:
  fn = eval(expr)
  fn = rlang::as_function(fn)
  .autocrate_fn(fn, ...)
}


.autocrate_fn = function(fn, ...) {
  fmls = formals(fn)
  body = body(fn)

  # named parameters are crated and hence can be used
  supplied = unname(sapply(names(rlang::list2(...)), as.name))
  # parameters passed to the function will be available
  fml_nms = unname(sapply(names(fmls), as.name))

  tmp = .qualify_expression(body, defined = c(supplied, fml_nms))
  # browser()
  # print(tmp$cl)
  if (length(tmp$unqual) > 0) {
    stop(
      "References to undefined globals in function: ",
      paste0(sapply(tmp$unqual, rlang::as_label), collapse = ",")
    )
  }
  body(fn) = tmp$cl
  def = deparse(fn)
  expr2 = parse(text = def, keep.source = FALSE)

  # recursively crate function globals...
  dots = rlang::list2(...)
  dots2 = lapply(dots, function(x) {
    if (is.function(x) && !carrier::is_crate(x)) .autocrate_fn(x) else x
  })
  # browser()
  do.call(
    .shim_crate,
    c(list(expr2), dots2)
  )
}

.shim_crate = function(expr, ...) {
  expr = rlang::enexpr(expr)
  carrier::crate(eval(expr), ...)
}


#' Qualify the namespaces of everything in an expression
#'
#' recursively search for a name that evaluates to a function that is
#' not part of a qualified reference and replaces it with a qualified
#' version. The function will be resolved in the current environment
#' do dplyr::filter will be picked over stats::filter if the tidyverse
#' is loaded at the point of evaluation.
#'
#' @param expr an expression or component of expressions.
#' @param defined a list of qualified references known at this point
#' @param unqual a list of unqualified references
#'
#' @returns a list of the qualified expression, plus updated defined and unqualified lists
#' @keywords internal
#'
#' @unit
#' testthat::expect_equal(
#'   .qualify_expression(expression(12))$cl,
#'   expression(12)
#' )
#'
#' testthat::expect_equal(
#'   .qualify_expression(expression(rnorm))$cl,
#'   expression(stats::rnorm)
#' )
#'
#' tmp = .qualify_expression(expression({new_var = 12}))
#' testthat::expect_equal(tmp$defined, list(as.name("new_var")))
#'
#' tmp2 = .qualify_expression(expression({list(new_var = 12, old_var = 14)}))
#' testthat::expect_equal(tmp2$defined, list())
#'
#' fn = function(x) {
#'   y <- runif(x)
#'   stats::rnorm(x,y)
#'   x[1] <- 12
#'   x[2] <- !!z # inlining happens during expression parsing...?
#' }
#'
#' tmp3 = .qualify_expression(body(fn), defined=list(as.name("x")))
#' testthat::expect_equal(tmp3$unqual, list(as.name("z")))
#' testthat::expect_equal(tmp3$defined, list(as.name("x"),as.name("y")))
#' testthat::expect_equal(format(tmp3$cl), c(
#'   "{",
#'   "    y <- stats::runif(x)",
#'   "    stats::rnorm(x, y)",
#'   "    x[1] <- 12",
#'   "    x[2] <- !!z",
#'   "}"
#' ))
#'
#' library(ggplot2)
#' tmp4 = .qualify_expression(expression(diamonds$cut))
#' testthat::expect_equal(tmp4$cl, expression(ggplot2::diamonds$cut))
#'
#' tmp5 = .qualify_expression(expression(diamonds$imaginary(rnorm(1000))))
#' testthat::expect_equal(
#'   tmp5$cl,
#'   expression(ggplot2::diamonds$imaginary(stats::rnorm(1000)))
#' )
.qualify_expression = function(expr, defined = list(), unqual = list()) {
  # lst = as.list(expr)
  if (is.expression(expr) || is.call(expr)) {
    # input is a complex call or expression
    # even if the expression wraps a length 1 list we need to
    out = list()
    lst = expr
    # browser()
    for (i in seq_along(lst)) {
      cl = lst[[i]]
      if (inherits(cl, "=") || inherits(cl, "<-")) {
        # the call defines an assignation
        assignee = cl[[2]]
        if (is.name(assignee)) {
          # a simple name is being assigned something.
          defined = c(defined, assignee)
        } else {
          # a complex expression. like x[1] = something
          # I don't think we can conclude that this defines the value of the LHS
          # e.g. a[1] = 1 causes an error if `a` does not exist
        }
        # Descend into RHS anyway more to expand the RHS which might need
        # qualification.
        tmp = .qualify_expression(cl[[3]], defined, unqual)
        defined = unique(c(defined, tmp$defined))
        unqual = unique(c(unqual, tmp$unqual))
        cl[[3]] = tmp$cl
        out = c(out, cl)
      } else if (length(cl) > 1) {
        if (cl[[1]] == as.name("::")) {
          # call is a qualified reference to something.
          # we do not need to descend into it.
          out = c(out, cl)
        } else if (cl[[1]] == as.name("$")) {
          # call is a dollar reference to something. the object is unknowable
          # the LHS subject may be a qualifiable item like .data
          # we need to descend into the LHS. The RHS is only interpretable in
          # the context of that list and we can't know (unless it is a function
          # call).
          tmp = .qualify_expression(cl[[2]], defined, unqual)
          defined = unique(c(defined, tmp$defined))
          unqual = unique(c(unqual, tmp$unqual))
          cl[[2]] = tmp$cl
          out = c(out, cl)
        } else {
          # A multi component call. we need to descend into it
          tmp = .qualify_expression(cl, defined, unqual)
          defined = unique(c(defined, tmp$defined))
          unqual = unique(c(unqual, tmp$unqual))
          out = c(out, tmp$cl)
        }
      } else {
        #  a length 1 component e.g. name or numeric
        tmp = .qualify_expression(cl, defined, unqual)
        defined = unique(c(defined, tmp$defined))
        unqual = unique(c(unqual, tmp$unqual))
        out = c(out, tmp$cl)
      }
    }
    out = if (is.expression(expr)) as.expression(out) else as.call(out)
  } else if (is.name(expr)) {
    # input is very simple
    cl = expr

    # an unqualified name:
    if (any(sapply(defined, function(x) cl == x))) {
      # could be reference to global function or variable
      out = cl
    } else if (.is_pkg_item(cl)) {
      # A function in a package or a variable like dplyr::.data.
      out = .qual(cl)
    } else {
      # this is where we end up for references to global
      # functions or variables that have not been correctly imported
      unqual = c(unqual, cl)
      out = cl
    }
  } else {
    # maybe a numeric or other constant?
    out = expr
  }

  names(out) = names(expr)

  return(list(
    cl = out,
    defined = defined,
    unqual = unqual
  ))
}


# internal - does a name resolve to a function or value from a package?
.is_pkg_item = function(name) {
  # is a function
  tmp = try(gsub("package:", "", find(as.character(name))[[1]]), silent = TRUE)
  if (inherits(tmp, "try-error")) {
    return(FALSE)
  }
  return(tmp != ".GlobalEnv")

  # isTRUE(try(is.function(eval(name)), silent = TRUE)) &&
  # resolves to a package function
  # !identical(rlang::fn_env(eval(name)), globalenv())
}

# internal - generate an qualified expression from a function (or variable) name
.qual = function(name, pkg = NULL) {
  name = as.name(name)
  if (is.environment(pkg)) {
    pkg = unname(environmentName(pkg))
  }

  # what is the namespace of the function represented by name?
  # tmp = getNamespaceName(rlang::fn_env(eval(name)))
  # tmp = find(as.character(name))[[1]]
  if (is.null(pkg)) {
    pkg = gsub("package:", "", find(as.character(name))[[1]])
  }
  pkg = as.name(pkg)

  if (pkg == as.name("base") || pkg == as.name("tools:prettycode")) {
    # only don't qualify base functions
    return(name)
  }
  # generate qualified alternative
  as.call(list(as.symbol("::"), pkg, name))
}
