# Crating functions ----

#' Crate a function definition that may not have all the namespaces qualified.
#'
#' Regular crating requires everything be namespaced and which is easy to
#' forget, particularly for `stats::` and `dplyr::`.
#'
#' @inheritParams carrier::crate
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
#'   "References to undefined globals in function `.fn`: tmpvar,tmpfn",
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
.autocrate = function(
  .fn,
  ...,
  .parent_env = baseenv(),
  .error_arg = ".fn",
  .error_call = environment()
) {
  expr = rlang::enexpr(.fn)
  # splicing happens here so we don;t need to know about it.

  # evaluate the function expression:
  fn = eval(expr)
  fn = try(rlang::as_function(fn), silent = TRUE)
  if (!rlang::is_function(fn)) {
    rlang::abort(
      sprintf("`%s` must evaluate to a function", .error_arg),
      call = .error_call
    )
  }

  .autocrate_fn(
    fn,
    ...,
    .parent_env = .parent_env,
    .error_arg = .error_arg,
    .error_call = .error_call
  )
}


#' Crate an existing function that may not have all the namespaces qualified.
#'
#' Regular crating requires everything be namespaced and which is easy to
#' forget, particularly for `stats::` and `dplyr::`.
#'
#' @inheritParams carrier::crate
#' @returns a crated function
#' @keywords internal
.autocrate_fn = function(
  .fn,
  ...,
  .parent_env = baseenv(),
  .error_arg = ".fn",
  .error_call = environment()
) {
  if (carrier::is_crate(.fn)) {
    return(.fn)
  }

  env <- new.env(parent = .parent_env)
  dots <- rlang::list2(...)
  if (!all(nzchar(rlang::names2(dots)))) {
    rlang::abort("All `...` arguments must be named")
  }
  dots2 = lapply(dots, function(x) {
    if (is.function(x) && !carrier::is_crate(x)) {
      .autocrate_fn(
        x,
        .parent_env = .parent_env,
        .error_arg = .error_arg,
        .error_call = .error_call
      )
    } else {
      x
    }
  })

  for (nm in names(dots)) {
    assign(nm, dots2[[nm]], envir = env)
  }

  fmls = formals(.fn)
  body = body(.fn)

  # named parameters are crated and hence can be used
  supplied = unname(sapply(names(dots), as.name))
  # parameters passed to the function will be available
  fml_nms = unname(sapply(names(fmls), as.name))

  tmp = .qualify_expression(body, defined = c(supplied, fml_nms))
  # browser()
  # print(tmp$cl)
  if (length(tmp$unqual) > 0) {
    rlang::abort(
      sprintf(
        "References to undefined globals in function `%s`: %s",
        .error_arg,
        paste0(sapply(tmp$unqual, rlang::as_label), collapse = ",")
      ),
      call = .error_call
    )
  }

  out = as.function(c(fmls, tmp$cl), envir = env)
  attr(out, "class") = "crate"

  return(out)
}

#' Crate a pre-existing function with or without checks
#'
#' @param .fn a reference to a function
#' @inheritParams carrier::crate
#' @param .nsqualify check for undefined global references and qualify namespaces?
#'
#' @returns a crated function
#' @keywords internal
#'
#' @unit
#' .crate_fn
.crate_fn = function(
  .fn,
  ...,
  .parent_env = baseenv(),
  .error_arg = ".fn",
  .error_call = environment(),
  .nsqualify = FALSE
) {
  if (.nsqualify) {
    .autocrate_fn(
      .fn,
      ...,
      ,
      .parent_env = .parent_env,
      .error_arg = .error_arg,
      .error_call = .error_call
    )
  } else {
    carrier::crate(
      rlang::set_env(.fn),
      ...,
      .parent_env = .parent_env,
      .error_arg = .error_arg,
      .error_call = .error_call
    )
  }
}

#' Create an environment of shared crates.
#'
#' Crated functions that share the same data in a single dependency free
#' environment.
#'
#' @param .fns a named list of functions
#' @inheritParams carrier::crate
#' @param .nsqualify check for undefined global references and qualify namespaces?
#'
#' @returns an environment containing crated functions
#' @keywords internal
#' @unit
#' msg = "hello"
#' e = .super_crate(
#'   list(
#'     plusx = function(z) {
#'      return(z+x)
#'     },
#'     plusy = function(z) {
#'      return(z+y+plusx(z))
#'     }
#'     # splicing does not work
#'     # hi = function() {
#'     #   print(!!msg)
#'     # }
#'   ),
#'   x = 1:10,
#'   y = 11:20
#' )
#'
#' testthat::expect_equal(
#'   e$plusy(1),
#'   c(14, 16, 18, 20, 22, 24, 26, 28, 30, 32)
#' )
.super_crate = function(
  .fns,
  ...,
  .parent_env = baseenv(),
  .error_arg = ".fn",
  .error_call = environment(),
  .nsqualify = FALSE
) {
  env <- new.env(parent = .parent_env)

  dots <- rlang::list2(...)

  if (!all(nzchar(rlang::names2(dots)))) {
    stop("All `...` arguments must be named")
  }

  for (nm in names(dots)) {
    x = dots[[nm]]
    if (is.function(x) && !carrier::is_crate(x)) {
      x = .crate_fn(
        x,
        .parent_env = .parent_env,
        .error_arg = .error_arg,
        .error_call = .error_call,
        .nsqualify = .nsqualify
      )
    }
    assign(nm, x, envir = env)
  }

  # all functions are available and
  # named parameters are crated and hence can be used
  supplied = c(
    unname(sapply(names(.fns), as.name)),
    unname(sapply(names(dots), as.name))
  )

  for (nm in names(.fns)) {
    fn = eval(.fns[[nm]])
    fmls = formals(fn)
    body = body(fn)

    if (.nsqualify) {
      # parameters passed to the function will be available

      fml_nms = unname(sapply(names(fmls), as.name))
      tmp = .qualify_expression(body, defined = c(supplied, fml_nms))

      if (length(tmp$unqual) > 0) {
        stop(
          "References to undefined globals in function: ",
          paste0(sapply(tmp$unqual, rlang::as_label), collapse = ",")
        )
      }

      body = tmp$cl
    }

    assign(nm, as.function(c(fmls, body), envir = env), envir = env)
    attr(env[[nm]], "class") = "crate"
  }
  return(env)
}

# Navigation utils ----

#' Explore crated data for a `crate`
#'
#' @param x a crated function
#' @param y item to retrieve
#' @returns crated data from `x`
#' @export
#' @concept crate_utils
#' @keywords internal
#' @name at.crate
`$.crate` = function(x, y) {
  if (is.character(y)) {
    ylab = y
  } else {
    ylab = deparse(substitute(y))
  }
  return(rlang::fn_env(x)[[ylab]])
}

#' Support for auto suggests on `crate`s
#'
#' @param x a crated function
#' @param pattern a regular expression
#' @returns the names of crated data items
#' @exportS3Method utils::.DollarNames crate
#' @concept crate_utils
#' @keywords internal
.DollarNames.crate = function(x, pattern) {
  return(utils::.DollarNames(rlang::fn_env(x), pattern))
}


# Namespace qualification ----

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
#'
#' tmp6 = .qualify_expression(expression(function(n,mu) rnorm(n,mu)))
#' testthat::expect_equal(tmp6$defined, list())
#' testthat::expect_equal(
#'   format(tmp6$cl),
#'   "expression(function(n, mu) stats::rnorm(n, mu))"
#' )
#'
#' tmp7 = .qualify_expression(expression(Species = NULL))
#' testthat::expect_equal(length(tmp7$cl), 1L)
#' testthat::expect_equal(tmp7$cl, expression(Species = NULL))
#'
.qualify_expression = function(expr, defined = list(), unqual = list()) {
  # TODO: this is very slow.
  # TODO: is it possible to do this with scrRefs attached?
  # lst = as.list(expr)
  if (is.expression(expr) || is.call(expr)) {
    # input is a complex call or expression
    # even if the expression wraps a length 1 list we need to
    out = list()
    lst = expr
    # browser()
    for (i in seq_along(lst)) {
      cl = lst[[i]]
      if (missing(cl)) {
        out = c(out, dplyr::expr())
      } else if (inherits(cl, "=") || inherits(cl, "<-")) {
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
        } else if (cl[[1]] == as.name("function")) {
          # A function definition
          # pick up names from function formals:
          localdefined = c(defined, lapply(names(cl[[2]]), as.name))
          # descend into body:
          tmp = .qualify_expression(cl[[3]], localdefined, unqual)
          # don't inherit defined variables in function
          unqual = unique(c(unqual, tmp$unqual))
          cl[[3]] = tmp$cl
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
    if (class(expr) == "NULL") {
      out = list(expr)
    } else {
      out = expr
    }
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
  tmp = try(
    gsub("package:", "", utils::find(as.character(name))[[1]]),
    silent = TRUE
  )
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
    pkg = gsub("package:", "", utils::find(as.character(name))[[1]])
  }
  pkg = as.name(pkg)

  if (pkg == as.name("base") || pkg == as.name("tools:prettycode")) {
    # only don't qualify base functions
    return(name)
  }
  # generate qualified alternative
  as.call(list(as.symbol("::"), pkg, name))
}
