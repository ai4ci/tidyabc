# Crate a function definition that may not have all the namespaces qualified.

Regular crating requires everything be namespaced and which is easy to
forget, particularly for `stats::` and `dplyr::`.

## Usage

``` r
.autocrate(
  .fn,
  ...,
  .parent_env = baseenv(),
  .error_arg = ".fn",
  .error_call = environment()
)
```

## Arguments

- .fn:

  A fresh formula or function. "Fresh" here means that they should be
  declared in the call to `crate()`. See examples if you need to crate a
  function that is already defined. Formulas are converted to purrr-like
  lambda functions using
  [`rlang::as_function()`](https://rlang.r-lib.org/reference/as_function.html).

- ...:

  Named arguments to declare in the environment of `.fn`.

- .parent_env:

  The default of [`baseenv()`](https://rdrr.io/r/base/environment.html)
  ensures that the evaluation environment of the crate is isolated from
  the search path. Specifying another environment such as the global
  environment allows this condition to be relaxed (but at the expense of
  no longer being able to rely on a local run giving the same results as
  one in a different process).

- .error_arg:

  An argument name as a string. This argument will be mentioned in error
  messages as the input that is at the origin of a problem.

- .error_call:

  The execution environment of a currently running function, e.g.
  `caller_env()`. The function will be mentioned in error messages as
  the source of the error. See the `call` argument of
  [`abort()`](https://rlang.r-lib.org/reference/abort.html) for more
  information.

## Value

a crated function

## Unit tests



    z = 13
    tmp = .autocrate(function(x) {
      y <- runif(x)
      stats::rnorm(x,y)
      x[1] <- 12
      x[2] <- !!z # inlining happens during expression parsing...?
    })

    # unqualified reference is qualified
    testthat::expect_equal(format(body(tmp)), c(
      "{",
      "    y <- stats::runif(x)",
      "    stats::rnorm(x, y)",
      "    x[1] <- 12",
      "    x[2] <- 13",
      "}"
    ))

    tmpfn = function(x) {
      "temp"
    }
    tmpvar = 10

    testthat::expect_error(
      {
        .autocrate(function(y) {
          if (tmpvar == 10) tmpfn(y)
        })
      },
      "References to undefined globals in function `.fn`: tmpvar,tmpfn",
      fixed = TRUE
    )

    with_ref = .autocrate(
      function(y) {
        if (tmpvar == 10) tmpfn(y)
      },
      tmpvar = tmpvar,
      tmpfn=tmpfn
    )

    testthat::expect_equal(
      ls(rlang::fn_env(with_ref)),
      c("tmpfn", "tmpvar")
    )

    library(dplyr)
    # Global definitions in dplyr must be explicitly assigned like in .
    a = .autocrate(function(x) {
      Species = NULL
      x 
    })

    testthat::expect_equal(a(iris), structure(
      c(1L, 1L, 1L, 1L, 1L),
      levels = c("setosa", "versicolor", "virginica"),
      class = "factor"
    ))

    # Autocrating is recursive on functions supplied as globals:
    not_crated_fn = function(x) {rnorm(x); return("success")}

    crated_fn = .autocrate(function(x) {
       generator(x)
    }, generator = not_crated_fn)

    tmp = format(body(rlang::fn_env(crated_fn)$generator))

    testthat::expect_equal(
      tmp,
      c("{", "    stats::rnorm(x)", "    return(\"success\")", "}")
    )
