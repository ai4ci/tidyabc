# Crate an existing function that may not have all the namespaces qualified.

Regular crating requires everything be namespaced and which is easy to
forget, particularly for `stats::` and `dplyr::`.

## Usage

``` r
.autocrate_fn(
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
