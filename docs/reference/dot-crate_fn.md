# Crate a pre-existing function with or without checks

Crate a pre-existing function with or without checks

## Usage

``` r
.crate_fn(
  .fn,
  ...,
  .parent_env = baseenv(),
  .error_arg = ".fn",
  .error_call = environment(),
  .nsqualify = FALSE
)
```

## Arguments

- .fn:

  a reference to a function

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

- .nsqualify:

  check for undefined global references and qualify namespaces?

## Value

a crated function

## Unit tests


    .crate_fn
