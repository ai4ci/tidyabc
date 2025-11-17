# Apply a function to each element of a vector returning a `test_list`

Analogous to
[`purrr::map_dbl()`](https://purrr.tidyverse.org/reference/map.html)

## Usage

``` r
map_test(.x, .f, ..., .progress = FALSE)
```

## Arguments

- .x:

  A list or atomic vector.

- .f:

  a function to apply that returns a `test` S3 object (usually an
  [`as.test()`](https://ai4ci.github.io/tidyabc/reference/as.test.md)
  call)

- ...:

  Additional arguments passed on to the mapped function.

  We now generally recommend against using `...` to pass additional
  (constant) arguments to `.f`. Instead use a shorthand anonymous
  function:

      # Instead of
      x |> map(f, 1, 2, collapse = ",")
      # do:
      x |> map(\(x) f(x, 1, 2, collapse = ","))

  This makes it easier to understand which arguments belong to which
  function and will tend to yield better error messages.

- .progress:

  Whether to show a progress bar. Use `TRUE` to turn on a basic progress
  bar, use a string to give it a name, or see
  [progress_bars](https://purrr.tidyverse.org/reference/progress_bars.html)
  for more details.

## Value

a `test_list`

## See also

[`purrr::map()`](https://purrr.tidyverse.org/reference/map.html)
