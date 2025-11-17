# Map over two inputs returning a `dist_fns_list`

Analogous to
[`purrr::map2_dbl()`](https://purrr.tidyverse.org/reference/map2.html)

## Usage

``` r
map2_dist_fns(.x, .y, .f, ..., .progress = FALSE)
```

## Arguments

- .x, .y:

  A pair of vectors, usually the same length. If not, a vector of length
  1 will be recycled to the length of the other.

- .f:

  a function to apply to each `.x`, `.y` pair that returns a `dist_fns`
  S3 object (usually an
  [`as.dist_fns()`](https://ai4ci.github.io/tidyabc/reference/as.dist_fns.md)
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

a `dist_fns_list`

## See also

[`purrr::map2()`](https://purrr.tidyverse.org/reference/map2.html)
