# Map over multiple inputs returning a `link_fns_list`

Analogous to
[`purrr::pmap_dbl()`](https://purrr.tidyverse.org/reference/pmap.html)

## Usage

``` r
pmap_link_fns(.l, .f, ..., .progress = FALSE)
```

## Arguments

- .l:

  A list of vectors. The length of `.l` determines the number of
  arguments that `.f` will be called with. Arguments will be supply by
  position if unnamed, and by name if named.

  Vectors of length 1 will be recycled to any length; all other elements
  must be have the same length.

  A data frame is an important special case of `.l`. It will cause `.f`
  to be called once for each row.

- .f:

  a function to apply to each `.l` item (usually an
  [`as.link_fns()`](https://ai4ci.github.io/tidyabc/reference/as.link_fns.md)
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

a `link_fns_list`

## See also

[`purrr::map()`](https://purrr.tidyverse.org/reference/map.html)
