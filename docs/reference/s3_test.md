# Create a `test` S3 object

TODO: document the purpose of this class, and its internal structure.
including any functions attached to it.

## Arguments

- x:

  a `test` S3 object

- ...:

  passed onto methods

## Details

`test` and `test_list` objects support `$` access for fields and `@`
access for attributes. `test_list`s can be made with the
[`c()`](https://rdrr.io/r/base/c.html) or
[`rep()`](https://rdrr.io/r/base/rep.html) functions, or with the
`purrr` style map functions, and they support subsetting. Individual
`test` members of lists can be accessed with `[[`.
