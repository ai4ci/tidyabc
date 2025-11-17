# Sort a `link_fns_list`

These boilerplate functions allow generic list behaviour from
`link_fns_list` classes allowing `link_fns` S3 objects to be used in
lists or dataframes.

## Usage

``` r
# S3 method for class 'link_fns_list'
sort(x, decreasing = FALSE, ...)
```

## Arguments

- x:

  a `link_fns_list` S3 object

- decreasing:

  reverse the sort order

- ...:

  passed onto methods

## Value

an ordered `link_fns_list` S3 object
