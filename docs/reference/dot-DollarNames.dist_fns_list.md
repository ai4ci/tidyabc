# Support for auto suggests on `dist_fns_list`s

These boilerplate functions allow generic list behaviour from
`dist_fns_list` classes allowing `dist_fns` S3 objects to be used in
lists or dataframes.

## Usage

``` r
# S3 method for class 'dist_fns_list'
.DollarNames(x, pattern)
```

## Arguments

- x:

  a `dist_fns_list` S3 object

## Value

the names of the children
