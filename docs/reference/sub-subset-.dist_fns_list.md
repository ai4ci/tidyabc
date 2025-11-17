# set a single value in a `dist_fns_list`

These boilerplate functions allow generic list behaviour from
`dist_fns_list` classes allowing `dist_fns` S3 objects to be used in
lists or dataframes.

## Usage

``` r
# S3 method for class 'dist_fns_list'
x[[...]] <- value
```

## Arguments

- x:

  a `dist_fns_list` S3 object

- ...:

  passed onto methods

- value:

  the value

## Value

the updated `dist_fns_list` S3 object
