# Assign a subset to a `link_fns_list`

These boilerplate functions allow generic list behaviour from
`link_fns_list` classes allowing `link_fns` S3 objects to be used in
lists or dataframes.

## Usage

``` r
# S3 method for class 'link_fns_list'
x[...] <- value
```

## Arguments

- x:

  a `link_fns_list` S3 object

- ...:

  passed onto methods

- value:

  the value as `link_fns_list` or `link_fns` S3 objects

## Value

the updated `link_fns_list` S3 object
