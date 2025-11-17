# Extract named item(s) from a `link_fns_list`

These boilerplate functions allow generic list behaviour from
`link_fns_list` classes allowing `link_fns` S3 objects to be used in
lists or dataframes.

## Usage

``` r
# S3 method for class 'link_fns_list'
x@y
```

## Arguments

- x:

  a `link_fns_list` S3 object

- y:

  attribute to retrieve

## Value

a vector or list of the underlying `link_fns` attribute values
