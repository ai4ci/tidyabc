# Cast to a list of `link_fns` S3 objects

This function wraps `link_fns` and unwraps plain lists such that the
result is a flat `link_fns_list` containing `link_fns` objects only

## Usage

``` r
as.link_fns_list(x)
```

## Arguments

- x:

  a `link_fns_list` S3 object

## Value

a `link_fns_list` S3 object
