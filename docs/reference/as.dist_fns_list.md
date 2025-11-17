# Cast to a list of `dist_fns` S3 objects

This function wraps `dist_fns` and unwraps plain lists such that the
result is a flat `dist_fns_list` containing `dist_fns` objects only

## Usage

``` r
as.dist_fns_list(x)
```

## Arguments

- x:

  a `dist_fns_list` S3 object

## Value

a `dist_fns_list` S3 object
