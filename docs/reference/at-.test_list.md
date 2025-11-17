# Extract named item(s) from a `test_list`

These boilerplate functions allow generic list behaviour from
`test_list` classes allowing `test` S3 objects to be used in lists or
dataframes.

## Usage

``` r
# S3 method for class 'test_list'
x@y
```

## Arguments

- x:

  a `test_list` S3 object

- y:

  attribute to retrieve

## Value

a vector or list of the underlying `test` attribute values
