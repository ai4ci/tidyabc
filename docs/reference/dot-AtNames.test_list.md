# Support for auto suggests on `test_list`s

These boilerplate functions allow generic list behaviour from
`test_list` classes allowing `test` S3 objects to be used in lists or
dataframes.

## Usage

``` r
# S3 method for class 'test_list'
.AtNames(x, pattern)
```

## Arguments

- x:

  a `test_list` S3 object

## Value

the names of the attributes
