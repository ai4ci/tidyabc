# set a single value in a `test_list`

These boilerplate functions allow generic list behaviour from
`test_list` classes allowing `test` S3 objects to be used in lists or
dataframes.

## Usage

``` r
# S3 method for class 'test_list'
x[[...]] <- value
```

## Arguments

- x:

  a `test_list` S3 object

- ...:

  passed onto methods

- value:

  the value

## Value

the updated `test_list` S3 object
