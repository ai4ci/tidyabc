# Create a new `link_fns` S3 object

Default is an identity link

## Usage

``` r
new_link_fns(
  trans = ~.x,
  inv = ~.x,
  support = c(-Inf, Inf),
  ddxtrans = NULL,
  ddxinv = NULL,
  name = NULL
)
```

## Arguments

- trans:

  the forward transforming function

- inv:

  the reverse inverting function

- support:

  the support or support of the link

- range:

  the range of the link (typically `-Inf` to `Inf`)

## Value

a new `link_fns` S3 object
