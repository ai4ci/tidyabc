# Merge aesthetic mappings and deduplicate the result

Merge aesthetic mappings and deduplicate the result

## Usage

``` r
.gg_merge_aes(..., allow_override = TRUE)
```

## Arguments

- ...:

  a set of `name=value`, and `aes(...)` specifications

- allow_override:

  are the original mappings kept and `name=value` pairs act as default
  and may be overridden by user supplied aesthetics

## Value

a single deduplicated set.

## Unit tests


    m1 = ggplot2::aes(x=a,y=b,colour=class)
    .gg_merge_aes(x=A,y=BB,m1)
