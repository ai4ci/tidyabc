# Generate a distribution from a link transform of another

Generates a new distribution by applying a link transformation to an
existing distribution `dist`. If \\X \sim \text{dist}\\ and \\h\\ is the
link function, this function returns the distribution of \\Y =
h^{-1}(X)\\. The CDF \\F_Y\\, quantile function \\Q_Y\\, PDF \\f_Y\\,
and RNG \\R_Y\\ of the transformed distribution are derived from the
original distribution's functions \\F_X\\, \\Q_X\\, \\f_X\\, \\R_X\\ and
the link function \\h\\ and its inverse \\h^{-1}\\ as follows: \$\$
F_Y(y) = F_X(h(y)) \$\$ \$\$ Q_Y(p) = h^{-1}(Q_X(p)) \$\$ \$\$ f_Y(y) =
f_X(h(y)) \cdot \|h'(y)\| \$\$ \$\$ R_Y(n) = h^{-1}(R_X(n)) \$\$ where
\\h'(y)\\ is the derivative of the link function. This function
implements these transformations for the `p`, `q`, `d`, and `r`
functions of the resulting `dist_fns` object.

## Usage

``` r
transform(link, dist, ..., name = NULL)
```

## Arguments

- link:

  a link function (or name of a link function)

- dist:

  distribution(s) as a name, function or a `dist_fn` S3 object

- ...:

  parameters for the underlying distribution if `dist` is a name or
  function.

- name:

  a name for the link function

## Value

a `dist_fn` or `dist_fn_list` holding the transformed distribution(s)

## Unit tests




    t = transform("log","norm")
    ps = seq(0,1,0.1)
    qs = 0:6
    testthat::expect_equal(t$q(ps),qlnorm(ps))
    testthat::expect_equal(t$p(qs),plnorm(qs))
    testthat::expect_equal(t$d(qs),dlnorm(qs))

    t2 = transform("log","norm", 0.4, 0.1)
    testthat::expect_equal(t2$q(ps),qlnorm(ps,0.4,0.1))

## Examples

``` r
n = as.dist_fns("norm",mean=0.5, sd=0.1)
t = transform("log",n)

plot(t)+ggplot2::geom_function(fun=~ dlnorm(.x, 0.5, 0.1), linetype="dashed")

```
