# Create a `link_fns` S3 object

The link function class allows forwards and backwards transformation.
Link functions can be defined by name or using a statistical
distribution in which case the forward link is a logit of the cumulative
probability and the reverse is the quantile of the expit.

## Arguments

- x:

  a `link_fns` S3 object

- ...:

  passed onto methods

## Details

A `link_fns` S3 object encapsulates a monotonic transformation function
\\h\\, its inverse \\h^{-1}\\, and their derivatives \\h'\\ and
\\(h^{-1})'\\. It also defines the support (domain) \\\[a, b\]\\ of the
original space and the range \\\[h(a), h(b)\]\\ of the transformed
space.

The function dispatches based on the input `x`:

- `character`: Selects standard links (e.g., "log", "logit", "probit",
  "identity"). For example, "log" defines \\h(x) = \log(x)\\ with
  support \\(0, \infty)\\.

- `dist_fns`: Defines the link as the logit of the CDF and the quantile
  of the expit: \\h(x) = \text{logit}(F(x))\\, \\h^{-1}(z) =
  F^{-1}(\text{expit}(z))\\, where \\F\\ and \\F^{-1}\\ are the CDF and
  quantile functions from the `dist_fns` object. The support is
  determined by the quantile function's range (e.g., \\\[Q(0),
  Q(1)\]\\).

- `family` (from `stats`): Uses the link function and its inverse from
  the GLM family object.

- `numeric`: If length 2, interprets as a support range \\\[a, b\]\\ and
  creates a logit-like transformation mapping this range to \\(-\infty,
  \infty)\\: \\h(x) = \text{logit}(\frac{x-a}{b-a})\\. If all values are
  finite and length \> 2, creates a standardization link: \\h(x) =
  \frac{x-\mu}{\sigma}\\, where \\\mu\\ and \\\sigma\\ are the mean and
  standard deviation of the input vector.
