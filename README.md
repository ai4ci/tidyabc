# tidyabc <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/ai4ci/tidyabc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ai4ci/tidyabc/actions/workflows/R-CMD-check.yaml)
[![R-universe version](https://ai4ci.r-universe.dev/tidyabc/badges/version)](https://ai4ci.r-universe.dev/tidyabc)
[![EPSRC badge](https://img.shields.io/badge/EPSRC%20grant-EP%2FY028392%2F1-05acb5)](https://gtr.ukri.org/projects?ref=EP%2FY028392%2F1)
[![DOI](https://zenodo.org/badge/1081063713.svg)](https://doi.org/10.5281/zenodo.17635872)
<!-- badges: end -->

## Approximate Bayesian Computation Made Tidy

`tidyabc` provides a flexible, intuitive, and modern framework for
Approximate Bayesian Computation (ABC). Designed for data-driven
science, it allows you to perform Bayesian inference on complex models where the
likelihood function is intractable or unknowable — common in fields like
epidemiology, ecology, systems biology, and social science.

Unlike traditional ABC tools, `tidyabc` integrates with the
`tidyverse` ecosystem. Define your simulation model and summary statistics using
familiar R functions, leverage powerful `dist_fns` objects to represent priors
and posteriors, and use `dplyr`-style workflows to manage your ABC iterations.
Whether you need a simple rejection sampler or the efficiency of Sequential
Monte Carlo (SMC) or Adaptive ABC, `tidyabc` gives you the tools to build,
diagnose, and interpret your inference within a clear, reproducible, and
publication-ready workflow.

## Installation

`tidyabc` is hosted on the [AI4CI r-universe](https://ai4ci.r-universe.dev/).
Installation from there is as follows:

``` r
options(repos = c(
  "ai4ci" = 'https://ai4ci.r-universe.dev/',
  CRAN = 'https://cloud.r-project.org'))

# Download and install tidyabc in R
install.packages("tidyabc")
```

You can install the development version of `tidyabc` from
[GitHub](https://github.com/ai4ci/tidyabc) with:

``` r
# install.packages("devtools")
devtools::install_github("ai4ci/tidyabc")
```

## Funding

The authors gratefully acknowledge the support of the UK Research and Innovation
AI programme of the Engineering and Physical Sciences Research Council [EPSRC
grant EP/Y028392/1](https://gtr.ukri.org/projects?ref=EP%2FY028392%2F1).
