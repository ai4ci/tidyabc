# Example generators

These are a set of internally cached functions to support examples. They
are exported as internal functions so that the examples can run
correctly and cache their output to prevent excessive repetition of the
code examples.

## Usage

``` r
example_sim_fn(mean, sd1, sd2)

example_scorer_fn(simdata, obsdata)

example_obs()

example_truth()

example_obsdata()

example_priors_list()

example_rejection_fit()

example_smc_fit()

example_adaptive_fit()
```

## Arguments

- mean:

  the example mean

- sd1:

  the example normal sd (A)

- sd2:

  the example gamma sd (B)

- simdata:

  output of example_sim_fn()

- obsdata:

  e.g. example_obsdata()

## Value

a list of A and B of samples from a normal and gamma

a list of scores

output of
[`test_simulation()`](https://ai4ci.github.io/tidyabc/reference/test_simulation.md)

example `truth` parameter

example `obsdata` parameter

example `priors_list` parameter

example `abc_rejection` output

example `abc_smc` output

example `abc_adaptive` output

## Functions

- `example_sim_fn()`: Example simulation function

- `example_scorer_fn()`: Example simulation function

- `example_obs()`: Example output of
  [`test_simulation()`](https://ai4ci.github.io/tidyabc/reference/test_simulation.md)

- `example_truth()`: Example for `truth` parameter

- `example_obsdata()`: Example for `obsdata` parameter

- `example_priors_list()`: Example for `priors_list` parameter

- `example_rejection_fit()`: Example `abc_rejection` output

- `example_smc_fit()`: Example `abc_smc` output

- `example_adaptive_fit()`: Example `abc_adaptive` output

## Examples

``` r
example_sim_fn(3,2,1) %>% lapply(head,10)
#> $A
#>  [1] 2.3869939 6.2155504 4.0387713 3.9976387 5.1475440 0.8287624 1.3667290
#>  [8] 0.5171063 6.4018627 0.8649127
#> 
#> $B
#>  [1] 2.681213 3.487453 4.559011 2.044231 2.607585 6.032914 3.245918 4.874561
#>  [9] 3.193990 3.650157
#> 
example_scorer_fn(
  example_obsdata(),
  example_sim_fn(3,2,1)
)
#> $A
#> [1] 1.042563
#> 
#> $B
#> [1] 1.966807
#> 
example_obs() %>% rapply(head,n=10, how="replace")
#> $obsdata
#> $obsdata$A
#>  [1] 3.879049 4.539645 8.117417 5.141017 5.258575 8.430130 5.921832 2.469878
#>  [9] 3.626294 4.108676
#> 
#> $obsdata$B
#>  [1] 3.963790 4.935891 5.197293 3.450608 5.984255 5.150334 7.583826 5.601786
#>  [9] 4.467521 4.746891
#> 
#> 
#> $truth
#> mean  sd1  sd2 
#>    5    2    1 
#> 
#> $sim_fn
#> [1] "function (mean, sd1, sd2) "                                                    
#> [2] "{"                                                                             
#> [3] "    return(list(A = stats::rnorm(1000, mean, sd1), B = tidyabc::rgamma2(1000, "
#> [4] "        mean, sd2)))"                                                          
#> [5] "}"                                                                             
#> 
#> $obsscores
#> $obsscores$A
#> [1] 0
#> 
#> $obsscores$B
#> [1] 0
#> 
#> 
#> $scorer_fn
#> [1] "function (simdata, obsdata) "                                      
#> [2] "{"                                                                 
#> [3] "    return(list(A = tidyabc::calculate_wasserstein(obsdata$A, "    
#> [4] "        simdata$A), B = tidyabc::calculate_wasserstein(obsdata$B, "
#> [5] "        simdata$B)))"                                              
#> [6] "}"                                                                 
#> 
example_truth()
#> mean  sd1  sd2 
#>    5    2    1 
example_obsdata() %>% lapply(head,10)
#> $A
#>  [1] 3.879049 4.539645 8.117417 5.141017 5.258575 8.430130 5.921832 2.469878
#>  [9] 3.626294 4.108676
#> 
#> $B
#>  [1] 3.963790 4.935891 5.197293 3.450608 5.984255 5.150334 7.583826 5.601786
#>  [9] 4.467521 4.746891
#> 
example_priors_list()
#> Parameters: 
#> * mean: unif(min = 0, max = 10)
#> * sd1: unif(min = 0, max = 5)
#> * sd2: unif(min = 0, max = 5)
#> Constraints:
#> * mean > sd2
example_rejection_fit()
#> ABC rejection, 1 wave.
#> ABC rejection fit: single wave
#> Parameter estimates:
#> # A tibble: 3 × 4
#> # Groups:   param [3]
#>   param mean_sd       median_95_CrI           ESS
#>   <chr> <chr>         <chr>                 <dbl>
#> 1 mean  4.992 ± 0.172 5.004 [4.587 — 5.362]  73.1
#> 2 sd1   2.236 ± 0.515 2.168 [1.367 — 3.804]  73.1
#> 3 sd2   1.230 ± 0.289 1.192 [0.754 — 1.923]  73.1
example_smc_fit()
#> ABC SMC fit: 8 waves - (converged)
#> Parameter estimates:
#> # A tibble: 3 × 4
#> # Groups:   param [3]
#>   param mean_sd       median_95_CrI           ESS
#>   <chr> <chr>         <chr>                 <dbl>
#> 1 mean  4.973 ± 0.034 4.974 [4.892 — 5.041]  249.
#> 2 sd1   1.995 ± 0.071 1.992 [1.802 — 2.225]  249.
#> 3 sd2   0.989 ± 0.035 0.988 [0.900 — 1.080]  249.
example_adaptive_fit()
#> ABC-Adaptive
#> Adaptive waves:  ■■■■■■■■                          24% | wave 2 ETA:  4s
#> Adaptive waves:  ■■■■■■■■■■■■■■■■■■■■              62% | wave 5 ETA:  2s
#> Effective sample size has reduced below 200.
#> Attempting recovery with larger acceptance rate: 33.333%
#> Adaptive waves:  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | wave 8 ETA:  0s
#> ABC adaptive fit: 8 waves - (not yet converged)
#> Parameter estimates:
#> # A tibble: 3 × 4
#> # Groups:   param [3]
#>   param mean_sd       median_95_CrI           ESS
#>   <chr> <chr>         <chr>                 <dbl>
#> 1 mean  4.980 ± 0.036 4.973 [4.894 — 5.076]  263.
#> 2 sd1   2.091 ± 0.064 2.094 [1.769 — 2.215]  263.
#> 3 sd2   1.015 ± 0.034 1.012 [0.925 — 1.167]  263.
```
