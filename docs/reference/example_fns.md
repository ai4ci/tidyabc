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
#>  [1] 3.9187221 3.8085454 4.5898907 1.0419118 3.6491363 1.1034679 3.2735285
#>  [8] 0.9694824 4.2431560 1.8132811
#> 
#> $B
#>  [1] 3.747569 2.396480 2.231745 2.966214 2.772659 1.675842 3.028974 1.277203
#>  [9] 1.203135 3.125411
#> 
example_scorer_fn(
  example_obsdata(),
  example_sim_fn(3,2,1)
)
#> $A
#> [1] 0.9916786
#> 
#> $B
#> [1] 2.04281
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
#> 1 mean  4.967 ± 0.180 4.943 [4.592 — 5.448]  76.3
#> 2 sd1   2.348 ± 0.578 2.265 [1.377 — 3.850]  76.3
#> 3 sd2   1.176 ± 0.316 1.094 [0.666 — 1.976]  76.3
example_smc_fit()
#> ABC SMC fit: 8 waves - (converged)
#> Parameter estimates:
#> # A tibble: 3 × 4
#> # Groups:   param [3]
#>   param mean_sd       median_95_CrI           ESS
#>   <chr> <chr>         <chr>                 <dbl>
#> 1 mean  4.977 ± 0.034 4.977 [4.896 — 5.076]  249.
#> 2 sd1   1.982 ± 0.077 1.987 [1.773 — 2.210]  249.
#> 3 sd2   0.994 ± 0.038 0.992 [0.901 — 1.097]  249.
example_adaptive_fit()
#> ABC-Adaptive
#> Adaptive waves:  ■■■■■■■■                          24% | wave 2 ETA:  4s
#> Adaptive waves:  ■■■■■■■■■■■■                      36% | wave 3 ETA:  3s
#> Converged on wave: 8
#> Adaptive waves:  ■■■■■■■■■■■■■■■■■■■■■■■■■■        84% | wave 7 ETA:  1s
#> ABC adaptive fit: 8 waves - (converged)
#> Parameter estimates:
#> # A tibble: 3 × 4
#> # Groups:   param [3]
#>   param mean_sd       median_95_CrI           ESS
#>   <chr> <chr>         <chr>                 <dbl>
#> 1 mean  4.974 ± 0.031 4.969 [4.888 — 5.074]  363.
#> 2 sd1   2.065 ± 0.057 2.075 [1.779 — 2.214]  363.
#> 3 sd2   1.009 ± 0.037 1.011 [0.892 — 1.118]  363.
```
