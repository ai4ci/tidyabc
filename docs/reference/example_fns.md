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
#>  [1] -0.2242239  5.1467233  2.4749420  1.9289318  0.7444798  4.6197233
#>  [7]  2.1224445 -1.9099227  3.4123414  2.9349354
#> 
#> $B
#>  [1] 2.495385 5.192140 4.253605 3.195637 1.943004 2.109413 3.193004 6.963725
#>  [9] 5.179482 3.322220
#> 
example_scorer_fn(
  example_obsdata(),
  example_sim_fn(3,2,1)
)
#> $A
#> [1] 0.9791697
#> 
#> $B
#> [1] 1.9886
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
#> 1 mean  5.000 ± 0.176 4.991 [4.593 — 5.417]  74.0
#> 2 sd1   2.309 ± 0.559 2.200 [1.362 — 3.848]  74.0
#> 3 sd2   1.219 ± 0.289 1.223 [0.664 — 2.054]  74.0
example_smc_fit()
#> ABC SMC fit: 8 waves - (converged)
#> Parameter estimates:
#> # A tibble: 3 × 4
#> # Groups:   param [3]
#>   param mean_sd       median_95_CrI           ESS
#>   <chr> <chr>         <chr>                 <dbl>
#> 1 mean  4.977 ± 0.033 4.976 [4.890 — 5.065]  255.
#> 2 sd1   1.993 ± 0.072 1.995 [1.824 — 2.161]  255.
#> 3 sd2   0.986 ± 0.034 0.984 [0.899 — 1.081]  255.
example_adaptive_fit()
#> ABC-Adaptive
#> Effective sample size has reduced below 200.
#> Attempting recovery with larger acceptance rate: 62.500%
#> Adaptive waves:  ■■■■■■■■                          23% | wave 2 ETA:  4s
#> Adaptive waves:  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■      89% | wave 7 ETA:  1s
#> ABC adaptive fit: 8 waves - (not yet converged)
#> Parameter estimates:
#> # A tibble: 3 × 4
#> # Groups:   param [3]
#>   param mean_sd       median_95_CrI           ESS
#>   <chr> <chr>         <chr>                 <dbl>
#> 1 mean  4.964 ± 0.129 4.958 [4.717 — 5.261] 1519.
#> 2 sd1   2.124 ± 0.233 2.078 [1.522 — 2.905] 1519.
#> 3 sd2   1.139 ± 0.173 1.134 [0.753 — 1.531] 1519.
```
