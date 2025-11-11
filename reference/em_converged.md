# Convergence Test for EM-Algorithm

Convergence Test for EM-Algorithm

## Usage

``` r
em_converged(loglik, previous_loglik, tol = 1e-04, check.increased = FALSE)
```

## Arguments

- loglik:

  numeric. Current value of the log-likelihood function.

- previous_loglik:

  numeric. Value of the log-likelihood function at the previous
  iteration.

- tol:

  numerical. The tolerance of the test. If \|LL(t) - LL(t-1)\| / avg \<
  tol, where avg = (\|LL(t)\| + \|LL(t-1)\|)/2, then algorithm has
  converged.

- check.increased:

  logical. Check if likelihood has increased.

## Value

A logical statement indicating whether EM algorithm has converged. if
`check.increased = TRUE`, a vector with 2 elements indicating the
convergence status and whether the likelihood has decreased.

## See also

[dfms-package](https://sebkrantz.github.io/dfms/reference/dfms-package.md)

## Examples

``` r
em_converged(1001, 1000)
#> [1] FALSE
em_converged(10001, 10000)
#> [1] TRUE
em_converged(10001, 10000, check = TRUE)
#> converged  decrease 
#>      TRUE     FALSE 
em_converged(10000, 10001, check = TRUE)
#> converged  decrease 
#>      TRUE      TRUE 
```
