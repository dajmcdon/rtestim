# Rt estimation algorithm configuration

Rt estimation algorithm configuration

## Usage

``` r
configure_rt_admm(
  rho = -1,
  alpha = 0.5,
  gamma = 0.9,
  tolerance = 1e-04,
  maxiter_newton = 50L,
  maxiter_line = 20L,
  verbose = 0,
  ...
)
```

## Arguments

- rho:

  Double. An ADMM parameter; coefficient of augmented term in the
  Lagrangian function.

- alpha:

  Double. A parameter adjusting upper bound in line search algorithm in
  `prox_newton` algorithm.

- gamma:

  Double. A parameter adjusting step size in line search algorithm in
  `prox_newton` algorithm.

- tolerance:

  Double. Tolerance of ADMM convergence.

- maxiter_newton:

  Integer. Maximum number of iterations for the outer Newton iteration.

- maxiter_line:

  Integer. Maximum number of iterations for the linesearch algorithm in
  the proximal Newton method.

- verbose:

  Integer.

- ...:

  space for future extensions

## Value

a list of model parameters with class `rt_admm_configuration`

## Examples

``` r
configure_rt_admm()
#> 
#> ── An Rt ADMM configuration ──
#> 
#> rho: -1
#> alpha: 0.5
#> gamma: 0.9
#> tolerance: 1e-04
#> maxiter_newton: 50
#> maxiter_line: 20
#> verbose: 0
configure_rt_admm(tolerance = 1e-6, verbose = 1L)
#> 
#> ── An Rt ADMM configuration ──
#> 
#> rho: -1
#> alpha: 0.5
#> gamma: 0.9
#> tolerance: 1e-06
#> maxiter_newton: 50
#> maxiter_line: 20
#> verbose: 1
```
