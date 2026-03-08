# Add confidence bands to estimated Rt or incidence curves

Create an approximate confidence band for the Rt or incidence estimate.
Note that the variance computation is approximate.

## Usage

``` r
confband(object, lambda, level = 0.95, type = c("Rt", "Yt"), ...)
```

## Arguments

- object:

  a `poisson_rt` or `cv_poisson_rt` object.

- lambda:

  the selected lambda. May be a scalar value, or in the case of
  `cv_poisson_rt` objects, `"lambda.min"` or `"lambda.max"`.

- level:

  the desired confidence level(s). These will be sorted if necessary.

- type:

  the type `Rt` or `Yt` for confidence intervals of fitted Rt or fitted
  incident cases

- ...:

  additional arguments for methods. Unused.

## Value

A `data.frame` containing the estimates `Rt` or `Yt` at the chosen
`lambda`, and confidence limits corresponding to `level`

## Examples

``` r
y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
out <- estimate_rt(y, nsol = 10)
head(confband(out, out$lambda[2]))
#> An `rt_confidence_band` object.
#> 
#> * type = Rt 
#> * lambda = 112.578 
#> * degrees of freedom = 4 
#> 
#> # A tibble: 6 × 3
#>     fit `2.5%` `97.5%`
#>   <dbl>  <dbl>   <dbl>
#> 1 0.849  0        2.31
#> 2 0.894  0        1.94
#> 3 0.939  0        2.03
#> 4 0.983  0        1.98
#> 5 1.03   0.135    1.92
#> 6 1.07   0.156    1.98
head(confband(out, out$lambda[2], level = c(0.95, 0.8, 0.5)))
#> An `rt_confidence_band` object.
#> 
#> * type = Rt 
#> * lambda = 112.578 
#> * degrees of freedom = 4 
#> 
#> # A tibble: 6 × 7
#>     fit `2.5%` `10.0%` `25.0%` `75.0%` `90.0%` `97.5%`
#>   <dbl>  <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#> 1 0.849  0       0       0.352    1.35    1.80    2.31
#> 2 0.894  0       0.216   0.538    1.25    1.57    1.94
#> 3 0.939  0       0.227   0.565    1.31    1.65    2.03
#> 4 0.983  0       0.335   0.643    1.32    1.63    1.98
#> 5 1.03   0.135   0.447   0.722    1.33    1.61    1.92
#> 6 1.07   0.156   0.475   0.757    1.38    1.66    1.98

cv <- cv_estimate_rt(y, nfold = 3, nsol = 30)
head(confband(cv, "lambda.min", c(0.5, 0.9)))
#> An `rt_confidence_band` object.
#> 
#> * type = Rt 
#> * lambda = 165.972 
#> * degrees of freedom = 4 
#> 
#> # A tibble: 6 × 5
#>     fit   `5%` `25%` `75%` `95%`
#>   <dbl>  <dbl> <dbl> <dbl> <dbl>
#> 1 0.843 0      0.354  1.33  2.04
#> 2 0.891 0.0284 0.540  1.24  1.75
#> 3 0.939 0.0403 0.573  1.31  1.84
#> 4 0.987 0.161  0.650  1.32  1.81
#> 5 1.03  0.292  0.731  1.34  1.78
#> 6 1.08  0.324  0.772  1.39  1.84
```
