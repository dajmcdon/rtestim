# Predict observed data using estimated Rt

Given an object of class `poisson_rt` produced with
[`estimate_rt()`](https://dajmcdon.github.io/rtestim/reference/estimate_rt.md),
calculate predicted observed cases for the estimated Rt values. Note:
This function is not intended for "new x" or to produce forecasts, but
rather to examine how Rt relates to observables.

## Usage

``` r
# S3 method for class 'poisson_rt'
predict(object, lambda = NULL, ...)
```

## Arguments

- object:

  An object of class `poisson_rt` produced with
  [`estimate_rt()`](https://dajmcdon.github.io/rtestim/reference/estimate_rt.md).

- lambda:

  Select which lambdas from the object to use. If not provided (the
  default), all are returned. Note that new lambdas not originally used
  in the estimation procedure may be provided, but the results will be
  calculated by linearly interpolating the estimated Rt's.

- ...:

  Not used.

## Value

A vector or matrix of predicted case counts.

## Examples

``` r
y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
out <- estimate_rt(y, nsol = 10)
preds <- predict(out)
plot(y)
matlines(preds, lty = 1)
```
