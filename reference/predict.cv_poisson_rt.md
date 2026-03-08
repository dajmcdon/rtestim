# Predict observed data using estimated Rt

Given an object of class `poisson_rt` produced with
[`estimate_rt()`](https://dajmcdon.github.io/rtestim/reference/estimate_rt.md),
calculate predicted observed cases for the estimated Rt values. Note:
This function is not intended for "new x" or to produce forecasts, but
rather to examine how Rt relates to observables.

## Usage

``` r
# S3 method for class 'cv_poisson_rt'
predict(object, which_lambda = c("lambda.min", "lambda.1se"), ...)
```

## Arguments

- object:

  result of cross validation of type `cv_poisson_rt`

- which_lambda:

  Select which lambdas from the object to use. If not provided, all Rt's
  are returned. Note that new lambdas not originally used in the
  estimation procedure may be provided, but the results will be
  calculated by linearly interpolating the estimated Rt's.

  The strings `lambda.min` or `lambda.1se` are allowed to choose either
  the lambda that minimizes the cross validation score or the largest
  lambda whose corresponding cross validation score is within 1 standard
  error of the minimal cross validation score.

- ...:

  not used.

## Value

A vector or matrix of predicted case counts.

## Examples

``` r
y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
cv <- cv_estimate_rt(y, korder = 3, nfold = 3, nsol = 30)
p <- predict(cv)
p <- predict(cv, which_lambda = cv$lambda[1])
p <- predict(cv, which_lambda = "lambda.1se")
p <- predict(cv, which_lambda = NULL)
plot(y)
matlines(p, lty = 2)
```
