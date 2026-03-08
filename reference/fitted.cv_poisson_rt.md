# Fitted cv_poisson_rt

Fitted cv_poisson_rt

## Usage

``` r
# S3 method for class 'cv_poisson_rt'
fitted(object, which_lambda = c("lambda.min", "lambda.1se"), ...)
```

## Arguments

- object:

  result of cross validation of type `cv_poisson_rt`

- which_lambda:

  select which Rt's to output. If not provided, all Rt's are returned.
  If provided a list of lambda,the corresponding Rt estimation will be
  returned.

  If provided a string, it must be either one of `lambda.min` or
  `lambda.1se`.

  - If provided `lambda.min`, return Rt which is generated from the
    lambda that minimizes the cross validation score.

  - If provided `lambda.1se`, return Rt which is generated from the
    lambda whose corresponding cross validation score is 1 standard
    error away of the minimal cross validation score.

- ...:

  not used.

## Value

Rt's estimated from provided lambda

## Examples

``` r
y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
cv <- cv_estimate_rt(y, korder = 3, nfold = 3, nsol = 30)
f <- fitted(cv)
f <- fitted(cv, which_lambda = cv$lambda[1])
f <- fitted(cv, which_lambda = "lambda.1se")
f <- fitted(cv, which_lambda = NULL)
```
