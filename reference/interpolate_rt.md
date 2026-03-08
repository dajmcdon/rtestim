# Interpolate (or extrapolate) Rt estimates to intermediate design points

Interpolate (or extrapolate) Rt estimates to intermediate design points

## Usage

``` r
interpolate_rt(object, xout, ...)

# S3 method for class 'cv_poisson_rt'
interpolate_rt(object, xout, which_lambda = c("lambda.min", "lambda.1se"), ...)

# S3 method for class 'poisson_rt'
interpolate_rt(object, xout, lambda = NULL, ...)
```

## Arguments

- object:

  A fitted object produced by
  [`estimate_rt()`](https://dajmcdon.github.io/rtestim/reference/estimate_rt.md)
  or
  [`cv_estimate_rt()`](https://dajmcdon.github.io/rtestim/reference/cv_estimate_rt.md).

- xout:

  a vector of new positions at which Rt should be produced, but where
  counts may not have been observed.

- ...:

  additional arguments passed to methods.

- which_lambda:

  Select which lambdas from the object to use. If not provided, all Rt's
  are returned. Note that new lambdas not originally used in the
  estimation procedure may be provided, but the results will be
  calculated by linearly interpolating the estimated Rt's.

  The strings `lambda.min` or `lambda.1se` are allowed to choose either
  the lambda that minimizes the cross validation score or the largest
  lambda whose corresponding cross validation score is within 1 standard
  error of the minimal cross validation score.

- lambda:

  Vector. A user supplied sequence of tuning parameters which determines
  the balance between data fidelity and smoothness of the estimated Rt;
  larger `lambda` results in a smoother estimate. The default, `NULL`
  results in an automatic computation based on `nlambda`, the largest
  value of `lambda` that would result in a maximally smooth estimate,
  and `lambda_min_ratio`. Supplying a value of `lambda` overrides this
  behaviour. It is likely better to supply a decreasing sequence of
  `lambda` values than a single (small) value. If supplied, the
  user-defined `lambda` sequence is automatically sorted in decreasing
  order.

## Value

A vector or matrix of interpolated Rt estimates.

## Examples

``` r
y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
out <- estimate_rt(y)

# originally estimated at
out$x
#>   [1]   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
#>  [19]  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36
#>  [37]  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54
#>  [55]  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72
#>  [73]  73  74  75  76  77  78  79  80  81  82  83  84  85  86  87  88  89  90
#>  [91]  91  92  93  94  95  96  97  98  99 100 101

# get the Rt at 3 new points (for all estimated lambdas)
int <- interpolate_rt(out, c(10.5, 11.5, 12.5))

# get the Rt at a single value of lambda
interpolate_rt(out, c(10.5, 11.5, 12.5), lambda = out$lambda[20])
#> [1] 1.255998 1.277763 1.299562

y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
out <- estimate_rt(y, nsol = 10)
interpolate_rt(out, xout = c(1.5, 2.5))
#>           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#> [1,] 0.9658908 0.9408460 0.9254087 0.9106554 0.8459683 0.6989246 0.5736986
#> [2,] 1.0028122 0.9939503 0.9863701 0.9767920 0.9267339 0.8335073 0.7535534
#>           [,8]      [,9]     [,10]
#> [1,] 0.4744588 0.5088818 0.5234261
#> [2,] 0.6796563 0.5025401 0.3769817
```
