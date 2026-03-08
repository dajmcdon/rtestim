# Plot cv_poisson_rt

Plot cv_poisson_rt

## Usage

``` r
# S3 method for class 'cv_poisson_rt'
plot(x, which_lambda = c("cv_scores", "lambda.min", "lambda.1se"), ...)
```

## Arguments

- x:

  result of cv_estimate_rt of class `cv_poisson_rt`

- which_lambda:

  select which Rt's to plot.

  If not provided, the cross validation score will be plotted. If
  provided a list of lambda, the corresponding Rt estimation will be
  plotted.

  If provided a string, it must be either one of `lambda.min`,
  `lambda.1se`, or `cv_scores`.

  - If provided `lambda.min`, plot Rt which is generated from the lambda
    that minimizes the cross validation score.

  - If provided `lambda.1se`, plot Rt which is generated from the lambda
    whose corresponding cross validation score is 1 standard error away
    of the minimal cross validation score.

  - If provided `cv_scores`, plot the cross validation score.

  - If NULL, all estimated Rt values are plotted.

- ...:

  Not used.

## Value

a [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)

## Examples

``` r
y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
cv <- cv_estimate_rt(y, korder = 1, nfold = 3, nsol = 30)
plot(cv)

plot(cv, which_lambda = cv$lambda[1])

plot(cv, which_lambda = "lambda.min")

plot(cv, which_lambda = "lambda.1se")

plot(cv, NULL)
```
