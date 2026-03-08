# Plot estimated Rt values from a `poisson_rt` object

Produces a figure showing some or all estimated Rt values for different
values of the penalty. The result is a
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html).
Additional user modifications can be added as desired.

## Usage

``` r
# S3 method for class 'poisson_rt'
plot(x, lambda = NULL, ...)
```

## Arguments

- x:

  output of the function
  [`estimate_rt()`](https://dajmcdon.github.io/rtestim/reference/estimate_rt.md)
  of class `poisson_rt`

- lambda:

  select which Rt's to plot. If not provided, all Rt's are plotted.

- ...:

  Not used.

## Value

a [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)

## Examples

``` r
y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
out <- estimate_rt(y, lambda = log(c(1.1, 1.3, 1.5)))
plot(out)
```
