# Plot estimated confidence bands for an estimate of Rt

Produces a figure showing a single estimated Rt value along with
approximate confidence bands. The result is a
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html).
Additional user modifications can be added as desired.

## Usage

``` r
# S3 method for class 'rt_confidence_band'
plot(x, colour = "#3A448F", ...)
```

## Arguments

- x:

  An object of class `rt_confidence_band` as produced by
  [`confband()`](https://dajmcdon.github.io/rtestim/reference/confband.md).

- colour:

  The colour of the desired plot

- ...:

  Not used.

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html).

## Examples

``` r
y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
out <- estimate_rt(y, nsol = 10)
cb <- confband(out, out$lambda[2], level = c(0.95, 0.8, 0.5))
plot(cb)

cb_y <- confband(out, out$lambda[2], level = c(0.95, 0.8, 0.5), type = "Yt")
plot(cb_y)
```
