
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rtestim

<!-- badges: start -->

[![R-CMD-check](https://github.com/dajmcdon/rtestim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dajmcdon/rtestim/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This package uses l1-regularized nonparametric regression (trend
filtering) to estimate the effective reproductive number, Rt. This value
roughly says “how many new infections will result from each new
infection today”. Values larger than 1 indicate that an epidemic is
growing while those less than 1 indicate decline.

## Installation

You can install the development version of rtestim from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dajmcdon/rtestim")
```

## Quick example

Here we create so data that “looks” like a typical wave in an epidemic.
Because the model uses regularized regression, we estimate the model at
a range of tuning parameters simultaneously.

``` r
set.seed(12345)
library(rtestim)
library(ggplot2)
dat <- data.frame(
  Time = 1:101, 
  incident_cases = c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
)
ggplot(dat, aes(Time, incident_cases)) +
  geom_point(colour = "cornflowerblue") +
  theme_bw()
```

<img src="man/figures/README-plot-data-1.png" width="100%" />

We fit the model and visualize the resulting estimated sequence of
$R_t$:

``` r
mod <- estimate_rt(observed_counts = dat$incident_cases, nsol = 20)
plot(mod)
```

<img src="man/figures/README-full-fit-1.png" width="100%" />

The additional parameter `nsol = 20` specify the number of $\lambda$s
for which $R_t$ is estimated.

A built in function for cross-validation can be used to select the
tuning parameter.

``` r
mod_cv <- cv_estimate_rt(dat$incident_cases, nsol = 20)
plot(mod_cv, which_lambda = 'lambda.1se')
```

<img src="man/figures/README-cv-estimate-1.png" width="100%" />
