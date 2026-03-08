# Compute the discretized density function for gamma distribution

The serial interval distribution expresses the probability of the
symptom onset of a secondary infection occurred a given number of days
after the primary infection. The serial interval distribution is
commonly represented by a discretized Gamma distribution in literature,
parametrized by the shape and scale parameters.

## Usage

``` r
discretize_gamma(x, shape = 2.5, scale = 2.5, rate = 1/scale)
```

## Arguments

- x:

  locations (times) where cases are observed. Must be nonnegative.

- shape, scale:

  shape and scale parameters. Must be positive, `scale` strictly.

- rate:

  an alternative way to specify the scale.

## Value

probability mass of the discretized gamma distribution

## Examples

``` r
discretize_gamma(1:30, shape = 1, scale = 1)
#>  [1] 6.321206e-01 2.325442e-01 8.554821e-02 3.147143e-02 1.157769e-02
#>  [6] 4.259195e-03 1.566870e-03 5.764193e-04 2.120528e-04 7.800987e-05
#> [11] 2.869823e-05 1.055749e-05 3.883883e-06 1.428801e-06 5.256264e-07
#> [16] 1.933671e-07 7.113580e-08 2.616940e-08 9.627183e-09 3.541643e-09
#> [21] 1.302898e-09 4.793092e-10 1.763280e-10 6.486749e-11 2.386338e-11
#> [26] 8.778840e-12 3.229589e-12 1.188076e-12 4.370671e-13 1.607855e-13
```
