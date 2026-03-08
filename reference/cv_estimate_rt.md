# Leave-kth-out cross validation for choosing a optimal parameter lambda

Leave-kth-out cross validation for choosing a optimal parameter lambda

## Usage

``` r
cv_estimate_rt(
  observed_counts,
  korder = 3L,
  dist_gamma = c(2.5, 2.5),
  nfold = 3L,
  error_measure = c("deviance", "mse", "mae"),
  x = 1:n,
  lambda = NULL,
  maxiter = 1000000L,
  delay_distn = NULL,
  delay_distn_periodicity = NULL,
  regular_splits = FALSE,
  invert_splits = FALSE,
  ...
)
```

## Arguments

- observed_counts:

  vector of the observed daily infection counts

- korder:

  Integer. Degree of the piecewise polynomial curve to be estimated. For
  example, `korder = 0` corresponds to a piecewise constant curve.

- dist_gamma:

  Vector of length 2. These are the shape and scale for the assumed
  serial interval distribution. Roughly, this distribution describes the
  probability of an infectious individual infecting someone else after
  some period of time after having become infectious. As in most
  literature, we assume that this interval follows a gamma distribution
  with some shape and scale.

- nfold:

  Integer. This number of folds to conduct the leave-kth-out cross
  validation. For leave-kth-out cross validation, every kth
  observed_counts and their corresponding position (evenly or unevenly
  spaced) are placed into the same fold. The first and last
  observed_counts are not assigned to any folds. Smallest allowable
  value is `nfold = 2`.

- error_measure:

  Metric used to calculate cross validation scores. Must be choose from
  `mse`, `mae`, and `deviance`. `mse` calculates the mean square error;
  `mae` calculates the mean absolute error; `deviance` calculates the
  deviance

- x:

  a vector of positions at which the counts have been observed. In an
  ideal case, we would observe data at regular intervals (e.g. daily or
  weekly) but this may not always be the case. May be numeric or Date.

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

- maxiter:

  Integer. Maximum number of iterations for the estimation algorithm.

- delay_distn:

  in the case of a non-gamma delay distribution, a vector or matrix (or
  [`Matrix::Matrix()`](https://rdrr.io/pkg/Matrix/man/Matrix.html)) of
  delay probabilities may be passed here. For a vector, these will be
  coerced to sum to 1, and padded with 0 in the right tail if necessary.
  If a time-varying delay matrix, it must be lower-triangular. Each row
  will be silently coerced to sum to 1. See also
  `vignette("delay-distributions")`.

- delay_distn_periodicity:

  Controls the relationship between the spacing of the computed delay
  distribution and the spacing of `x`. In the default case, `x` would be
  regular on the sequence `1:length(observed_cases)`, and this would
  be 1. But if `x` is a `Date` object or spaced irregularly, the
  relationship becomes more complicated. For example, weekly data when
  `x` is a date in the form `YYYY-MM-DD` requires specifying
  `delay_distn_periodicity = "1 week"`. Or if `observed_cases` were
  reported on Monday, Wednesday, and Friday, then
  `delay_distn_periodicity = "1 day"` would be most appropriate.

- regular_splits:

  Logical. If `TRUE`, the folds for k-fold cross-validation are chosen
  by placing every kth point into the same fold. The first and last
  points are not included in any fold and are always included in
  building the predictive model. As an example, with 15 data points and
  `kfold = 4`, the points are assigned to folds in the following way:
  \$\$ 0 \\ 1 \\ 2 \\ 3 \\ 4 \\ 1 \\ 2 \\ 3 \\ 4 \\ 1 \\ 2 \\ 3 \\ 4 \\
  1 \\ 0 \$\$ where 0 indicates no assignment. Therefore, the folds are
  not random and running `cv_estimate_rt()` twice will give the same
  result.

- invert_splits:

  Logical. Typical K-fold CV would use K-1 folds for the training set
  while reserving 1 fold for evaluation (repeating the split K times).
  Setting this to true inverts this process, using a much smaller
  training set with a larger evaluation set. This tends to result in
  larger values of `lambda` that minimize CV.

- ...:

  additional parameters passed to
  [`estimate_rt()`](https://dajmcdon.github.io/rtestim/reference/estimate_rt.md)
  function

## Value

An object with S3 class `"cv_poisson_rt"`. Among the list components:

- `full_fit` An object with S3 class `"poisson_rt"`, fitted with all
  `observed_counts` and `lambda`

- `cv_scores` leave-kth-out cross validation scores

- `cv_se` leave-kth-out cross validation standard error

- `lambda.min` lambda which achieved the optimal cross validation score

- `lambda.1se` lambda that gives the optimal cross validation score
  within one standard error.

- `lambda` the value of `lambda` used in the algorithm.

## Examples

``` r
y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
cv <- cv_estimate_rt(y, korder = 3, nfold = 3, nsol = 30)
cv
#> 
#> Call: cv_estimate_rt(observed_counts = y, korder = 3, nfold = 3, nsol = 30)
#> 
#> Degree of the estimated piecewise polynomial curve: 3 
#> 
#> Summary of cross validation across lambda:
#>                 lambda index cv_scores  cv_se dof
#> Max Lambda   357.20215     1     1.772 0.2181   4
#> 1se Lambda   357.20215     1     1.772 0.2181   4
#> CV Minimizer   1.17532    19     1.681 0.1887  10
#> Min Lambda     0.03572    30     3.471 1.0993  25
#> 
```
