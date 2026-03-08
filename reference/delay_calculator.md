# Calculate the total infectiousness at each observed time point.

The total infectiousness at each observed time point is calculated by
\\\sum\_{s=1}^t I\_{t-s}w_s\\, where \\I\\ denotes the vector containing
observed incidence, and \\w\\ denotes the generation interval
distribution. Typically, the generation interval is challenging to
estimate from data, so the serial interval is used instead. The serial
interval distribution expresses the probability of a secondary infection
caused by a primary infection which occurred \\s\\ days earlier.

## Usage

``` r
delay_calculator(
  observed_counts,
  x = NULL,
  dist_gamma = c(2.5, 2.5),
  delay_distn = NULL,
  delay_distn_periodicity = NULL,
  xout = x
)
```

## Arguments

- observed_counts:

  vector of the observed daily infection counts

- x:

  a vector of positions at which the counts have been observed. In an
  ideal case, we would observe data at regular intervals (e.g. daily or
  weekly) but this may not always be the case. May be numeric or Date.

- dist_gamma:

  Vector of length 2. These are the shape and scale for the assumed
  serial interval distribution. Roughly, this distribution describes the
  probability of an infectious individual infecting someone else after
  some period of time after having become infectious. As in most
  literature, we assume that this interval follows a gamma distribution
  with some shape and scale.

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

- xout:

  a vector of positions at which the results should be returned. By
  default, this will be the same as `x`, but in the case that
  observations are unequally spaced, alternatives may be desired. Note
  that `xout` must satisfy `min(x) <= min(xout)` and
  `max(x) >= max(xout)`.

## Value

A vector containing the total infectiousness at each point `xout`.

## Examples

``` r
delay_calculator(c(3, 2, 5, 3, 1), dist_gamma = c(2.5, 2.5))
#> [1] 3.000000 2.811312 2.828111 3.022142 3.119169
```
