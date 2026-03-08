# Changelog

## rtestim (development version)

## rtestim 1.0.2

- Rerun
  [`Rcpp::compileAttributes()`](https://rdrr.io/pkg/Rcpp/man/compileAttributes.html)
  to guard against suspicious `Rf_error()` (email from
  [Rcpp](https://www.rcpp.org) team).

## rtestim 1.0.1

CRAN release: 2025-10-24

- Added a `NEWS.md` file to track changes to the package.
- Account for double interval censoring in
  [`discretize_gamma()`](https://dajmcdon.github.io/rtestim/reference/discretize_gamma.md)
  [\#81](https://github.com/dajmcdon/rtestim/issues/81).
- Bug fix for `korder=0` introduced by the change to
  [tvdenoising](https://github.com/glmgen/tvdenoising).

## rtestim 1.0.0

CRAN release: 2025-07-05

- Initial CRAN submission
