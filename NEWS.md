# rtestim (development version)

* Correct minor issue in the construction of the time-varying delay example
in `vignettes/articles/delay-distributions.Rmd`

# rtestim 1.0.2

* Rerun `Rcpp::compileAttributes()` to guard against suspicious `Rf_error()` (email from `{Rcpp}` team).

# rtestim 1.0.1

* Added a `NEWS.md` file to track changes to the package.
* Account for double interval censoring in `discretize_gamma()` #81.
* Bug fix for `korder=0` introduced by the change to `{tvdenoising}`.

# rtestim 1.0.0

* Initial CRAN submission
