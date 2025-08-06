#' Compute the discretized density function for gamma distribution
#'
#' The serial interval distribution expresses the probability of the
#' symptom onset of a secondary infection occurred a given number of days after
#' the primary infection. The serial interval distribution is commonly
#' represented by a discretized Gamma distribution in literature, parametrized
#' by the shape and scale parameters.
#'
#' @param x locations (times) where cases are observed. Must be nonnegative.
#' @inheritParams stats::pgamma
#'
#'
#'
#' @return probability mass of the discretized gamma distribution
#' @export
#'
#' @examples
#' discretize_gamma(1:30, shape = 1, scale = 1)
discretize_gamma <- function(x, shape = 2.5, scale = 2.5, rate = 1 / scale) {
  assert_number(shape, lower = 0, finite = TRUE)
  assert_number(scale, lower = 0, finite = TRUE)
  assert_number(rate, lower = 0, finite = TRUE)
  if (!missing(scale) && !missing(rate)) {
    if (abs(rate * scale - 1) < 1e-15) {
      cli_warn("specify `rate` or `scale` but not both")
    } else {
      cli_abort("specify `rate' or `scale` but not both")
    }
  }
  assert_numeric(x, lower = 0, any.missing = FALSE)
  if (is.unsorted(x, strictly = TRUE)) {
    cli_abort("`x` must be sorted in increasing order.")
  }
  pgmr <- stats::pgamma(x + 1, shape = shape, rate = rate)
  pgml <- stats::pgamma(x - 1, shape = shape, rate = rate)
  pgm <- pgmr - pgml
  pgm / sum(pgm)
}
