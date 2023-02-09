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
  stopifnot(shape > 0, scale > 0, all(x >= 0))
  pgm <- stats::pgamma(x, shape = shape, scale = scale)
  pgm <- c(0, pgm)
  pgm <- diff(pgm)
  pgm / sum(pgm)
}
