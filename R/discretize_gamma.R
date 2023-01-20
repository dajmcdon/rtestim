#' Density function for the discretized gamma function with the given shape and scale parameter
#'
#' @details The serial interval distribution expresses the probability of the symptom onset of a secondary infection
#' occurred a given number of days after the primary infection. The serial interval distribution is commonly represented
#' by a discretized Gamma distribution in literature, parametrized by the shape and scale parameters.
#'
#' @usage discretize_gamma(n, shape=2.5, scale=2.5)
#' @param n, length of the observed case count
#' @param shape shape parameter of the discretized Gamma distribution
#' @param scale scale parameter of the discretized Gamma distribution
#'
#' @return probability mass of the discretized gamma distribution in a vector of size n
#' @export
#'
#' @examples
#' discretize_gamma(30, shape = 2.5, scale = 2.5)
discretize_gamma <- function(n, shape = 2.5, scale = 2.5){
  stopifnot(shape > 0, scale > 0, n >= 1, rlang::is_intergerish(n))
  x <- 1:n
  pgm <- pgamma(x, shape = shape, scale = scale)
  pgm <- c(0, pgm)
  pgm <- diff(pgm)
  pgm/sum(pgm)
}
