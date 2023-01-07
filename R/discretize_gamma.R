
#' discretized Gamma distribution given shape and scale parameters. This is intended to be the serial interval distribution of the infectious disease.
#'
#' @param n, length of past observations
#' @param shape shape parameter of the Gamma distribution
#' @param scale scale parameter of the Gamma distribution
#'
#' @return probability mass of discretized gamma distribution in a vector of size n-1
#' @export
#'
#' @examples discretize_gamma(30, shape = 2, scale = 2)
discretize_gamma <- function(n, shape, scale){
  stopifnot(shape > 0, scale > 0, n >= 1)
  x = 1:n
  pgm <- pgamma(x, shape = shape, scale = scale)
  pgm <- c(0, pgm)
  pgm <- diff(pgm)
  pgm <- pgm/sum(pgm)
  return(pgm)
}
