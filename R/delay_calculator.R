#' Calculate the total infectiousness at each observed time point.
#'
#' @details The total infectiousness at each observed time point is calculated by \eqn{\sum_{a=1}^t I_{t-s}w_s}
#' , where \eqn{I} denotes the vector containing all observed case counts, and \eqn{w} denotes the serial interval distribution.
#' The serial interval distribution expresses the probability of the symptom onset of a secondary infection occurred a given
#' number of days after the primary infection
#'
#' @usage delay_calculator(current_counts, dist_gamma)
#'
#' @param current_counts A vector of size n, containing all observed case counts
#' @param dist_gamma A vector of size 2, representing the shape and scale parameter of the discretized Gamma distribution
#'
#' @return A vector of size n-1, containing the total infectiousness at each observed time point
#' @export
#'
#' @examples
#' delay_calculator(c(3,2,5,3,1), c(2.5, 2.5))
delay_calculator <- function(current_counts, dist_gamma){
  n <- length(current_counts)
  w <- discretize_gamma(n, dist_gamma[1], dist_gamma[2])
  cw <- cumsum(w)
  convolved_seq <- convolve(current_counts, rev(w))[1:n]/cw
  c(convolved_seq[1], convolved_seq[1:(n-1)])
}
