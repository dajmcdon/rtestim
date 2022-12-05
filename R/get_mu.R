#' Compute mu: the largest singular value of D multiplied by two times lambda
#'
#' @param n length of the current observed daily counts
#' @param degree degree of the piecewise polynomial curve to be fitted,
#' e.g., degree = 0 corresponds to piecewise constant curves
#' @param lambda a parameter to balance the data fidelity and graphical
#' smoothness of fitted curves, e.g., a greater lambda results in better
#' smoothness
#'
#' @return mu
#' @export
#'
#' @examples get_mu(5, 1, 0.1)
get_mu <- function(n, degree, lambda) {
  mu <- 2 * lambda * (4^degree - 1/n)
  mu
}
