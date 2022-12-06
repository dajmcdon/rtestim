#' Generate a divided difference matrix of order (degree - 1) when degree > 0;
#' generate an identity matrix when degree = 0
#'
#' @param n column dimension
#' @param degree order of the divided difference matrix, i.e., degree of the
#' piecewise polynomial curve to be fitted
#'
#' @return a sparse matrix of dimension (n - degree) * n
generate_D <- function(n, degree = 1) {
  stopifnot(n > 0, n > degree + 1, degree >= 0)
  stopifnot(rlang::is_integerish(n), rlang::is_integerish(degree))
  D <- Matrix::Diagonal(n)
  if (degree > 0) D <- Matrix::diff(D, differences = degree)
  D
}
