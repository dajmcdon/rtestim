#' Generate an identity matrix
#'
#' @param n matrix dimension
#'
#' @return a sparse square matrix of dimension n*n
generate_I <- function(n) {
  stopifnot(n > 0, rlang::is_integerish(n))
  Matrix::bandSparse(n,
    k = 0,
    diagonals = list(rep(1, n))
  )
}

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
  D <- generate_I(n)
  if (degree > 0) D <- Matrix::diff(D, differences = degree)
  D
}
