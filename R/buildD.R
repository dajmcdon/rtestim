#' Generate an identity matrix
#'
#' @param n matrix dimension
#'
#' @return a sparse matrix
generate_I <- function(n) {
  stopifnot(n > 0, rlang::is_integerish(n))
  Matrix::bandSparse(n,
    k = 0,
    diagonals = list(rep(1, n))
  )
}

#' Generate a divided difference matrix of order zero
#'
#' @param n number of matrix rows
#'
#' @return a (n-1) * n sparse matrix
generate_D <- function(n) {
  stopifnot(n > 0, rlang::is_integerish(n))
  Matrix::bandSparse(n,
    k = c(0, 1),
    diagonals = list(rep(-1, n), rep(1, n - 1))
  )[-n, ]
}

#' Generate a divided difference matrices of order k
#'
#' @param n row dimension
#' @param k order of divided difference
#'
#' @return a (n-k-1) * n sparse matrix
generate_Dk <- function(n, k = 1) {
  stopifnot(n > 0, n > k + 1, k >= 0)
  stopifnot(rlang::is_integerish(n), rlang::is_integerish(k))
  D0 <- generate_D(n)
  if (k > 0) D0 <- Matrix::diff(D0, differences = k)
  D0
}

#' Compute mu
#'
#' @param n length of observed sequence
#' @param k order of divided difference
#' @param lambda tuning parameter
#'
#' @return mu
#' @export
#'
#' @examples get_mu(5, 1, 0.1)
get_mu <- function(n, k, lambda) {
  if (k == 0) {
    mu <- lambda * 2
  } else {
    Dk <- generate_Dk(n, k - 1)
    mu <- RSpectra::eigs(spam::crossprod.spam(Dk),
      k = 1,
      which = "LM"
    )$values * lambda * 2
  }
  mu
}
