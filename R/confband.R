
confband <- function(object, lambda = NULL, level = 0.95, ...) {
  nbd <- function(n, ord) {
    if (n < ord + 1) return(Diagonal(n, x = 0))
    buildD(n, ord)
  }

  y <- object$observed_counts
  n <- length(y)
  rhat <- fitted(object, lambda)
  yhat <- predict(object, lambda)
  k <- object$degree
  knots <- find_knots(object, lambda)
  Ds <- Matrix::bdiag(lapply(knots$lens, nbd, ord = k))
  # still not sure about this line
  covs <- diag(solve(Matrix::Diagonal(n, yhat) + lambda * crossprod(Ds))) * rhat^2
  return(list(
    Rt = rhat,
    lower = pmax(rhat + qt(1 - level / 2, n - knots$dof) * sqrt(covs), 0),
    upper = pmax(rhat + qt(level / 2, n - knots$dof) * sqrt(covs), 0)
  ))
}



interpolate_mat <- function(mat, lam_list, lambda) {
  m <- length(lambda)
  log_mat <- log(mat)
  out <- log_mat[,lam_list$left, drop = FALSE] %*% diag(lam_list$frac, m, m) +
    log_mat[,lam_list$right, drop = FALSE] %*% diag(1 - lam_list$frac, m, m)
  drop(exp(out))
}

