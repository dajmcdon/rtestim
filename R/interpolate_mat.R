interpolate_mat <- function(mat, lam_list, lambda) {
  m <- length(lambda)
  log_mat <- log(mat)
  out <- log_mat[,lam_list$left, drop = FALSE] %*% diag(lam_list$frac, m, m) +
    log_mat[,lam_list$right, drop = FALSE] %*% diag(1 - lam_list$frac, m, m)
  drop(exp(out))
}
