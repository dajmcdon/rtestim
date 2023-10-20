interpolate_mat <- function(mat, lam_list, lambda, take_log = FALSE) {
  m <- length(lambda)
  if (take_log) mat <- log(mat)
  out <- mat[, lam_list$left, drop = FALSE] %*% diag(lam_list$frac, m, m) +
    mat[, lam_list$right, drop = FALSE] %*% diag(1 - lam_list$frac, m, m)
  if (take_log) out <- exp(out)
  drop(out)
}
