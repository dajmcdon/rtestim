interpolate_lambda <- function(lambda, new_lambda) {
  r <- range(lambda)
  rr <- range(new_lambda)
  if (rr[1] < r[1]) {
    cli::cli_warn(
      c("Supplied `lambda` has values smaller than those used to estimate Rt originally.",
        i = "You may want to refit at these values.",
        i = "Using the smallest available `lambda`'s available."))
  }
  if (rr[2] > r[2]) {
    cli::cli_warn(
      c("Supplied `lambda` has values larger than those used to estimate Rt originally.",
        i = "You may want to refit at these values.",
        i = "Using the largest available `lambda`'s available."))
  }
  if (length(lambda) == 1) {# degenerate case of only one lambda
    nums <- length(new_lambda)
    left <- rep(1, nums)
    right <- left
    sfrac <- rep(1, nums)
  } else {
    k <- length(lambda)
    sfrac <- (lambda[1] - new_lambda) / (lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda) / (lambda[1] - lambda[k])
    sfrac[sfrac < min(lambda)] <- min(lambda)
    sfrac[sfrac > max(lambda)] <- max(lambda)
    coord <- approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac <- (sfrac - lambda[right]) / (lambda[left] - lambda[right])
    sfrac[left==right] <- 1
    sfrac[abs(lambda[left] - lambda[right]) < .Machine$double.eps] <- 1
  }
  list(left=left,right=right,frac=sfrac)
}
