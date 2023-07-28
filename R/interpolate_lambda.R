interpolate_lambda <- function(lambda, new_lambda) {
  if (min(new_lambda) < min(lambda)) {
    cli::cli_warn(
      c("Requested `lambda` has values smaller than those used to estimate Rt.",
        i = "You may want to refit at these values.",
        i = "Using the smallest `lambda`'s available."))
  }
  if (max(new_lambda) > max(lambda)) {
    cli::cli_warn(
      c("Requested `lambda` has values larger than those used to estimate Rt.",
        i = "You may want to refit at these values.",
        i = "Using the largest `lambda`'s available."))
  }
  if (length(lambda) == 1) {
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
    coord <- stats::approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac <- (sfrac - lambda[right]) / (lambda[left] - lambda[right])
    sfrac[left == right] <- 1
    sfrac[abs(lambda[left] - lambda[right]) < .Machine$double.eps] <- 1
  }
  list(left = left, right = right, frac = sfrac)
}
