#' Find the knots for an estimated `poisson_rt` object.
#'
#' Produces the lengths of the segments between knots as well as the number
#' of knots.
#'
#' @param object A fitted `poisson_rt` object.
#' @param lambda A scalar value of the penalty parameter. If this value is not
#'   present in the `object$lambda` sequence, the result will be approximate
#'   based on interpolating linearly between the available values.
#'
#' @return A list containing the segment lengths and the number of segments.
#' @keywords internal
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
#' out <- estimate_rt(y)
#' kn <- find_knots(out, out$lambda[20])
#' @noRd
find_knots <- function(object, lambda, ...) {
  rlang::check_dots_empty()
  assert_number(lambda, lower = min(object$lambda), upper = max(object$lambda))

  n <- length(object$observed_counts)
  lam_list <- interpolate_lambda(object$lambda, lambda)
  alp <- interpolate_mat(object$alp, lam_list, lambda, take_log = FALSE)
  dof <- object$dof
  dof <- dof[lam_list$left] * lam_list$frac +
    dof[lam_list$right] * (1 - lam_list$frac)
  knots <- which(abs(alp) > 1e-10)
  r <- c(knots + object$korder, n)
  l <- c(1, knots + object$korder + 1)
  xpieces <- map2(l, r, function(a, b) object$x[a:b])
  structure(enlist(knots, xpieces, dof, l, r), class = "rtestim_knots")
}


