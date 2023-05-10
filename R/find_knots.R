
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
#' @export
#' @keywords internal
find_knots <- function(object, lambda, ...) {
  UseMethod("find_knots")
}

#' @export
find_knots.poisson_rt <- function(object, lambda, ...) {
  rlang::check_dots_empty()
  arg_is_scalar(lambda)
  arg_is_numeric(lambda)

  n <- length(object$observed_counts)
  lam_list <- interpolate_lambda(object$lambda, lambda)
  alp <- interpolate_mat(object$alp, lam_list, lambda, take_log = FALSE)
  dof <- object$dof
  dof <- dof[lam_list$left] * lam_list$frac +
    dof[lam_list$right] * (1 - lam_list$frac)
  knots <- which(abs(alp) > object$tolerance)
  r <- c(knots + object$degree, n)
  l <- c(1, knots + object$degree + 1)
  pieces <- map2(l, r, function(a, b) object$x[a:b])
  list(knots = knots, pieces = pieces, dof = dof)
}
