#' Interpolate (or extrapolate) Rt estimates to intermediate design points
#'
#' @param object A fitted object produced by `estimate_rt()` or `cv_estimate_rt()`.
#' @param xout a vector of new positions at which Rt should be produced,
#'   but where counts may not have been observed.
#' @param ... additional arguments passed to methods.
#'
#' @return A vector or matrix of interpolated Rt estimates.
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
#' out <- estimate_rt(y)
#'
#' # originally estimated at
#' out$x
#'
#' # get the Rt at 3 new points (for all estimated lambdas)
#' int <- interpolate_rt(out, c(10.5, 11.5, 12.5))
#'
#' # get the Rt at a single value of lambda
#' interpolate_rt(out, c(10.5, 11.5, 12.5), lambda = out$lambda[20])
#'
#' @export
interpolate_rt <- function(object, xout, ...) {
  if (inherits(xout, "Date")) xout <- as.numeric(xout)
  arg_is_numeric(xout)
  UseMethod("interpolate_rt")
}

#' @export
interpolate_rt.default <- function(object, xout, ...) {
  cli_abort(
    "No interpolation methods exist for objects of class {.cls {class(object)}}."
  )
}
