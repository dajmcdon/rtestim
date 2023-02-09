#' @exportS3Method summary poisson_rt
summary.poisson_rt <- function(object, ...) {
  "this method is not yet implemented"
}

#' @export
plot.poisson_rt <- function(x, ...) {
  plot(x$x, x$Rt, ty = "l")
}
