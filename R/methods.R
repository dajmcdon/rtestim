#' @exportS3Method summary poisson_rt
summary.poisson_rt <- function(object, ...) {
  "this method is not yet implemented"
}

#' Plot summary of `poisson_rt` models
#'
#' @method plot summary.poisson_rt
#' @param x summary of `poisson_rt` models
#' @param ... .
#'
#' @return a figure
#' @export
plot.poisson_rt <- function(x, ...) {
  plot(x$x, x$Rt, ty = "l")
}
