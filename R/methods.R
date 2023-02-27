
#' @exportS3Method summary poisson_rt
summary.poisson_rt <- function(object, ...) {
  "this method is not yet implemented"
}

#' @export
plot.poisson_rt <- function(x, ...) {
  plot(x$x, x$Rt, ty = "l")
}



#' Plot cv_result
#'
#' @param x result of cv_estimate_rt of class `cv_result`
#' @param ...
#'
#' @return plot of cv scores
#' @exportS3Method
#'
#' @examples
#' cv <- cv_estimate_rt(c(1:20), degree = 2, fold = 2, nsol=100)
#' plot(cv)
plot.cv_result <- function(x, ...) {
  par(mar = c(5,4,4,2), mfrow=c(1,2))
  plot(x$cv_scores, type = "l")
  plot(x$optimal_Rt)
}
