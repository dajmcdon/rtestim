#' Summarize results of a poisson_rt model
#'
#' @method summary poisson_rt
#'
#' @param object a fitted model of class `poisson_rt`
#' @param ... .
#'
#' @return a data table of estimates and a logical value of convergence.
#' The data table columns include the current observed daily counts (Signal),
#' the estimated reproduction rate (R_rate), and the estimated Poisson mean
#' parameter (pois_mean)
#'
#' @export
#'
#' @examples
#' y <- c(rev(seq(2, 6, by = 1)), seq(2, 6, by = 1))
#' mod <- estimate_rt(
#'   observed_counts = y, degree = 1L, lambda = c(.1,.5),
#'   algo = "linear_admm",
#'   init = rt_admm_configuration(y, degree = 1L)
#' )
#' summary(mod)
summary.poisson_rt <- function(object, ...){
  n <- length(object$observed_counts)
  Results <- cbind(
    Observed_counts = object$observed_counts,
    Estimated_mean <- object$Rt * object$weighted_past_counts
  )
  num = dim(Results)[2]-1
  colnames(Results)[-1] <- paste("Estimated_mean_",1:num, sep="")

  lst = list(Results=as.data.frame(Results), Convergence=object$convergence)
  class(lst) = "summary.poisson_rt"
  return(lst)
}

#' Plot summary of `poisson_rt` models
#' (We may change it to draw multiple curves(estimates) simultaneously.)
#'
#' @method plot summary.poisson_rt
#' @param x summary of `poisson_rt` models
#' @param ... .
#'
#' @return a figure
#' @export
#'
#' @examples
#' y <- c(rev(seq(2, 6, by = 1)), seq(2, 6, by = 1))
#' mod <- estimate_rt(y, dist_gamma=c(1,1), degree = 2L, lambda = c(.1,.5),
#'   algo = "linear_admm", init = rt_admm_configuration(y, degree = 2L)
#' )
#' plot(summary(mod))
plot.summary.poisson_rt <- function(x, ...){
  n = dim(x$Results)[1]
  fig <- x$Results %>%
    ggplot(aes(x = 1:n)) +
    geom_point(aes(y = .data$Observed_counts)) +
    geom_line(aes(y = .data$Estimated_mean_1), col = "#08519C") +
    labs(x = "Time", y = "Daily infection counts (on dots)",
         title = "The estimated piecewise polynomial curve (in line)") +
    theme_bw()
  print(fig)
}
