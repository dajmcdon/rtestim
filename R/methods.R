#' Summarize results of a admm_rr model
#'
#' @method summary admm_rr
#'
#' @param object a fitted model of class `admm_rr`
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
#' mod = admm_solver(
#'   current_counts = y, weighted_past_counts = rep(1, 10), degree = 1,
#'   init = admm_initializer(current_counts = y,
#'   weighted_past_counts = rep(1, 10), degree = 1)
#' )
#' summary(mod)
summary.admm_rr <- function(object, ...){
  n = length(object$current_counts)
  res <- data.frame(
    Time = 1:n,
    Signal = object$current_counts,
    R_rate = object$R_rate,
    pois_mean = object$R_rate * object$weighted_past_counts
  )

  lst = list(Results = res, Convergence = object$convr)
  class(lst) = "summary.admm_rr"
  return(lst)
}

#' Plot summary of `admm_rr` models
#'
#' @method plot summary.admm_rr
#' @param x summary of `admm_rr` models
#' @param ... .
#'
#' @return a figure
#' @export
#'
#' @examples
#' y <- c(rev(seq(2, 6, by = 1)), seq(2, 6, by = 1))
#' mod = admm_solver(
#'   current_counts = y, weighted_past_counts = rep(1, 10), degree = 1,
#'   init = admm_initializer(current_counts = y,
#'   weighted_past_counts = rep(1, 10), degree = 1)
#' )
#' plot(summary(mod))
plot.summary.admm_rr <- function(x, ...){
  fig <- x$Results %>%
    ggplot(aes(x = .data$Time)) +
    geom_point(aes(y = .data$Signal)) +
    geom_line(aes(y = .data$pois_mean), col = "#08519C") +
    labs(x = "Time", y = "Daily infection counts (on dots)",
         title = "The estimated piecewise polynomial curve (in line)") +
    theme_bw()
  print(fig)
}
