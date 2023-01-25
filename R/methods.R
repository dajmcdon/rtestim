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
#' TODO: Need to change example
#' y <- c(rev(seq(2, 6, by = 1)), seq(2, 6, by = 1))
#' mod = admm_solver(
#'   current_counts = y, weighted_past_counts = rep(1, 10), degree = 1,
#'   init = admm_initializer(current_counts = y,
#'   weighted_past_counts = rep(1, 10), degree = 1)
#' )
#' summary(mod)
summary.poisson_rt <- function(object, ...){
  n = length(object$observed_counts)
  res <- data.frame(
    Time = 1:n,
    Signal = object$observed_counts,
    R_rate = object$Rt,
    pois_mean = object$Rt * object$weighted_past_counts
  )

  lst = list(Results = res)
  class(lst) = "summary.poisson_rt"
  return(lst)
}



#' Plot summary of `poisson_rt` models
#'
#' @method plot summary.poisson_rt
#' @param x summary of `poisson_rt` models
#' @param ... .
#'
#' @return a figure
#' @export
#'
#' @examples
#' TODO: change this example
#' y <- c(rev(seq(2, 6, by = 1)), seq(2, 6, by = 1))
#' mod = admm_solver(
#'   current_counts = y, weighted_past_counts = rep(1, 10), degree = 1,
#'   init = admm_initializer(current_counts = y,
#'   weighted_past_counts = rep(1, 10), degree = 1)
#' )
#' plot(summary(mod))
plot.summary.poisson_rt <- function(summary, ...){
  fig_cases <- summary$Results %>%
    ggplot(aes(x = .data$Time)) +
    geom_point(aes(y = .data$Signal)) +
    geom_line(aes(y = .data$pois_mean), col = "#08519C") +
    labs(x = "Time", y = "Daily infection counts (on dots)",
         title = "The estimated piecewise polynomial curve (in line)") +
    theme_bw()

  fig_rt <- summary$Results %>%
    ggplot(aes(x = .data$Time)) +
    geom_line(aes(y = .data$R_rate), col = "#14754C") +
    labs(x = "Time", y = "Estimated Rt",
         title = "The estimated Rt") +
    theme_bw()

  fig <- ggpubr::ggarrange(fig_cases, fig_rt, ncol=2)

  print(fig)
}

