
#' Summary of the `poisson_rt` object
#'
#' @param object
#' @param ...
#'
#' @return
#' @exportS3Method summary poisson_rt
#'
#' @examples
summary.poisson_rt <- function(object, ...) {
  res <- list(
    Time = object$x, # This need to change, it shows (x - x[1]) / diff(range(x)) * n
    Observed_time_points = object$observed_counts,
    R_rate = object$Rt,
    Predicted_cases = object$Rt * object$weighted_past_counts
  )

  class(res) = "summary.poisson_rt"
  return(res)
}

#' Plot `poisson_rt` object
#'
#' @param x output of function `estimate_rt` of class `poisson_rt`
#' @param ...
#'
#' @return
#' @exportS3Method
#'
#' @examples
#' y <- c(rev(seq(2, 6, by = 1)), seq(2, 6, by = 1))
#' mod <- estimate_rt(observed_counts = y, degree = 1, lambda = 0.001)
#' plot(mod)
plot.poisson_rt <- function(x, ...) {
  observed_counts <- x$observed_counts
  Rt <- x$Rt
  weighted_past_counts <- x$weighted_past_counts
  obs_freq <- x$x
  signal <- Rt*weighted_past_counts

  fig_cases <- ggplot2::ggplot(data.frame(obs_freq = obs_freq,
                                          observed_counts = observed_counts,
                                          signal = signal))+
    ggplot2::geom_point(aes(x = obs_freq, y = observed_counts))+
    ggplot2::geom_line(aes(x = obs_freq, y = signal))+
    ggplot2::labs(x = "Time", y = "Daily infection counts (on dots)",
                  title = "Estimated piecewise polynomial curve (in line)") +
    ggplot2::theme_bw()

  fig_rt <- ggplot(data.frame(obs_freq = obs_freq,
                              Rt = Rt))+
    ggplot2::geom_line(aes(x = obs_freq, y = Rt), col = "#14754C") +
    ggplot2::labs(x = "Time", y = "Estimated Rt",
                  title = "Estimated Rt") +
    ggplot2::theme_bw()

  fig <- ggpubr::ggarrange(fig_cases, fig_rt, ncol=2)

  print(fig)
}
