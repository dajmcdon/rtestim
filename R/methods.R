#' Summarize results of a admm_rr model
#'
#' @method summary admm_rr
#'
#' @param object a fitted admm_rr model
#' @param ... .
#'
#' @return a data table of estimates and a logical value of convergence
#'
#' @rdname summary.admm_rr
#'
#' @export
summary.admm_rr <- function(object, ...){
  n = length(object$y)
  res <- data.table::data.table(
    Signal = object$y,
    R_rate = as.vector(object$R_rate),
    pois_para = as.vector(object$R_rate * object$x)
  )

  lst = list(Results = res, Convergence = object$convr)
  class(lst) = "summary.admm_rr"
  return(lst)
}

#' Plot a `admm_rr' class object
#'
#' @param x time
#' @param y signal; estimated poisson mean
#' @param ... .
#'
#' @method plot admm_rr
#'
#' @return a plot
#'
#' @rdname plot.admm_rr
#'
#' @export
#'
#' @examples
#' y <- c(rev(seq(2, 6, by = 1)), seq(2, 6, by = 1))
#' mod = admm_solver(
#'   y = y, x = rep(1, 10), k = 1,
#'   init = admm_initializer(y = y, x = rep(1, 10), k = 1)
#' )
#' res <- summary(mod)$Results
#' plot(res$pois_para, y)
plot.admm_rr <- function(x, y, ...) {
  res <- data.table::data.table(Signal = y,
                                Time = 1:length(y),
                                pois_para = x)
  fig <- res %>%
    ggplot2::ggplot(ggplot2::aes(x = Time)) +
    ggplot2::geom_point(ggplot2::aes(y = Signal)) +
    ggplot2::geom_line(ggplot2::aes(y = pois_para)) +
    ggplot2::labs(x = "Time", y = "Signal") +
    ggplot2::scale_colour_manual(values = RColorBrewer::brewer.pal(n = 2,
                                   name = "Paired")) +
    ggplot2::theme_bw()
  fig
  class(fig) = "plot.admm_rr"
}
