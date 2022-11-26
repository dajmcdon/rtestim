#' Summarize results of a admm_rr model
#'
#' @method summary admm_rr
#'
#' @param mod a fitted admm_rr model
#'
#' @return a data table of estimates and a logical value of convergence
#'
#' @rdname summary.admm_rr
#'
#' @export
summary.admm_rr <- function(mod){
  n = length(mod$y)
  res <- data.table::data.table(
    Time = 1:n,
    Signal = mod$y,
    R_rate = as.vector(mod$R_rate),
    pois_para = as.vector(mod$R_rate * mod$x)
  )
  lst = list(Results = res, Convergence = mod$convr)
  class(lst) = "summary.admm_rr"
  return(lst)
}

#' Plot a `admm_rr' class object
#'
#' @method plot admm_rr
#'
#' @param mod a fitted admm_rr model
#'
#' @return a plot
#'
#' @rdname plot.admm_rr
#'
#' @export
plot.admm_rr <- function(mod) {
  res <- summary.admm_rr(mod)$Results
  fig <- res %>%
    tidyr::pivot_longer(R_rate:pois_para, names_to = "Type",
                 values_to = "Est") %>%
    dplyr::group_by(Type) %>%
    ggplot2::ggplot(ggplot2::aes(x = Time)) +
    ggplot2::geom_point(ggplot2::aes(y = Signal)) +
    ggplot2::geom_line(ggplot2::aes(y = Est, colour = Type, linetype = Type)) +
    ggplot2::labs(x = "Time", y = "Signal") +
    ggplot2::scale_colour_manual(values = RColorBrewer::brewer.pal(n = 4,
                                   name = "Paired")) +
    ggplot2::theme_bw()
  fig
}
