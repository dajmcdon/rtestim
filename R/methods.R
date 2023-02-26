

#' Print `poisson_rt` object
#'
#' @param x output of function `estimate_rt` of class `poisson_rt`
#' @param ... further arguments passed to or from other methods.
#'
#' @usage \method{print}{poisson_rt}(x, \dots)
#' @return status of the `poisson_rt` object
#' @exportS3Method print poisson_rt
#'
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
#' out <- estimate_rt(y, lambda = log(c(1.1,1.3,1.5)))
#' out
print.poisson_rt <- function(x, ...) {
  cat("Algorithm terminated\n")
  if (all(x$convergence)) cat("All runs converged!\n")
  cat("Degree of the piecewise polynomial curve fitted:", x$degree, "\n")
}

#' Summary of the `poisson_rt` object
#'
#' @param object output of function `estimate_rt` of class `poisson_rt`
#' @param ... further arguments passed to or from other methods.
#'
#' @usage \method{summary}{poisson_rt}(object, \dots)
#' @return summary of the `poisson_rt` object in a table
#' @exportS3Method summary poisson_rt
#'
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
#' out <- estimate_rt(y, lambda = log(c(1.1,1.3,1.5)))
#' summary(out)
summary.poisson_rt <- function(object, ...) {
  lambda <- object$lambda
  nsol <- length(lambda)
  tab <- matrix(0, nsol, 2)
  tab[, 1] <- lambda
  tab[, 2] <- object$convergence
  colnames(tab) <- c("lambda", "convergence")

  tab <- as.table(tab)
  class(tab) <- "summary.poisson_rt"

  cat("Degree of the piecewise polynomial curve fitted:", object$degree, "\n")
  cat("\n")
  print(tab)
}

#' Plot `poisson_rt` object
#'
#' @param x output of function `estimate_rt` of class `poisson_rt`
#' @param which_lambda select which Rt's to plot. If no lambdas are provided,
#' all Rt's are plotted. Lambdas provided must match the lambda used in
#' generating the Rt's.
#' @param ... further arguments passed to or from other methods.
#'
#' @return plot of all or selected Rt from an object of class `poisson_rt`
#' @exportS3Method
#'
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
#' out <- estimate_rt(y, lambda = log(c(1.1,1.3,1.5)))
#' plot(out)
plot.poisson_rt <- function(x, which_lambda = NULL, ...) {
  arg_is_numeric(which_lambda, allow_null = TRUE)
  Rt <- x$Rt
  lambda <- x$lambda

  if (!is.null(which_lambda)) {
    if (!all(which_lambda %in% lambda))
      cli::cli_abort("Can only plot for lambda that used to generate Rt in
                     `estimate_rt`")
    idx <- which(lambda %in% which_lambda)
    Rt <- Rt[, idx]
    lambda <- lambda[idx]
  }

  plt_color <- c(1:length(lambda))

  graphics::matplot(Rt, type = "l", lty = 1, col = plt_color, main = "Estimated Rt",
          xlab = "Time")
  graphics::legend("topright",
         legend = round(lambda, 3), title = "Lambda", lty = 1,
         col = plt_color)
  plt <- grDevices::recordPlot()
  plt
}
