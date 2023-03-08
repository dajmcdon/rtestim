#' @method summary poisson_rt
#' @export
summary.poisson_rt <- function(object, ...) {
  rlang::check_dots_empty()
  ns <- length(object$lambda)
  if (ns > 5) {
    xlam <- round(stats::quantile(1:ns))
    names(xlam) <- rev(c("Max.", "3rd Qu.", "Median", "1st Qu.", "Min."))
  } else {
    xlam <- seq_len(ns)
    names(xlam) <- paste0("s", seq_len(ns))
  }
  tab <- with(object, data.frame(
    lambda = lambda[xlam],
    index = xlam,
    approx_dof = dof[xlam],
    niterations = niter[xlam]))
  lambda <- object$lambda
  rownames(tab) <- names(xlam)
  out <- structure(
    list(call = object$call, table = tab, degree = object$degree, nlam = ns),
    class = "summary.poisson_rt")
  out
}

#' @method print summary.poisson_rt
#' @export
print.summary.poisson_rt <- function(x,
                                     digits = max(3, getOption("digits") - 3),
                                     ...) {
  rlang::check_dots_empty()
  cat("\nCall: ", deparse(x$call), fill = TRUE)
  cat("\nDegree of the estimated piecewise polynomial curve:", x$degree, "\n")
  cat("\nSummary of the", x$nlam, "estimated models:\n")
  print(x$tab, digits = digits)
  cat("\n")
}

#' @method print poisson_rt
#' @export
print.poisson_rt <- function(x, digits = min(3, getOption("digits") - 3), ...) {
  rlang::check_dots_empty()
  print(summary(x), digits = digits)
}

#' Plot estimated Rt values from a `poisson_rt` object
#'
#' Produces a figure showing some or all estimated Rt values for different
#' values of the penalty. The result is a [ggplot2::ggplot()]. Additional user
#' modifications can be added as desired.
#'
#'
#' @param x output of the function [estimate_rt()] of class `poisson_rt`
#' @param which_lambda select which Rt's to plot. If not provided,
#'   all Rt's are plotted.
#' @param ... Not used.
#'
#' @export
#'
#' @importFrom rlang .data
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
#' out <- estimate_rt(y, lambda = log(c(1.1,1.3,1.5)))
#' plot(out)
plot.poisson_rt <- function(x, which_lambda = NULL, ...) {
  arg_is_numeric(which_lambda, allow_null = TRUE)

  n <- length(x$observed_counts)
  if (is.null(which_lambda)) {
    Rt <- x$Rt
    lambda <- x$lambda
  } else {
    Rt <- fitted(x, lambda = which_lambda)
    lambda <- which_lambda
  }
  k <- length(lambda)

  df <- data.frame(
    Rt = c(Rt),
    lambda = rep(lambda, each = n),
    Time = rep(x$x, k)
  )

  plt <- ggplot2::ggplot(
    df,
    ggplot2::aes(.data$Time, .data$Rt, colour = .data$lambda,
                 group = .data$lambda)) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_viridis_c(trans = "log10")
  if (k == 1) plt <- plt + ggplot2::theme(legend.position= "none")
  plt
}

#' @importFrom stats fitted
#' @export
fitted.poisson_rt <- function(object, lambda = NULL, ...) {
  rlang::check_dots_empty()

  if (is.null(lambda)) return(object$Rt)

  lam_list <- interpolate_lambda(object$lambda, lambda)
  k <- length(lambda)
  log_r <- log(object$Rt)
  ret <- log_r[ ,lam_list$left, drop = FALSE] %*% diag(lam_list$frac, k, k) +
    log_r[ ,lam_list$right, drop = FALSE] %*% diag(1 - lam_list$frac, k, k)
  drop(exp(ret))
}

#' @importFrom stats coef
#' @export
coef.poisson_rt <- fitted.poisson_rt

#' Predict observed data using estimated Rt
#'
#' Given an object of class `poisson_rt` produced with [estimate_rt()],
#' calculate predicted observed cases for the estimated Rt values.
#' Note: This function is not indented for "new" or to produce forecasts, but
#' rather to examine how Rt relates to observables.
#'
#' @param object An object of class `poisson_rt` produced with [estimate_rt()].
#' @param lambda Select which lambdas from the object to use. If not provided
#'   (the default), all are returned. Note that new lambdas not originally
#'   used in the estimation procedure may be provided, but the results will be
#'   calculated by linearly interpolating the estimated Rt's.
#' @param ... Not used.
#'
#' @return A vector or matrix of predicted case counts.
#' @export
#'
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
#' out <- estimate_rt(y, nsol = 10)
#' preds <- predict(out)
#' plot(y)
#' matlines(preds, lty = 1)
predict.poisson_rt <- function(object, lambda = NULL, ...) {
  rlang::check_dots_empty()

  Rt <- fitted(object, lambda)
  Rt * object$weighted_past_counts
}
