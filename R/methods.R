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
#'   all Rt's are plotted. Any lambdas provided must match the values used in
#'   generating the Rt's.
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
  Rt <- x$Rt
  lambda <- x$lambda

  if (!is.null(which_lambda)) {
    if (!all(which_lambda %in% lambda))
      cli::cli_abort("Can only plot for lambda that used to generate Rt in
                     `estimate_rt()`")
    idx <- which(lambda %in% which_lambda)
    Rt <- as.matrix(Rt[, idx])
    lambda <- lambda[idx]
  }

  df <- data.frame(
    Rt = c(Rt),
    lambda = rep(lambda, each = nrow(Rt)),
    Time = rep(x$x, ncol(Rt))
  )

  ggplot2::ggplot(
    df,
    ggplot2::aes(.data$Time, .data$Rt, colour = .data$lambda,
                 group = .data$lambda)) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_viridis_c(trans = "log10")
}

#' @method summary cv_result
#' @export
summary.cv_result <- function(object, ...) {
  rlang::check_dots_empty()
  cv_scores <- object$cv_scores
  ns <- length(cv_scores)

  xcv <- order(cv_scores)
  if (ns > 5) xcv <- xcv[1:5]

  tab <- with(object, data.frame(
    cv_scores = cv_scores[xcv],
    index = xcv,
    lambda = object$lambda[xcv])
    )
  out <- structure(
    list(call = object$call, table = tab, degree = object$degree, ncv = ns),
    class = "summary.cv_result")
  out
}


#' @method print summary.cv_result
#' @export
print.summary.cv_result <- function(x,
                                     digits = max(3, getOption("digits") - 3),
                                     ...) {
  rlang::check_dots_empty()
  cat("\nCall: ", deparse(x$call), fill = TRUE)
  cat("\n")
  cat("\nDegree of the estimated piecewise polynomial curve:", x$degree, "\n")
  cat("\nLambda =", x$tab$lambda[1], "gives the best CV score")
  cat("\n")
  cat("\nSummary of cross validation from", x$ncv, "lambdas:\n")
  print(x$tab, digits = digits)
  cat("\n")
}

#' Plot cv_result
#'
#' @param x result of cv_estimate_rt of class `cv_result`
#' @param plot_Rt plot cross-validation scores only if set to `FALSE`; plot
#' Rt and specify which Rt to generate with the `which_lambda` parameter if set
#' to `TRUE`
#' @param which_lambda select which Rt's to plot, if `plot_rt == TRUE`. If not
#' provided, plot R estimated with the optimal lambda. Any lambda provided must
#' match the values used in generating the Rt's.
#' @param ... Not used.
#'
#' @return plot of cv scores
#' @exportS3Method
#'
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
#' cv <- cv_estimate_rt(y, degree = 3, fold = 2, nsol=30)
#' plot(cv)
#' plot(cv, plot_Rt = TRUE)
plot.cv_result <- function(x, plot_Rt = FALSE, which_lambda = NULL, ...) {
  arg_is_numeric(which_lambda, allow_null = TRUE)
  arg_is_logical(plot_Rt)
  lambda <- x$lambda
  cv_scores <- x$cv_scores

  if (!plot_Rt && !is.null(which_lambda)) {
    cli::cli_warn("plot_Rt is False, ignoring input to `which_lambda`. Plot
                  cross validation score as default")
  }
  if (!plot_Rt) {
    df <- data.frame(
      cv_scores = cv_scores,
      lambda = lambda
    )
    plt_scores <- ggplot2::ggplot(df, ggplot2::aes(x=lambda, y = cv_scores))+
      ggplot2::geom_line() +
      ggplot2::theme_bw() +
      ggplot2::labs(title="Cross Validation Scores",
                    x = "Lambda", y = "CV scores")+
      ggplot2::scale_x_log10()
    return(plt_scores)

  } else {
    if (is.null(which_lambda))
      return(plot(x$full_fit, which_lambda = x$optimal_lambda))
    else
      return(plot(x$full_fit, which_lambda = which_lambda))
  }
}
