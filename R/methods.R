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
    lambda_idx <- match_lambda(which_lambda, lambda)
    Rt <- as.matrix(Rt[, lambda_idx])
    lambda <- lambda[lambda_idx]
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

#' @method summary cv_poisson_rt
#' @export
summary.cv_poisson_rt <- function(object, ...) {

  rlang::check_dots_empty()
  cv_scores <- object$cv_scores
  lambda <- object$lambda

  lambda_output <- c(min(lambda), object$optimal_lambda, object$lambda_1se,
                     max(lambda))
  lambda_idx <- match(lambda_output, lambda)

  tab <- with(object, data.frame(
    lambda = lambda_output,
    index = lambda_idx,
    cv_scores = cv_scores[lambda_idx],
    cv_se = object$cv_se[lambda_idx],
    dof = object$dof[lambda_idx])
    )
  rownames(tab) <- c("Min Lambda", "CV Minimizer", "1se Lambda", "Max Lambda")

  out <- structure(
    list(call = object$call, table = tab, degree = object$degree),
    class = "summary.cv_poisson_rt")
  out
}


#' @method print summary.cv_poisson_rt
#' @export
print.summary.cv_poisson_rt <- function(x,
                                     digits = max(3, getOption("digits") - 3),
                                     ...) {

  rlang::check_dots_empty()

  lambda_warning = NULL
  if (x$table$index[2] == 1) lambda_warning = "smallest"
  if (x$table$index[3] == x$table$index[4]) lambda_warning = "largest"

  cat("\nCall: ", deparse(x$call), fill = TRUE)
  cat("\n")
  cat("\nDegree of the estimated piecewise polynomial curve:", x$degree, "\n")
  cat("\n")
  if (!is.null(lambda_warning)) {
    cat("Warning: the CV minimum occurred at the", lambda_warning,
        "lambda in the path.\n\n")
  }
  cat("\nSummary of cross validation from", x$ncv, "lambdas:\n")
  print(x$tab, digits = digits)
  cat("\n")
}

#' Plot cv_poisson_rt
#'
#' @param x result of cv_estimate_rt of class `cv_poisson_rt`
#' @param plot_Rt plot cross-validation scores only if set to `FALSE`; plot
#' Rt and specify which Rt to generate with the `which_lambda` parameter if set
#' to `TRUE`
#' @param which_lambda select which Rt's to plot, if `plot_Rt == TRUE`. If not
#' provided, plot R estimated with the optimal lambda. Any lambda provided must
#' match the values used in generating the Rt's.
#' @param ... Not used.
#'
#' @return plot of cv scores
#' @exportS3Method
#'
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
#' cv <- cv_estimate_rt(y, degree = 3, nfold = 2, nsol=30)
#' plot(cv)
#' plot(cv, plot_Rt = TRUE)
#' plot(cv, plot_Rt = TRUE, which_lambda = cv$lambda[1])
plot.cv_poisson_rt <- function(x, plot_Rt = FALSE, which_lambda = NULL, ...) {

  arg_is_numeric(which_lambda, allow_null = TRUE)
  arg_is_logical(plot_Rt)

  lambda <- x$lambda
  lambda_1se <- x$lambda_1se
  optimal_lambda <- x$optimal_lambda
  cv_scores <- x$cv_scores
  cv_se <- x$cv_se

  cv_lambda_min <- cv_scores[which.min(lambda)]
  cv_lambda_1se <- cv_scores[match(lambda_1se, lambda)]

  if (!plot_Rt && !is.null(which_lambda)) {
    cli::cli_warn("plot_Rt is False, ignoring input to `which_lambda`. Plot
                  cross validation score as default")
  }

  if (!plot_Rt) {
    df <- data.frame(
      cv_scores = cv_scores,
      lambda = lambda,
      upper = cv_scores + x$cv_se,
      lower = cv_scores - x$cv_se,
      cv_se = cv_se
    )

    plt_scores <- ggplot2::ggplot(df, ggplot2::aes(x=lambda, y = cv_scores))+
      ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper,
                                          width = 0.1))+
      ggplot2::geom_line() +
      ggplot2::geom_point(color="darkblue")+
      ggplot2::geom_vline(xintercept = optimal_lambda, linetype='dotted')+
      ggplot2::geom_vline(xintercept = lambda_1se, linetype='dotted')+
      ggplot2::annotate("text", x = optimal_lambda,
                        y = c(max(cv_scores) + 1.1 * max(cv_se)),
                        label = "CV minimizer")+
      ggplot2::annotate("text", x = lambda_1se,
                        y = c(max(cv_scores) + 1.3 * max(cv_se)),
                        label = "1se Lambda")+
      ggplot2::theme_bw() +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(title="Cross Validation Scores",
                    x = "Lambda", y = "CV scores")+
      ggplot2::ylim(c(min(cv_scores) - 2* max(cv_se)),
                    max(cv_scores) + 2*max(cv_se))+
      ggplot2::scale_x_log10()

    return(plt_scores)

  } else {
    if (is.null(which_lambda))
      return(plot(x$full_fit, which_lambda = x$optimal_lambda))
    else
      return(plot(x$full_fit, which_lambda = which_lambda))
  }
}

fitted.cv_poisson_rt <- function(object, which_lambda, ...) {
  arg_is_numeric_or_char(which_lambda)
  rlang::check_dots_empty()

  full_fit <- object$full_fit
  lambda <- object$lambda
  lambda_idx <- match_lambda(which_lambda, lambda)

  if (is.numeric(which_lambda)) {
    full_fit$Rt[match(which_lambda, lambda)]
  }
  optimal_Rt <- object$optimal_lambda * object$weighted_past_counts
}


#' Find indices of elements in smaller list in a bigger list
#'
#' @param lambda list where the indices of elements are to be found from
#' @param which_lambda list of elements whose indices are to be found from `lambda`
#'
#' @return list of index of the elements in `which_lambda` from `lambda`
#' @export
#'
#' @examples
#' lambda <- c(1:20)
#' which_lambda <- c(1,3,6)
#' which_lambda_notin(lambda, which_lambda)
match_lambda <- function(which_lambda, lambda) {
  lambda_idx <- match(which_lambda, lambda)
  n <- length(lambda_idx)
  na_loc <- which(is.na(lambda_idx))

  if (length(na_loc) > 0) {
    if (length(na_loc) == n) {
      cli::cli_abort("No lambda is used to generate Rt in `estimate_rt()`")
    } else {
      warn_msg <- paste0("The [", toString(na_loc),
                         "]-th lambdas provided are not used to generate Rt",
                         " in `estimate_rt()`")
      cli::cli_warn(warn_msg)
      lambda_idx <- lambda_idx[-na_loc]
    }
  }
  return(lambda_idx)
}


