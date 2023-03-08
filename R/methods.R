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
    Rt <- Rt[, lambda_idx, drop = FALSE]
    lambda <- lambda[lambda_idx]
  }

  df <- data.frame(
    Rt = c(Rt),
    lambda = rep(lambda, each = nrow(Rt)),
    Time = rep(x$x, ncol(Rt))
  )

  plt <- ggplot2::ggplot(
    df,
    ggplot2::aes(.data$Time, .data$Rt, colour = .data$lambda,
                 group = .data$lambda)) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_viridis_c(trans = "log10")
  if (ncol(Rt) == 1) plt <- plt + ggplot2::theme(legend.position= "none")

  print(plt)
}

#' @method summary cv_poisson_rt
#' @export
summary.cv_poisson_rt <- function(object, ...) {

  rlang::check_dots_empty()
  cv_scores <- object$cv_scores
  lambda <- object$lambda

  lambda_output <- c(min(lambda), object$lambda.min, object$lambda.1se,
                     max(lambda))
  lambda_idx <- almost_match(lambda_output, lambda)
  print(object$lambda.min)
  print(object$lambda.1se)
  tab <- with(object, data.frame(
    lambda = lambda_output,
    index = lambda_idx,
    cv_scores = cv_scores[lambda_idx],
    cv_se = object$cv_se[lambda_idx],
    dof = object$full_fit$dof[lambda_idx])
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
  if (x$table$index[2] == x$table$index[4]) lambda_warning = "largest"

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
#' @param which_lambda select which Rt's to plot.
#'
#' If not provided, the cross validation score will be plotted. If provided a
#' list of lambda, the corresponding Rt estimation will be plotted.
#'
#' If provided a string, it
#' must be either one of `lambda.min`, `lambda.1se`, or `cv_scores`.
#'
#'  * If provided `lambda.min`, plot Rt which is generated from the lambda that
#'  minimizes the cross validation score.
#'
#'  * If provided `lambda.1se`, plot Rt which is generated from the lambda whose
#'  corresponding cross validation score is 1 standard error away of the
#'  minimal cross validation score.
#'
#'  * If provided `cv_scores`, plot the cross validation score.
#'
#' @param ... Not used.
#'
#' @return plot of cv scores
#' @exportS3Method
#'
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
#' cv <- cv_estimate_rt(y, degree = 1, nfold = 2, nsol = 30)
#' plot(cv)
#' plot(cv, which_lambda = cv$lambda[1])
#' plot(cv, which_lambda = "lambda.min")
#' plot(cv, which_lambda = "lambda.1se")
plot.cv_poisson_rt <- function(x,
                               which_lambda = c("cv_scores",
                                                "lambda.min",
                                                "lambda.1se"), ...) {

  rlang::check_dots_empty()
  if (is.character(which_lambda))
    which_lambda <- match.arg(which_lambda)
  else arg_is_numeric(which_lambda, allow_null = TRUE)

  lambda.1se <- x$lambda.1se
  lambda.min <- x$lambda.min

  print(which_lambda)

  if (is.numeric(which_lambda)) {
    return(plot(x$full_fit, which_lambda = which_lambda))

  } else if (which_lambda == "cv_scores" || all(is.null(which_lambda))) {
      df <- data.frame(
        cv_scores = x$cv_scores,
        lambda = x$lambda,
        cv_se = x$cv_se,
        upper = x$cv_scores + x$cv_se,
        lower = x$cv_scores - x$cv_se
      )

      plt_scores <- ggplot2::ggplot(df)+
        ggplot2::geom_errorbar(ggplot2::aes(x = .data$lambda,
                                            y = .data$cv_scores,
                                            ymin = .data$lower,
                                            ymax = .data$upper,
                                            width = 0.1))+
        ggplot2::geom_point(ggplot2::aes(x = .data$lambda, y = .data$cv_scores),
                            color="darkblue")+
        ggplot2::geom_vline(xintercept = lambda.min, linetype='dotted')+
        ggplot2::geom_vline(xintercept = lambda.1se, linetype='dotted')+
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Lambda", y = "CV scores")+
        ggplot2::scale_x_log10()

      return(plt_scores)

    } else if (which_lambda == "lambda.1se") {
        return(plot(x$full_fit, which_lambda = x$lambda.1se))

    } else if (which_lambda == "lambda.min") {
        return(plot(x$full_fit, which_lambda = x$lambda.min))
  }
}


#' Fitted cv_poisson_rt
#'
#' @param object result of cross validation of type `cv_poisson_rt`
#' @param which_lambda select which Rt's to output. If not provided, all Rt's
#' are returned. If provided a list of lambda,the corresponding Rt estimation
#' will be returned.
#'
#' If provided a string, it must be either one of `lambda.min` or `lambda.1se`.
#'
#'  * If provided `lambda.min`, return Rt which is generated from
#'  the lambda that minimizes the cross validation score.
#'
#'  * If provided `lambda.1se`, return Rt which is generated from the lambda
#'  whose corresponding cross validation score is 1 standard error away of the
#'  minimal cross validation score.
#' @param ... not used.
#'
#' @return Rt's estimated from provided lambda
#' @exportS3Method
#'
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
#' cv <- cv_estimate_rt(y, degree = 3, nfold = 2, nsol = 30)
#' f <- fitted(cv)
#' f <- fitted(cv, which_lambda = cv$lambda[1])
#' f <- fitted(cv, which_lambda = "lambda.min")
#' f <- fitted(cv, which_lambda = "lambda.1se")
fitted.cv_poisson_rt <- function(object, which_lambda = NULL, ...) {
  if (is.character(which_lambda))
    which_lambda <- match.arg(which_lambda,
                              c("lambda.min", "lambda.1se"))
  else arg_is_numeric(which_lambda, allow_null = TRUE)
  rlang::check_dots_empty()

  full_fit <- object$full_fit
  lambda <- object$lambda

  if (is.numeric(which_lambda)) {
    lambda_idx <- match_lambda(which_lambda, lambda)
    return(full_fit$Rt[, almost_match(which_lambda, lambda)])

  } else if (is.null(which_lambda)) {
    return(full_fit$Rt)

  } else if (which_lambda == "lambda.min") {
    return(full_fit$Rt[, which.min(object$cv_scores)])

  } else if (which_lambda == "lambda.1se") {
    return(full_fit$Rt[, almost_match(object$lambda.1se, lambda)])
  }
}


#' Find indices of a list of elements from another list
#'
#' Find indices of every elements from `which_lambda` in `lambda`. Return
#' the indices if all elements in `which_lambda` are in `lambda`. Return
#' indices and warning if some elements in `which_lambda` are not in `lambda`.
#' Abort the program if no element in `which_lambda` is in `lambda`
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
#' match_lambda(lambda, which_lambda)
match_lambda <- function(which_lambda, lambda) {
  lambda_idx <- almost_match(which_lambda, lambda)
  n <- length(lambda_idx)
  no_match_loc <- which(lambda_idx == 0)

  if (length(no_match_loc) > 0) {
    if (length(no_match_loc) == n) {
      cli::cli_abort("No lambda is used to generate Rt in `estimate_rt()`")
    } else {
      warn_msg <- paste0("The [", toString(no_match_loc),
                         "]-th lambdas provided are not used to generate Rt",
                         " in `estimate_rt()`")
      cli::cli_warn(warn_msg)
      lambda_idx <- lambda_idx[-no_match_loc]
    }
  }
  return(lambda_idx)
}

