#' @method summary cv_poisson_rt
#' @export
summary.cv_poisson_rt <- function(object, ...) {

  rlang::check_dots_empty()

  tab <- with(object, data.frame(
    lambda = lambda,
    index = seq_along(lambda),
    cv_scores = cv_scores,
    cv_se = cv_se,
    dof = full_fit$dof
  ))
  n <- nrow(tab)
  if (n > 5) {
    l1 <- which(abs(object$lambda - object$lambda.min) < 1e-10)
    l2 <- which(abs(object$lambda - object$lambda.1se) < 1e-10)
    idx <- c(1, l1, l2, n)
    tab <- tab[idx, ]
    rownames(tab) <- c("Min Lambda", "CV Minimizer", "1se Lambda", "Max Lambda")
  }

  out <- structure(
    list(call = object$call, table = tab, degree = object$full_fit$degree),
    class = "summary.cv_poisson_rt")
  out
}


#' @method print summary.cv_poisson_rt
#' @export
print.summary.cv_poisson_rt <- function(
    x,
    digits = max(3, getOption("digits") - 3),
    ...) {

  rlang::check_dots_empty()

  lambda_warning = NULL
  if (x$table$index[2] == 1) lambda_warning = "smallest"
  if (x$table$index[2] == x$table$index[4]) lambda_warning = "largest"

  cat("\nCall:", deparse(x$call), fill = TRUE)
  cat("\nDegree of the estimated piecewise polynomial curve:", x$degree, "\n")
  if (!is.null(lambda_warning)) {
    cat("Warning: the CV minimum occurred at the", lambda_warning,
        "lambda in the path.\n\n")
  }
  cat("\nSummary of cross validation across lambda:\n")
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
#'  * If NULL, all estimated Rt values are plotted.
#'
#' @param ... Not used.
#'
#' @return plot of cv scores
#' @exportS3Method
#'
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
#' cv <- cv_estimate_rt(y, degree = 1, nfold = 3, nsol = 30)
#' plot(cv)
#' plot(cv, which_lambda = cv$lambda[1])
#' plot(cv, which_lambda = "lambda.min")
#' plot(cv, which_lambda = "lambda.1se")
#' plot(cv, NULL)
plot.cv_poisson_rt <- function(
    x, which_lambda = c("cv_scores", "lambda.min", "lambda.1se"), ...) {

  rlang::check_dots_empty()
  plt_scores <- FALSE
  if (is.character(which_lambda)) {
    which_lambda <- match.arg(which_lambda)
    if (which_lambda == "cv_scores") plt_scores <- TRUE
    else which_lambda <- x[[which_lambda]]
  } else {
    arg_is_numeric(which_lambda, allow_null = TRUE)
  }

  if (plt_scores) {
    df <- with(x, data.frame(
      cv_scores = cv_scores,
      lambda = lambda,
      cv_se = cv_se,
      upper = cv_scores + cv_se,
      lower = cv_scores - cv_se
    ))
    plt <- ggplot2::ggplot(df)+
      ggplot2::geom_errorbar(ggplot2::aes(x = .data$lambda,
                                          y = .data$cv_scores,
                                          ymin = .data$lower,
                                          ymax = .data$upper,
                                          width = 0.1))+
      ggplot2::geom_point(ggplot2::aes(x = .data$lambda, y = .data$cv_scores),
                          color="darkblue")+
      ggplot2::geom_vline(xintercept = x$lambda.min, linetype='dotted')+
      ggplot2::geom_vline(xintercept = x$lambda.1se, linetype='dotted')+
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Lambda", y = "CV scores")+
      ggplot2::scale_x_log10()
  } else {
    plt <- plot(x$full_fit, which_lambda = which_lambda)
  }

  return(plt)
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
#' cv <- cv_estimate_rt(y, degree = 3, nfold = 3, nsol = 30)
#' f <- fitted(cv)
#' f <- fitted(cv, which_lambda = cv$lambda[1])
#' f <- fitted(cv, which_lambda = "lambda.1se")
#' f <- fitted(cv, which_lambda = NULL)
fitted.cv_poisson_rt <- function(object,
                                 which_lambda = c("lambda.min", "lambda.1se"),
                                 ...) {
  rlang::check_dots_empty()
  if (is.character(which_lambda)) {
    which_lambda <- match.arg(which_lambda)
    which_lambda <- object[[which_lambda]]
  } else {
    arg_is_numeric(which_lambda, allow_null = TRUE)
  }
  fitted(object$full_fit, which_lambda)
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

