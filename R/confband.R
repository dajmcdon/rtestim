

#' Add confidence bands to estimated Rt curves
#'
#' Create an approximate confidence band for the Rt estimate. Note that the
#' variance computation is approximate.
#'
#' @param object a `poisson_rt` or `cv_poisson_rt` object.
#' @param lambda the selected lambda. May be a scalar value, or in the case of
#'  `cv_poisson_rt` objects, `"lambda.min"` or `"lambda.max"`.
#' @param level the desired confidence level(s). These will be sorted if
#'   necessary.
#' @param ... additional arguments for methods. Unused.
#'
#' @return A `data.frame` containing the estimated `Rt` at the chosen `lambda`,
#'  and confidence limits corresponding to `level`
#' @export
#'
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
#' out <- estimate_rt(y, nsol = 10)
#' head(confband(out, out$lambda[2]))
#' head(confband(out, out$lambda[2], level = c(0.95, 0.8, 0.5)))
#'
#' cv <- cv_estimate_rt(y, nfold = 3, nsol = 30)
#' head(confband(cv, "lambda.min", c(0.5, 0.9)))
confband <- function(object, lambda, level = 0.95, ...) {
  UseMethod("confband")
}

#' @export
confband.cv_poisson_rt <- function(
    object,
    lambda = c("lambda.min", "lambda.1se"),
    level = 0.95, ...) {
  rlang::check_dots_empty()
  arg_is_probabilities(level)
  if (is.character(lambda)) lambda <- object[[match.arg(lambda)]]
  else arg_is_numeric(lambda)
  confband(object$full_fit, lambda = lambda, level = level)
}

#' @export
confband.poisson_rt <- function(object, lambda, level = 0.95, ...) {
  rlang::check_dots_empty()
  arg_is_numeric(lambda)
  arg_is_scalar(lambda)
  arg_is_probabilities(level)
  level <- sort(level, decreasing = TRUE)

  nbd <- function(n, ord) {
    if (n < ord + 1) return(Matrix::Diagonal(n, x = 0))
    buildD(n, ord)
  }

  y <- object$observed_counts
  n <- length(y)
  Rt <- fitted(object, lambda)
  yhat <- predict(object, lambda)
  knots <- find_knots(object, lambda)
  Ds <- Matrix::bdiag(lapply(knots$lens, nbd, ord = object$degree))
  # The procedure for this approximation:
  # 0. Pretend we knew the knots (and lambda is fixed + known) --> Ds
  # 1. pretend we had used ||Ds r||_2^2 with the same lambda
  # 2. Apply multivariate delta method
  #   a. Var(y) = (I * exp(theta) + lambda (Ds'Ds))^(-1)
  #   b. g(theta) = exp(theta) / w
  #   c. g'(theta) = exp(theta) / w = r
  #   d. Need Var(y) * g'(theta)^2 -->
  covs <- diag(Matrix::solve(
    Matrix::Diagonal(n, yhat) + lambda * Matrix::crossprod(Ds)
  )) * Rt^2
  a <- (1 - level) / 2
  a <- c(a, rev(1 - a))
  cb <- outer(sqrt(covs), stats::qt(a, n - knots$dof))
  cb <- pmax(Rt + cb, 0)
  colnames(cb) <- fmt_perc(a)
  vctrs::vec_cbind(Rt = Rt, cb)
}

fmt_perc <- function(probs, digits = 3) {
  paste(
    format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
    "%")
}
