#' Add confidence bands to estimated Rt or incidence curves
#'
#' Create an approximate confidence band for the Rt or incidence estimate. Note
#' that the variance computation is approximate.
#'
#' @param object a `poisson_rt` or `cv_poisson_rt` object.
#' @param lambda the selected lambda. May be a scalar value, or in the case of
#'  `cv_poisson_rt` objects, `"lambda.min"` or `"lambda.max"`.
#' @param level the desired confidence level(s). These will be sorted if
#'   necessary.
#' @param type the type `Rt` or `Yt` for confidence intervals of fitted
#'   Rt or fitted incident cases
#' @param ... additional arguments for methods. Unused.
#'
#' @return A `data.frame` containing the estimates `Rt` or `Yt` at the chosen
#'  `lambda`, and confidence limits corresponding to `level`
#' @export
#'
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
#' out <- estimate_rt(y, nsol = 10)
#' head(confband(out, out$lambda[2]))
#' head(confband(out, out$lambda[2], level = c(0.95, 0.8, 0.5)))
#'
#' cv <- cv_estimate_rt(y, nfold = 3, nsol = 30)
#' head(confband(cv, "lambda.min", c(0.5, 0.9)))
confband <- function(object, lambda, level = 0.95, type = c("Rt", "Yt"), ...) {
  UseMethod("confband")
}

#' @export
confband.cv_poisson_rt <- function(
    object,
    lambda = c("lambda.min", "lambda.1se"),
    level = 0.95,
    type = c("Rt", "Yt"), ...) {
  rlang::check_dots_empty()
  assert_numeric(level, lower = 0, upper = 1)
  type <- rlang::arg_match(type)
  if (is.character(lambda)) {
    lambda <- object[[rlang::arg_match(lambda)]]
  } else {
    assert_numeric(lambda, len = 1L)
  }
  confband(object$full_fit, lambda = lambda, level = level, type = type)
}

#' @export
confband.poisson_rt <- function(object, lambda, level = 0.95, type = c("Rt", "Yt"), ...) {
  rlang::check_dots_empty()
  assert_numeric(lambda, len = 1L)
  assert_numeric(level, lower = 0, upper = 1)
  level <- sort(level, decreasing = TRUE)
  type <- rlang::arg_match(type)

  nbd <- function(piece, ord) {
    n <- length(piece)
    if (n <= ord + 1) {
      return(Matrix::Diagonal(n, x = 1))
    }
    get_D(ord, piece)
  }

  y <- object$observed_counts
  n <- length(y)
  Rt <- fitted(object, lambda)
  kn <- find_knots(object, lambda)
  wt <- object$weighted_past_counts
  # wt <- map2(kn$l, kn$r, function(a, b) object$weighted_past_counts[a:b])
  yhat <- predict(object, lambda)
  # yhat <- map2(kn$l, kn$r, function(a, b) yhat[a:b])

  # Ds <- lapply(kn$xpieces, nbd, ord = object$korder)
  # browser()
  # The procedure for this approximation:
  # 00. theta is natural parameter in exp family
  # 0. Pretend we knew the knots (and lambda is fixed + known) --> Ds
  # 1. pretend we had used ||Ds theta||_2^2 with the same lambda and normal
  #    likelihood.
  # 2. Apply multivariate delta method
  #   a. Var(\hat\mu) = (I * 1/theta^2 + lambda (Ds'Ds))^(-1), evaluated at
  #      \theta = \hat\mu
  #   b. g(mu) = mu / w
  #   c. g'(u) = 1 / w
  #   d. Need Var(\hat\mu) * g'(mu)^2 -->
  # Do it blockwise, try the easy way, use ginv if this fails
  # covs <- pmap(list(Ds, wt, yhat), function(D, w, yh) {
  #   n <- length(w)
  #   kernel <- Matrix::Diagonal(n, x = 1 / yh^2) + lambda * Matrix::crossprod(D)
  #   t1 <- Matrix::diag(Matrix::solve(kernel, Matrix::Diagonal(n, x = 1/w^2)))
  #   if (any(t1 < 0)) t1 <- diag(MASS::ginv(as.matrix(kernel))) / w^2
  #   t1
  # })
  # covs <- unlist(covs)
  # covs <- Matrix::diag(Matrix::solve(
  #   Matrix::Diagonal(n, 1 / yhat^2) + lambda * Matrix::crossprod(Ds)
  # )) / wt^2
  #
  D <- get_D(object$korder, object$x)
  kernel <- Matrix::Diagonal(x = 1 / yhat^2) + lambda * Matrix::crossprod(D)
  covs <- Matrix::diag(Matrix::solve(kernel))
  a <- (1 - level) / 2
  a <- c(a, rev(1 - a))

  if (type == "Rt") {
    covs <- covs / wt^2
    fit <- Rt
  } else if (type == "Yt") fit <- yhat
  cb <- outer(sqrt(covs), stats::qt(a, n - kn$dof))
  cb <- pmax(fit + cb, 0)
  colnames(cb) <- fmt_perc(a)

  tibble::new_tibble(
    type = type,
    vctrs::vec_cbind(fit = fit, cb),
    lambda = lambda,
    CIs = level,
    dof = kn$dof,
    xval = object$x,
    class = "rt_confidence_band"
  )
}

fmt_perc <- function(probs, digits = 3) {
  paste0(
    format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
    "%"
  )
}

#' @exportS3Method print rt_confidence_band
print.rt_confidence_band <- function(x, ...) {
  cat("An `rt_confidence_band` object.\n\n")
  cat(paste("* type =", attr(x, "type"), "\n"))
  cat(paste("* lambda =", round(attr(x, "lambda"), 3), "\n"))
  cat(paste("* degrees of freedom =", attr(x, "dof"), "\n"))
  cat("\n")
  NextMethod()
}

#' Plot estimated confidence bands for an estimate of Rt
#'
#' Produces a figure showing a single estimated Rt value along with approximate
#' confidence bands. The result is a [ggplot2::ggplot()]. Additional user
#' modifications can be added as desired.
#'
#' @param x An object of class `rt_confidence_band` as produced by `confband()`.
#' @param colour The colour of the desired plot
#' @param ... Not used.
#'
#'
#' @exportS3Method plot rt_confidence_band
#' @examples
#'
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
#' out <- estimate_rt(y, nsol = 10)
#' cb <- confband(out, out$lambda[2], level = c(0.95, 0.8, 0.5))
#' plot(cb)
#' cb_y <- confband(out, out$lambda[2], level = c(0.95, 0.8, 0.5), type = "Yt")
#' plot(cb_y)
plot.rt_confidence_band <- function(x, colour = "#3A448F", ...) {
  x$xval <- attr(x, "xval")
  CIs <- names(x)[grep("[0-9]", names(x))]
  xlab <- ifelse(inherits(attr(x, "xval"), "Date"), "Date", "Time")
  ylab <- paste(
    "Estimated", attr(x, "type"), "with",
    paste(fmt_perc(rev(attr(x, "CIs"))), collapse = ", "),
    "\nconfidence bands"
  )
  plt <- ggplot2::ggplot(x, ggplot2::aes(x = .data$xval)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$fit), colour = colour) +
    ggplot2::theme_bw() +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)
  ci_plot <- plot_cis(plt, CIs, colour)
  if (attr(x, "type") == "Rt") {
    ci_plot <- ci_plot + ggplot2::geom_hline(yintercept = 1)
  }
  ci_plot
}

plot_cis <- function(plot, CIs, fill = "#3A448F",
                     alpha = 0.6, linewidth = 0.05) {
  n <- length(CIs) / 2
  alpha <- alpha / (n - 1)

  for (i in 1:n) {
    bottom <- CIs[i]
    top <- rev(CIs)[i]
    if (i == 1) {
      plot <- plot +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = .data[[bottom]], ymax = .data[[top]]),
          alpha = 0.2, linewidth = linewidth, fill = fill
        )
    } else {
      plot <- plot +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = .data[[bottom]], ymax = .data[[top]]),
          fill = fill, alpha = alpha
        )
    }
  }
  plot
}
