#' Leave-kth-out cross validation for choosing a optimal parameter lambda
#'
#' @inheritParams estimate_rt
#' @param nfold Integer. This number of folds to conduct the leave-kth-out
#' cross validation. For leave-kth-out cross validation, every kth
#' observed_counts and their corresponding position (evenly or unevenly
#' spaced) are placed into the same fold. The first and last observed_counts are
#' not assigned to any folds. Smallest allowable value is `nfold = 2`.
#' @param error_measure Metric used to calculate cross validation scores.
#' Must be choose from `mse`, `mae`, and `deviance`.
#' `mse` calculates the mean square error; `mae` calculates the mean absolute error;
#' `deviance` calculates the deviance
#' @param ... additional parameters passed to `estimate_rt()` function

#' @return An object with S3 class `"cv_poisson_rt"`. Among the list components:
#' * `full_fit` An object with S3 class `"poisson_rt"`, fitted with all
#' `observed_counts` and `lambda`
#' * `cv_scores` leave-kth-out cross validation scores
#' * `cv_se` leave-kth-out cross validation standard error
#' * `lambda.min` lambda which achieved the optimal cross validation score
#' * `lambda.1se` lambda that gives the optimal cross validation score
#' within one standard error.
#' * `lambda` the value of `lambda` used in the algorithm.
#' @export
#'
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
#' cv <- cv_estimate_rt(y, korder = 3, nfold = 3, nsol = 30)
#' cv
cv_estimate_rt <- function(
    observed_counts,
    korder = 3L,
    dist_gamma = c(2.5, 2.5),
    nfold = 3L,
    error_measure = c("mse", "mae", "deviance"),
    x = 1:n,
    lambda = NULL,
    maxiter = 1e6L,
    delay_distn = NULL,
    ...) {

  arg_is_pos_int(nfold)
  n <- length(observed_counts)
  arg_is_length(n, x)
  xin <- x
  if (inherits(xin, "Date")) x <- as.numeric(x)
  arg_is_numeric(x)

  if (nfold == 1) cli_abort("nfold must be greater than 1")

  ## Run program one time to create lambda
  full_fit <- estimate_rt(
    observed_counts = observed_counts,
    korder = korder,
    dist_gamma = dist_gamma,
    x = xin,
    lambda = lambda,
    maxiter = maxiter,
    delay_distn = delay_distn,
    ...)

  if (is.null(lambda)) lambda <- full_fit$lambda

  foldid <- c(0, rep_len(1:nfold, n - 2), 0)
  cvall <- matrix(0, nfold, length(lambda))
  error_measure <- match.arg(error_measure)
  err_fun <- switch(
    error_measure,
    mse = function(y, m) (y - m)^2,
    mae = function(y, m) abs(y - m),
    deviance = function(y, m) {
      devr <- y * log(m) - m
      devy <- y * log(y) - y
      devy[y == 0] <- 0
      2 * (devr - devy) }
  )

  for (f in 1:nfold) {
    train_idx <- foldid != f
    test_idx <- foldid == f

    mod <- estimate_rt(
      observed_counts = observed_counts[train_idx],
      x = x[train_idx],
      dist_gamma = dist_gamma,
      korder = korder,
      lambda = lambda,
      maxiter = maxiter,
      delay_distn = delay_distn,
      ...)

    interp_rt <- interpolate_rt(mod, x[test_idx])

    wpc <- delay_calculator(
      observed_counts = observed_counts[train_idx],
      x = x[train_idx] - min(x[train_idx]) + 1,
      dist_gamma = dist_gamma,
      delay_distn = delay_distn,
      output_partial_seq = FALSE
    )


    pred_observed_counts <- interp_rt * wpc[test_idx]
    score <- colMeans(err_fun(observed_counts[test_idx], pred_observed_counts))
    tryCatch(cvall[f,] <- score,
             error = function(w) {
               print("Error in `cvall[f,] <- score`. Please increase the maximum iteration `maxiter`.")
               stop()
              })
  }

  ### Calculate CV summary
  cv_scores <- colMeans(cvall)
  cv_se <- apply(cvall, FUN = stats::sd, MARGIN = 2) / sqrt(nfold)
  i0 <- which.min(cv_scores)

  structure(
    enlist(
      full_fit,
      cv_scores,
      cv_se,
      lambda,
      lambda.min = lambda[i0],
      lambda.1se = max(
        lambda[cv_scores <= cv_scores[i0] + cv_se[i0]],
        na.rm = TRUE),
      call = match.call()
    ),
    class = "cv_poisson_rt"
  )
}
