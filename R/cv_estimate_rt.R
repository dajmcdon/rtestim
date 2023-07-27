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

#' @return An object with S3 class `"cv_result"`. Among the list components:
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
cv_estimate_rt <- function(observed_counts,
                           korder = 3L,
                           dist_gamma = c(2.5, 2.5),
                           nfold = 3,
                           error_measure = c("mse", "mae", "deviance"),
                           x = 1:n,
                           lambda = NULL,
                           maxiter = 1e6L,
                           ...) {

  arg_is_pos_int(nfold)
  n <- length(observed_counts)
  arg_is_length(n, x)

  if (nfold == 1) cli::cli_abort("nfold must be greater than 1")

  ## Run program one time to create lambda
  full_data_fit <- estimate_rt(
    observed_counts = observed_counts,
    korder = korder,
    x = x,
    lambda = lambda,
    maxiter = maxiter,
    ...)

  if (is.null(lambda)) lambda <- full_data_fit$lambda

  foldid = fold_calculator(n, nfold)
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
    ### train test splits
    train_idx <- foldid != f
    test_idx <- foldid == f

    ## Run solver with the training set
    mod <- estimate_rt(
      observed_counts = observed_counts[train_idx],
      x = x[train_idx],
      korder = korder,
      lambda = lambda,
      maxiter = maxiter,
      ...)

    pred_rt <- stats::approx(
      x[train_idx], observed_counts[train_idx], x[test_idx])$y
    wpc <- delay_calculator(y[train_idx], x[train_idx])


    pred_observed_counts <- pred_rt * weighted_past_counts[test_idx]
    score <- colMeans(err_fun(observed_counts[test_idx], pred_observed_counts))
    cvall[f,] <- score
  }

  ### Calculate CV summary
  cv_scores <- colMeans(cvall)
  cv_se <- apply(cvall, FUN = stats::sd, MARGIN = 2)/sqrt(nfold)
  i0 <- which.min(cv_scores)

  structure(
    list(
      full_fit = full_data_fit,
      cv_scores = cv_scores,
      cv_se = cv_se,
      lambda = lambda,
      lambda.min = lambda[which.min(cv_scores)],
      lambda.1se = max(lambda[cv_scores <= cv_scores[i0] + cv_se[i0]]),
      call = match.call()
    ),
    class = "cv_poisson_rt"
  )
}



#' Helper function. Calculate the fold index for each `observed_counts`
#'
#' @inheritParams cv_estimate_rt
#' @param n length of the sequence to partition from
#' @return a vector of fold index at which the counts are distributed into folds
#'
fold_calculator <- function(n, nfold) {
  c(0, rep_len(1:nfold, n - 2), 0)
}


#' Helper function. Calculate Rt at hold-out set
#'
#' @param rt Matrix. Rt estimation at observed time points in training set
#' for all values of lambdas
#' @param n Integer. Number of total observations
#' @param train_idx vector of Integers. Index of `observed_counts` and `x` to be
#' assigned to the training set
#' @param test_idx vector of Integers. Index of `observed_counts` and `x` to be
#' assigned to the testing set
#' @param train_x vector of Integers. Location of the `observed_counts` in the
#' training set
#' @param test_x vector of Integers. Location of the `observed_counts` in the
#' testing set
#'
#' @return Predicted Rt at the hold-out set
pred_kth_rt <- function(rt, n, train_idx, test_idx, train_x, test_x) {
  ilo <- which((seq(1, n) %in% (test_idx - 1))[train_idx])
  ihi <- which((seq(1, n) %in% (test_idx + 1))[train_idx])
  a <- (test_x - train_x[ilo])/(train_x[ihi] - train_x[ilo])
  rt[ilo,] * (1 - a) + rt[ihi,] * a
}
