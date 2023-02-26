#' Leave-kth-out cross validation for choosing a optimal parameter lambda
#'
#' @inheritParams estimate_rt
#' @param fold Integer. This the number of folds to conduct the leave-kth-out
#' cross validation. For leave-kth-out cross validation, every kth
#' observed_counts and their corresponding position (evenly or unevenly
#' spaced) are placed into the same fold. The first and last observed_counts are
#' not assigned to any folds. Smallest allowable value if `fold = 2`. It is
#' generally not recommended to set `fold` to a large number

#' @return An object with S3 class `"cv_result"`. Among the list components:
#' * `cv_scores` leave-kth-out cross validation scores
#' * `optimal_lambda` lambda which achieved the optimal cv_scores from  cross
#' validation
#' * `optimal_Rt` the estimated effective reproduction rate, estimated with the
#' optimal lambda
#' * `observed_counts` vector of the observed daily infection counts
#' * `weighted_past_counts` the weighted sum of past infections counts with
#'   corresponding serial interval functions (or its Gamma approximation) as
#'   weights
#'
#' @export
#'
#' @examples cv_estimate_rt(c(1:20), degree = 2, fold = 2, lambda = c(1:4))
cv_estimate_rt <- function(observed_counts,
                           degree = 3L,
                           dist_gamma = c(2.5, 2.5),
                           fold = 3,
                           x = NULL,
                           lambda = NULL,
                           lambdamin = NULL,
                           lambdamax = NULL,
                           lambda_min_ratio = 1e-4,
                           ...) {

  ## Run program once to create lambda

  full_data_fit <- estimate_rt(
    observed_counts = observed_counts,
    degree = degree,
    x = x,
    lambda = lambda,
    lambdamax = lambdamax,
    lambdamin = lambdamin,
    lambda_min_ratio = lambda_min_ratio,
    ...)

  ## Use values from the full data fit
  if (is.null(lambda)) lambda <- full_data_fit$lambda
  weighted_past_counts <- full_data_fit$weighted_past_counts
  x <- full_data_fit$x
  n <- length(observed_counts)


  # Cross validation (copied from glmgen cv.trendfilter)
  foldid = fold_calculator(n, fold)
  if (length(x) == 0) x <- 1:n
  cvall <- matrix(0, fold, length(lambda))

  for (f in 1:fold) {
    ### train test splits
    train_idx <- which(foldid != f)
    test_idx <- which(foldid == f)

    train_observed_counts <- observed_counts[train_idx]
    train_weighted_past_counts <- weighted_past_counts[train_idx]
    train_x <- x[train_idx]

    test_observed_counts <- observed_counts[test_idx]
    test_weighted_past_counts <- weighted_past_counts[test_idx]
    test_x <- x[test_idx]

    ### Run solver with the training set
    mod <- estimate_rt(
      observed_counts = train_observed_counts,
      weighted_past_counts = train_weighted_past_counts,
      x = train_x,
      degree = degree,
      lambda = lambda,
      lambdamin = lambdamin,
      lambdamax = lambdamax,
      lambda_min_ratio = lambda_min_ratio,
      ...)


    ### Predict training value ###
    pred_rt <- pred_kth_rt(mod$Rt,
                           n = n,
                           train_idx = train_idx,
                           test_idx = test_idx,
                           train_x = train_x,
                           test_x = test_x)

    ### Calculate CV scores (Poisson loss) and store ###
    pred_observed_counts <- pred_rt*test_weighted_past_counts
    # taking neg-log of density without the factorial term as score
    score <- -colMeans(test_observed_counts*log(pred_observed_counts)
                       - pred_observed_counts)
    cvall[f,] <- score # scores for all lambda for this left-out set
  }

  ### Notes: cvall is a k by n_lambda matrix, each row contains score for each
  ### lambda for the particular fold

  ### Calculate CV scores
  cv_scores <- colMeans(cvall)
  op_lambda <- lambda[which.min(cv_scores)]


  ### Re-train with optimal lambda from cv
  op_mod <- estimate_rt(
    observed_counts = observed_counts,
    weighted_past_counts = weighted_past_counts,
    x = x,
    degree = degree,
    lambda = op_lambda,
    lambdamin = lambdamin,
    lambdamax = lambdamax,
    lambda_min_ratio = lambda_min_ratio,
    ...)


  structure(
    list(
      weighted_past_counts = weighted_past_counts,
      observed_counts = observed_counts,
      cv_scores = cv_scores,
      optimal_Rt = op_mod$Rt,
      optimal_lambda = op_lambda,
      lambda = lambda,
      x = x
    ),
    class = "cv_result"
  )
}



#' Helper function. Calculate the fold index for each `observed_counts`
#'
#' @inheritParams cv_estimate_rt
#' @return a vector of fold index at which the counts are distributed into folds
#' @export
#'
#' @examples fold_calculator(20, 3)
fold_calculator <- function(n, fold) {
  c(0,rep(seq(1,fold),n-2)[seq(1,n-2)],0)
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
#' @export
#'
#' @examples
#' n <- 10
#' x <- 1:10
#' train_idx <- c(1, 3, 5, 7, 9, 10)
#' test_idx <- c(2, 4, 6, 8)
#' train_x <- x[train_idx]
#' test_x <- x[test_idx]
#' rt <- matrix(c(1.1, 1.3, 1.5, 1.7, 1.9, 2.0), nrow=6)
#'
#' pred_kth_rt(rt, n, train_idx, test_idx, train_x, test_x
#' # Should equal to c(1.2, 1.4, 1.6, 1.8)
pred_kth_rt <- function(rt, n, train_idx, test_idx, train_x, test_x) {
  ilo <- which((seq(1,n)%in%(test_idx-1))[train_idx])
  ihi <- which((seq(1,n)%in%(test_idx+1))[train_idx])
  a <- (test_x - train_x[ilo])/(train_x[ihi] - train_x[ilo])
  rt[ilo,]*(1-a) + rt[ihi,]*a
}
