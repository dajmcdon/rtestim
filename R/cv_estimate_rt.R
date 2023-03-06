#' Leave-kth-out cross validation for choosing a optimal parameter lambda
#'
#' @inheritParams estimate_rt
#' @param nfold Integer. This the number of folds to conduct the leave-kth-out
#' cross validation. For leave-kth-out cross validation, every kth
#' observed_counts and their corresponding position (evenly or unevenly
#' spaced) are placed into the same fold. The first and last observed_counts are
#' not assigned to any folds. Smallest allowable value if `fold = 2`. It is
#' generally not recommended to set `fold` to a large number
#' @param ... not used.

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
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15)*500 + 1))
#' cv <- cv_estimate_rt(y, degree = 3, nfold = 2, nsol=30)
#' cv
cv_estimate_rt <- function(observed_counts,
                           degree = 3L,
                           dist_gamma = c(2.5, 2.5),
                           nfold = 3,
                           x = NULL,
                           lambda = NULL,
                           ...) {

  arg_is_pos_int(nfold)
  if (nfold==1) cli::cli_abort("nfold must be greater than 1")

  ## Run program one time to create lambda
  full_data_fit <- estimate_rt(
    observed_counts = observed_counts,
    degree = degree,
    x = x,
    lambda = lambda,
    ...)

  ## Use values from the full data fit
  if (is.null(lambda)) lambda <- full_data_fit$lambda
  weighted_past_counts <- full_data_fit$weighted_past_counts
  x <- full_data_fit$x
  n <- length(observed_counts)


  # Cross validation (copied from glmgen cv.trendfilter)
  foldid = fold_calculator(n, nfold)
  cvall <- matrix(0, nfold, length(lambda))

  for (f in 1:nfold) {
    ### train test splits
    train_idx <- which(foldid != f)
    test_idx <- which(foldid == f)

    ### Run solver with the training set
    mod <- estimate_rt(
      observed_counts = observed_counts[train_idx],
      x = x[train_idx],
      degree = degree,
      lambda = lambda,
      ...)

    ### Predict training value ###
    pred_rt <- pred_kth_rt(mod$Rt,
                           n = n,
                           train_idx = train_idx,
                           test_idx = test_idx,
                           train_x = x[train_idx],
                           test_x = x[test_idx])

    ### Calculate CV scores (Poisson loss) and store ###
    pred_observed_counts <- pred_rt * weighted_past_counts[test_idx]
    # taking neg-log of density without the factorial term as score
    score <- -colMeans(observed_counts[test_idx] * log(pred_observed_counts)
                       - pred_observed_counts)
    cvall[f,] <- score # scores for all lambda for this left-out set
  }

  ### Notes: cvall is a k by n_lambda matrix, each row contains score for each
  ### lambda for the particular fold

  ### Calculate CV summary
  cv_scores <- colMeans(cvall)
  # calculate se for all lambda
  cv_se <- apply(cvall, FUN = stats::sd, MARGIN = 2)/sqrt(nfold)
  i0 <- which.min(cv_scores)

  op_lambda <- lambda[which.min(cv_scores)]
  lambda_1se <- max(lambda[cv_scores <= cv_scores[i0] + cv_se[i0]])

  structure(
    list(
      full_fit = full_data_fit,
      weighted_past_counts = weighted_past_counts,
      observed_counts = observed_counts,
      cv_scores = cv_scores,
      cv_se = cv_se,
      optimal_lambda = op_lambda,
      lambda = lambda,
      lambda_1se = lambda_1se,
      x = x,
      call = match.call(),
      degree = degree,
      dof = full_data_fit$dof
    ),
    class = "cv_poisson_rt"
  )
}



#' Helper function. Calculate the fold index for each `observed_counts`
#'
#' @inheritParams cv_estimate_rt
#' @param n length of the sequence to partition from
#' @return a vector of fold index at which the counts are distributed into folds
#' @export
#'
#' @examples fold_calculator(20, 3)
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
#' @export
#'
#' @examples
#' x <- 1:10
#' n <- 10
#' train_idx <- c(1, 3, 5, 7, 9, 10)
#' test_idx <- c(2, 4, 6, 8)
#' train_x <- x[train_idx]
#' test_x <- x[test_idx]
#' rt <- matrix(c(1.1, 1.3, 1.5, 1.7, 1.9, 2.0), nrow=6)
#' pred_kth_rt(rt, n, train_idx, test_idx, train_x, test_x)
pred_kth_rt <- function(rt, n, train_idx, test_idx, train_x, test_x) {
  ilo <- which((seq(1, n) %in% (test_idx - 1))[train_idx])
  ihi <- which((seq(1, n) %in% (test_idx + 1))[train_idx])
  a <- (test_x - train_x[ilo])/(train_x[ihi] - train_x[ilo])
  rt[ilo,] * (1 - a) + rt[ihi,] * a
}
