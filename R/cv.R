#' Leave-kth-out cross validation for choosing a optimal parameter lambda
#'
#' @param observed_counts vector of the observed daily infection counts
#' @param degree Integer. Degree of the piecewise polynomial curve to be
#'   estimated. Ror example, `degree = 0` corresponds to a piecewise constant
#'   curve.
#' @param dist_gamma Vector of length 2. These are the shape and scale for the
#'   assumed serial interval distribution. Roughly, this distribution describes
#'   the probability of an infectious individual infecting someone else after
#'   some period of time after having become infectious.
#'   As in most literature, we assume that this interval follows a gamma
#'   distribution with some shape and scale.
#' @param fold Integer. This the number of folds to conduct the leave-kth-out
#' cross validation. For leave-kth-out cross validation, every kth
#' observed_counts and their corresponding position (can be evenly or unevenly
#' spaced) are placed into the same fold. The first and last observed_counts are
#' not assigned to any folds.
#' @param x a vector of positions at which the counts have been observed. In an
#'   ideal case, we would observe data at regular intervals (e.g. daily or
#'   weekly) but this may not always be the case.
#' @param lambda Vector. A user supplied sequence of tuning parameters which
#'   determines the balance between data fidelity and
#'   smoothness of the estimated Rt; larger `lambda` results in a smoother
#'   estimate. The default, `NULL`
#'   results in an automatic computation based on `nlambda`, the largest value
#'   of `lambda` that would a maximally smooth estimate, and `lambdamin_ratio`.
#'   Supplying a value of `lambda` overrides
#'   this behaviour. It is likely better to supply a
#'   decreasing sequence of `lambda` values than a single (small) value. If
#'   supplied, the user-defined `lambda` sequence is automatically sorted in
#'   decreasing order.
#' @param nlambda Integer. The number of tuning parameters `lambda` at which to
#'   compute Rt.
#' @param lambdamin Optional value for the smallest `lambda` to use. This should
#'   be greater than zero.
#' @param lambdamax Optional value for the largest `lambda` to use.
#' @param lambda_min_ratio If neither `lambda` nor `lambdamin` is specified, the
#'   program will generate a lambdamin by lambdamax * lambda_min_ratio
#'   A multiplicative factor for the minimal lambda in the
#'   `lambda` sequence, where `lambdamin = lambdamin_ratio * lambdamax`.
#'   A very small value will lead to the solution `Rt = log(observed_counts)`.
#'   This argument has no effect if there is user-defined `lambda` sequence.
#' @param maxiter Integer. Maximum number of iterations for the estimation
#'   algorithm.
#' @param init a list of internal configuration parameters of class
#'   `rt_admm_configuration`.
#'
#' @return An object with S3 class `"cv_result"`. Among the list components:
#' * `cv_scores` leave-kth-out cross validation score
#' * `optimal_lambda` optimal lambda chosen from the cross validation
#' * `optimal_Rt` the estimated effective reproduction rate, estimated with the
#' optimal lambda
#'
#' @export
#'
#' @examples TODO: add this later
cv_estimate_rt <- function(observed_counts,
                           degree = 3L,
                           dist_gamma = c(2.5, 2.5),
                           fold = 3,
                           x = NULL,
                           lambda = NULL,
                           nlambda = 100L,
                           lambdamin = NULL,
                           lambdamax = NULL,
                           lambda_min_ratio = 1e-4,
                           maxiter = 1e4,
                           init = NULL) {

  # create weighted past case counts
  # TODO: weighted_past_counts here is incorrect. During cross validation, any
  # info regarding the test sets should not be contained in the training set.
  # weighted_past_counts for day t contains case count for day t-1 which could
  # be in the test fold.
  # Using this for not but need to change later
  weighted_past_counts <- delay_calculator(observed_counts, x, dist_gamma)
  if (is.null(init))
    init <- rt_admm_configuration(observed_counts, degree, weighted_past_counts)
  if (!inherits(init, "rt_admm_configuration"))
    cli::cli_abort("`init` must be created with `rt_admm_configuration()`.")
  if (is.null(init$primal_var)) {
    init <- rt_admm_configuration(
      observed_counts, init$degree, weighted_past_counts,
      auxi_var = init$auxi_var, dual_var = init$dual_var)
  }
  # validate maxiter
  maxiter <- as.integer(maxiter)
  n <- length(observed_counts)

  # (1) check that counts are non-negative, integer, vector
  if (any(observed_counts < 0)) cli::cli_abort("`observed_counts` must be non-negative")
  if (!is.vector(observed_counts)) cli::cli_abort("observed_counts must be a vector")

  # (2) checks on lambda, lambdamin, lambdamax
  lambda_size <- length(lambda)
  if (lambda_size > 0) {
    nsol <- lambda_size
    lambdamin <- min(lambda)
    lambdamax <- max(lambda)
  } else {
    msg <- "If lambda is not specified,"
    if (is.null(lambdamax))
      cli::cli_abort("{msg} lambdamax must be specified")
    if (lambda_min_ratio < 0 || lambda_min_ratio > 1)
      cli::cli_abort("{msg} lambda_min_ratio must be in [0,1]")
    if (lambdamax < 0)
      cli::cli_abort("{msg} lambdamax must be positive.")
    if (lambdamin < 0 && !is.null(lambdamin))
      cli::cli_abort("{msg} lambdamin must be positive.")
    if (lambdamin >= lambdamax && !is.null(lambdamin))
      cli::cli_abort("{msg} lambdamin must be < lambdamax.")

    # TODO: need to create lambda if not specified.
  }


  # (3) check that x is a double vector of length 0 or n
  if (!is.double(x)) x = as.double(x)
  if (!is.numeric(x)) cli::cli_abort("x must be a numeric vector")
  if (!(length(x) == n | length(x) == 0))
    cli::cli_abort("x must be length 0 or n")


  # Cross validation (copied from glmgen cv.trendfilter)
  foldid = fold_calculator(n, fold)
  if (length(x) == 0) x <- 1:n
  cvall <- matrix(0,fold,length(lambda))

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
    mod <- rtestim_path(
      train_observed_counts,
      train_x,
      train_weighted_past_counts,
      init$degree,
      lambda = lambda,
      nsol = nsol,
      rho = init$rho,
      maxiter = maxiter,
      tolerance = init$tolerance,
      verbose = init$verbose)

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

  ### Notes: cvall is a k by n_lambda matrix, each row contains score for each lambda
  ### for the particular fold

  ### Calculate CV scores
  cv_scores <- colMeans(cvall)
  op_lambda <- lambda[which.min(cv_scores)]



  ### Re-train with optimal lambda from cv
  op_mod <- rtestim_path(
    observed_counts,
    x,
    weighted_past_counts,
    init$degree,
    lambda = op_lambda,
    rho = init$rho,
    maxiter = maxiter,
    tolerance = init$tolerance,
    verbose = init$verbose)


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
#' @param n Integer. Total number of `observed_counts`
#' @param fold Integer. Number of folds for cross validation. In cross validation,
#' data are distributed into folds. Each fold during cross-validation is left as
#' the left-out data to be tested on after the model is fitted with the rest of
#' the set.
#'
#' @return a vector of fold index at which the counts are distributed into folds
#' @export
#'
#' @examples
fold_calculator <- function(n, fold) {
  c(0,rep(seq(1,fold),n-2)[seq(1,n-2)],0)
}


#' Helper function. Calculate Rt at hold-out fold
#'
#' @param rt an nlambda by (n-n_holdout) matrix. Rt estimated based on all lambdas for a
#' particular fold.
#' @param n
#' @param train_idx
#' @param test_idx
#' @param train_x
#' @param test_x
#'
#' @return
#' @export
#'
#' @examples
pred_kth_rt <- function(rt, n, train_idx, test_idx, train_x, test_x) {
  ilo <- which((seq(1,n)%in%(test_idx-1))[train_idx])
  ihi <- which((seq(1,n)%in%(test_idx+1))[train_idx])
  a <- (test_x - train_x[ilo])/(train_x[ihi] - train_x[ilo])
  rt[ilo,]*(1-a) + rt[ihi,]*a
}











