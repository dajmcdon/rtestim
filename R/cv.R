#' CV for admm
#'
#' @param observed_counts
#' @param degree
#' @param dist_gamma
#' @param fold
#' @param x
#' @param lambda
#' @param nsol
#' @param lambdamin
#' @param lambdamax
#' @param lambda_min_ratio
#' @param maxiter
#' @param init
#'
#' @return
#' @export
#'
#' @examples
cv_estimate_rt <- function(observed_counts,
                           degree = 3L,
                           dist_gamma = c(2.5, 2.5),
                           fold = 3,
                           x = NULL,
                           lambda = NULL,
                           nsol = 100L,
                           lambdamin = NULL,
                           lambdamax = NULL,
                           lambda_min_ratio = 1e-4,
                           maxiter = 1e4,
                           init = NULL) {

  # create weighted past cases
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
  # if (!all(rlang::is_intergerish(observed_counts))) not required
  #  cli::cli_abort("`observed_counts` must be integers")


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
  }


  # (3) check that x is a double vector of length 0 or n
  if (!is.double(x)) x = as.double(x)
  if (!is.numeric(x)) cli::cli_abort("x must be a numeric vector")
  if (!(length(x) == n | length(x) == 0))
    cli::cli_abort("x must be length 0 or n")


  # Cross validation (copied from glmgen cv.trendfilter)
  foldid = c(0,rep(seq(1,fold),n-2)[seq(1,n-2)],0)
  if (length(x) == 0) x <- 1:n
  cvall = matrix(0,fold,length(lambda))

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

    ### Run solver
    mod <- rtestim_path(
      train_observed_counts,
      train_x,
      train_weighted_past_counts,
      init$degree,
      lambda = lambda,
      lambdamax = lambdamax,
      lambdamin = lambdamin,
      nsol = nsol,
      rho = init$rho,
      maxiter = maxiter,
      tolerance = init$tolerance,
      lambda_min_ratio = lambda_min_ratio,
      verbose = init$verbose)

    ### Predict training value ###
    ilo <- which((seq(1,n)%in%(test_idx-1))[train_idx])
    ihi <- which((seq(1,n)%in%(test_idx+1))[train_idx])
    a <- (test_x - train_x[ilo])/(train_x[ihi] - train_x[ilo])
    pred_rt <- mod$Rt[ilo,]*(1-a) + mod$Rt[ihi,]*a
    score <- colMeans((test_observed_counts - (pred_rt*test_weighted_past_counts))^2)
    cvall[f,] <- score
  }
  output <- list()
  cv_scores <- colMeans(cvall)
  output$cv_scores <- cv_scores
  op_lambda <- lambda[which.min(cv_scores)]
  output$optimal_lambda <- op_lambda

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

  output$optimal_rt <- op_mod$Rt
  output$pred_count <- op_mod$Rt*weighted_past_counts

  output

}

