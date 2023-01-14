#' ADMM initialization
#'
#' We should convert this into an S3 method so that we can pass in a
#' vector of counts or an admm_initializer object (and overwrite)
#'
#' @param current_counts the current daily infection counts
#' @param weighted_past_counts the weighted sum of past infection counts with
#' corresponding serial interval functions (or its Gamma approximation) as
#' weights
#' @param degree degree of the piecewise polynomial curves to be fitted,
#' e.g., degree = 0 corresponds to piecewise constant curves
#' @param primal_var the initial value of log(Rt)
#' @param auxi_var auxiliary variable in the ADMM algorithm
#' @param dual_var dual variable in the ADMM algorithm
#'
#' @return a list of model parameters with class `admm_initializer`
#'
#' @export
admm_initializer <- function(current_counts,
                             degree,
                             weighted_past_counts = NULL,
                             primal_var = NULL,
                             auxi_var = NULL,
                             dual_var = NULL,
                             rho = -1,
                             rho_adjust = -1,
                             tolerance = 1e-4,
                             verbose = 0) {
  n <- length(current_counts)
  degree <- as.integer(degree)

  if (degree < 0) cli::cli_abort("`degree` must be non-negative. ")
  if (length(rho) != 1 || !is.numeric(rho)) {
    cli::cli_warn(c("`rho` must be a numeric scalar.",
                    i = "Resetting to default value."))
    rho = -1
  }
  if (length(rho_adjust) != 1 || !is.numeric(rho_adjust)) {
    cli::cli_warn(c("`rho_adjust` must be a numeric scalar.",
                    i = "Resetting to default value."))
    rho_adjust = -1
  }
  if (length(tolerance) != 1 || !is.numeric(tolerance) || tolerance <= 0) {
    cli::cli_warn(c("`tolerance` must be a positive scalar.",
                    i = "Resetting to default value."))
    tolerance = 1e-4
  }
  if (length(verbose) != 1 || !is.numeric(verbose) || verbose < 0) {
    cli::cli_warn(c("`verbose` must be a non-negative scalar.",
                    i = "Resetting to default value."))
    verbose = 0
  }

  if (is.null(primal_var)) {
    if (!is.null(weighted_past_counts)) {
      # should we divide by n?
      # what do we do when current_counts == 0
      primal_var <- log(current_counts / (n * weighted_past_counts))
    }
  } else {
    if (length(primal_var) != n) {
      cli::cli_abort("`primal_var` must have length {n}.")
    }
  }
  if (is.null(auxi_var)) auxi_var <- double(n - degree)
  else {
    if (length(auxi_var) != n - degree) {
      cli::cli_abort("`auxi_var` must have length {n - degree}.")
    }
  }
  if (is.null(dual_var)) dual_var <- double(n - degree)
  else {
    if (length(dual_var) != n - degree) {
      cli::cli_abort("`dual_var` must have length {n - degree}.")
    }
  }

  structure(
    list(
      degree = degree,
      primal_var = primal_var,
      auxi_var = auxi_var,
      dual_var = dual_var,
      rho = rho,
      rho_adjust = rho_adjust,
      tolerance = tolerance,
      verbose = verbose
    ),
    class = "admm_initializer"
  )
}

#' ADMM solver
#'
#' Compute the Effective Reproduction Number Rt of infectious diseases using the Alternative Direction Method of Multipliers method
#'
#' @details The Effective Reproduction Number Rt of an infectious disease is estimated by the smoothed ratio between the
#' observed daily infection counts and the total infectiousness each day of the past infected. The total infectiousness at day \eqn{t} is
#' calculated by \eqn{TI_{t}\sum_{a=1}^t I_{t-s}w_s}, i.e. the expected number of new cases at day \eqn{t}
#'
#' @param current_counts the current daily infection counts
#' @param degree degree of the piecewise polynomial curve to be fitted,
#' e.g., degree = 0 corresponds to a piecewise constant curve
#' @param lambda a parameter to balance the data fidelity and graphical
#' smoothness of fitted curves; a greater lambda results in a smoother curve
#' @param maxiter maximal number of iteration
#' @param init a list of model initialization of class `admm_initializer`
#' @param dist_gamma shape and scale parameter of the discretized Gamma distribution,
#' representing the serial interval distribution, which indicate how infectious someone is if they are infected a given days ago
#' @param x the observation position. Default to 1:n, indicating case counts are evenly spaced and collected.
#' @param nsol number of lambdas to generate, if lambda is not pre-determined
#' @param lambdamin If lambda is not pre-determined, the program will generate a sequence of lambda, with lambdamin being the smallest lambda value
#' @param lambdamax If lambda is not pre-determined, the program will generate a sequence of lambda, with lambdamin being the largest lambda value
#' @param lambda_min_ratio If lambda is not pre-deteremined, and if lambdamin is not pre-determined, the program will generate a lambdamin
#' by lambdamax * lambda_min_ratio
#'
#' @return current_counts the current daily infection counts
#' @return weighted_past_counts the weighted sum of past infection counts
#' @return R_rate: the estimated effective reproduction rate
#' @return convr: if the model converges `convr==TRUE` or not `convr==FALSE`
#'
#' @export
#'
#' @examples
#' y <- c(rev(seq(2, 6, by = 1)), seq(2, 6, by = 1))
#' admm_solver(
#'   current_counts = y, weighted_past_counts = rep(1, 10), degree = 1,
#'   init = admm_initializer(current_counts = y,
#'   weighted_past_counts = rep(1, 10), degree = 1)
#' )
admm_solver <- function(current_counts,
                        degree,
                        dist_gamma = c(1, 1),
                        x = double(0),
                        lambda = double(0),
                        nsol = 100L,
                        lambdamin = -1,
                        lambdamax = -1,
                        lambda_min_ratio = 1e-4,
                        maxiter = 1e4,
                        init = NULL) {
  # create weighted past cases
  weighted_past_counts <- delay_calculator(current_counts, dist_gamma)
  if (is.null(init)) {
    init <- admm_initializer(current_counts, degree, weighted_past_counts)
  }
  if (!inherits(init, "admm_initializer")) {
    cli::cli_abort("`init` must be created with `admm_initializer()`.")
  }
  if (is.null(init$primal_var)) {
    init <- admm_initializer(
      current_counts, init$degree, weighted_past_counts,
      auxi_var = init$auxi_var, dual_var = init$dual_var)
  }
  # validate maxiter
  maxiter <- as.integer(maxiter)

  n <- length(current_counts)


  # (1) check that counts are non-negative, integer
  if (any(counts < 0)) cli::cli_abort("counts must be non-negative")
  if (!all(rlang::is_intergerish(current_counts))) cli::cli_abort("counts must be integers")


  # (2) checks on lambda, lambdamin, lambdamax
  lambda_size = length(lambda)

  if (lambda_size > 0) {

    if (nsol_int != lambda_size || nsol%%1 != nsol) cli::cli_abort("nsol must be integer, and nsol must be equal to the size of lambda when the size of lambda is greater than 0")

  } else {

    if (lambda_min_ratio < 0 || lambda_min_ratio > 1) cli::cli_abort("lambda_min_ratio must be in [0,1]")
    if (lambdamax < 0) cli::cli_abort("If lambda is not determined, lambdamax must be specified and larger than 0")
    if (lambdamin > lambdamax) cli::cli_abort("If lambda is not determined, lambdamin must be smaller than lambdamax")
  }


  # (3) check that x is a double vector of length 0 or n
  if (!is.double(x)) cli::cli_abort("x must be a double vector")
  if (!(length(x) == n | length(x) == 0)) cli::cli_abort("x must be of size either 0 or n")

  # Based on the create_lambda helper function in utils.cpp
  # If lambda has length > 0:
  #   Doesn't matter what lambdamin, lambdamax, or lambda_min_ratio are, lambdamin, lambdamax are filled based on lambda in utils.cpp
  # If lambda has length == 0:
  #   1) lambdamin > 0 && lambdamax > 0 && lambdamin < lamdamax, ignore lambda_min_ratio
  #   2) lambdamax > 0 && 0 <= lambda_min_ratio <= 1, ignore lambdamin


  # (1) check that counts are non-negative, integer
  # (2) checks on lambda, lambdamin, lambdamax (don't need to adjust)
  #   * If lambda has length > 0, lambdamin and max should be negative
  #   * Need 0 <= lambda_min_ratio <= 1
  #   * need nsol > 0, integer, equal to length(lambda) when it has positive length
  #   * lambdamin < lambdamax or both negative
  # (3) check that x is a double vector of length 0 or n


  mod <- rtestim_path(
    current_counts,
    x,
    weighted_past_counts,
    init$degree,
    lambda = lambda,
    lambdamax = lambdamax,
    lambdamin = lambdamin,
    nsol = nsol,
    rho_adjust = init$rho_adjust,
    rho = init$rho,
    maxiter = maxiter,
    tolerance = init$tolerance,
    lambda_min_ratio = lambda_min_ratio,
    verbose = init$verbose)

  structure(
    list(
      current_counts = current_counts,
      x = x %||% 1:n,
      weighted_past_counts = weighted_past_counts,
      Rt = mod$Rt,
      lambda = mod$lambda,
      degree = mod$degree,
      niter = mod$niter
    ),
    class = "admm_rr"
  )
}
