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
#' We should rename this function.
#'
#' @param current_counts the current daily infection counts
#' @param weighted_past_counts the weighted sum of past infection counts with
#' corresponding serial interval functions (or its Gamma approximation) as
#' weights
#' @param degree degree of the piecewise polynomial curve to be fitted,
#' e.g., degree = 0 corresponds to a piecewise constant curve
#' @param lambda a parameter to balance the data fidelity and graphical
#' smoothness of fitted curves; a greater lambda results in a smoother curve
#' @param mu a parameter used in the algorithm; use mu = NULL to compute it
#' using a default method (`get_mu`)
#' @param tol tolerance of convergence of primal & dual residuals
#' @param maxiter maximal number of iteration
#' @param init a list of model initialization of class `admm_initializer`
#'
#' @return current_counts the current daily infection counts
#' @return weighted_past_counts the weighted sum of past infection counts
#' @return R_rate: the estimated reproduction rate
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
