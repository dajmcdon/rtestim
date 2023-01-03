#' ADMM initialization
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
#' @return a list of model initialization with class `admm_initializer`
#'
#' @export
admm_initializer <- function(current_counts, degree,
                             primal_var = NULL,
                             auxi_var = NULL, dual_var = NULL) {
  n <- length(current_counts)
  degree <- as.integer(degree)
  if (degree < 0) cli::cli_abort("`degree` must be non-negative. ")
  weighted_past_cases <- delay_calculator(current_counts, dist_gamma)

  if (is.null(primal_var)) {
    # should we divide by n?
    # what do we do when current_counts == 0
    primal_var <- log(current_counts / (n * weighted_past_counts))

  } else {
    if (length(primal_var) != n) {
      cli::cli_abort("`primal_var` must have length {n - degree}.")
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
    list(primal_var = primal_var,
         auxi_var = auxi_var,
         dual_var = dual_var),
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
                        lambda = NULL,
                        lambdamin = -1,
                        lambdamax = -1,
                        mu = NULL,
                        tol = 1e-3,
                        maxiter = 1e4,
                        init = NULL) {
  # create weighted past cases

  if (is.null(init)) {
    init <- admm_initializer(current_counts, degree, weighted_past_cases)

  }
  if (!inherits(init, "admm_initializer")) {
    cli::cli_abort("`init` must be created with `admm_initializer()`.")
  }
  maxiter <- as.integer(maxiter)
  n <- length(current_counts)
  rho <- lambda

  # (1) check that counts are non-negative, integer
  # (2) checks on lambda, lambdamin, lambdamax (don't need to adjust)
  #   * If lambda is non-null, lambdamin and max should be negative
  #   * Need 0 <= lambda_min_ratio <= 1
  #   * need nsol > 0, integer
  #   * lambdamin < lambdamax or both negative
  # (3) check that x is a double vector of length 0 or n

  if (is.null(mu)) mu <- get_mu(n, degree, lambda)


  mod <- rtestim_path(
    current_counts,
    x,
    weighted_past_cases,
    korder,
    lambda = lambda,
    lambdamax = lambdamax,
    lambdamin = lambdamin,
    nsol = nsol,
    rho_adjust = rho_adjust,
    rho = -1,
    maxiter = 1e5,
    tolerance = 1e-3,
    lambda_min_ratio = 1e-4,
    verbose = 0)

  mod <- rtestim_path(
    M = maxiter,
    y = current_counts,
    x = weighted_past_counts,
    n = n,
    theta = double(n),
    z = init$auxi_var,
    u = init$dual_var,
    lambda = lambda,
    rho = rho,
    mu = mu,
    D = init$D,
    tol = tol
  )

  structure(list(current_counts = current_counts,
                 weighted_past_counts = weighted_past_counts,
                 R_rate = as.vector(exp(mod$theta)),
                 ),
            class = "admm_rr")
}
