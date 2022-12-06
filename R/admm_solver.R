#' ADMM initialization
#'
#' @param current_counts the current daily infection counts
#' @param weighted_past_counts the weighted sum of past infection counts with
#' corresponding serial interval functions (or its Gamma approximation) as
#' weights
#' @param degree degree of the piecewise polynomial curves to be fitted,
#' e.g., degree = 0 corresponds to piecewise constant curves
#' @param auxi_var auxiliary variable in the ADMM algorithm
#' @param dual_var dual variable in the ADMM algorithm
#'
#' @return a list of model initialization with class `admm_initializer`
#'
#' @export
admm_initializer <- function(current_counts, weighted_past_counts, degree,
                             auxi_var = NULL, dual_var = NULL) {
  n <- length(current_counts)
  degree <- as.integer(degree)
  if (degree < 0) {
    cli::cli_abort("`degree` must be non-negative. ")
  } else {
    D <- generate_D(n, degree)
  }

  # lambda_max?

  if (!is.null(auxi_var)) {
    if (length(auxi_var) != n - degree) {
      cli::cli_abort("`auxi_var` must have length {n - degree}.")
    }
  } else {
    auxi_var <- double(n - degree)
  }
  if (!is.null(dual_var)) {
    if (length(dual_var) != n - degree) {
      cli::cli_abort("`dual_var` must have length {n - degree}.")
    }
  } else {
    dual_var <- double(n - degree)
  }

  structure(list(auxi_var = auxi_var, dual_var = dual_var,
                 D = D), # lambda_max = lambda_max
    class = "admm_initializer"
  )
}

#' ADMM solver
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
admm_solver <- function(current_counts, weighted_past_counts, degree,
                        lambda = 0.01, mu = NULL, tol = 1e-3,  maxiter = 1e3L,
                        init = admm_initializer(current_counts,
                                                weighted_past_counts,
                                                degree)) {
  if (!inherits(init, "admm_initializer")) {
    cli::cli_abort("`init` must be created with `admm_initializer()`.")
  }
  maxiter <- as.integer(maxiter)
  n <- length(current_counts)
  rho <- lambda

  if (is.null(mu)) mu <- get_mu(n, degree, lambda)

  mod <- admm(
    M = maxiter, y = current_counts, x = weighted_past_counts, n = n,
    theta = double(n), z = init$auxi_var, u = init$dual_var, lambda = lambda,
    rho = rho, mu = mu, D = init$D, tol = tol
  )
  convr <- (mod$iter_num < maxiter)

  structure(list(current_counts = current_counts,
                 weighted_past_counts = weighted_past_counts,
                 R_rate = as.vector(exp(mod$theta)),
                 convr = convr
                 ),
            class = "admm_rr")
}
