#' ADMM initialization
#'
#' @param y observed daily infection count
#' @param x scales of reproduction rate
#' @param k order of divided difference matrix D
#' @param z auxiliary variable
#' @param u dual variable
#'
#' @return a list of model initialization
#' @export
admm_initializer <- function(y, x, k, z = NULL, u = NULL) {
  n <- length(y)
  k <- as.integer(k)
  if (k < 0) {
    rlang::abort("`k` is at least 0. ")
  } else if (k == 0) {
    Dk <- generate_I(n)
  } else {
    Dk <- generate_Dk(n, k - 1)
  }

  # lambda_max?

  if (!is.null(z)) {
    if (length(z) != n - k) {
      rlang::abort("`z` must have length {n - k}.")
    }
  } else {
    z <- double(n - k)
  }
  if (!is.null(u)) {
    if (length(u) != n - k) {
      rlang::abort("`u` must have length {n - k}.")
    }
  } else {
    u <- double(n - k)
  }

  structure(list(z = z, u = u, Dk = Dk), # lambda_max = lambda_max
    class = "admm_initializer"
  )
}

#' ADMM solver
#'
#' @param y observed daily infection count
#' @param x scales of reproduction rate
#' @param k order of divided difference matrix
#' @param lambda Lagrangian multiplier
#' @param mu a hyper parameter depending on order, problem size, and lambda
#' @param tol tolerance of convergence of primal & dual residuals
#' @param maxiter maximal number of iteration
#' @param init ADMM initializer
#' @return y: observed daily infection count
#' @return x: scales of reproduction rate
#' @return R_rate: reproduction rate
#' @return convr: if the model converges (convr==TRUE) or not (convr==FALSE)
#' @export
#'
#' @examples
#' y <- c(rev(seq(2, 6, by = 1)), seq(2, 6, by = 1))
#' admm_solver(
#'   y = y, x = rep(1, 10), k = 1,
#'   init = admm_initializer(y = y, x = rep(1, 10), k = 1)
#' )
admm_solver <- function(y, x, k, lambda = 0.01, mu = NULL, tol = 1e-3,
                        maxiter = 1e3L,
                        init = admm_initializer(y, x, k)) {
  if (!inherits(init, "admm_initializer")) {
    rlang::abort("`init` should be created with `admm_initializer()`.")
  }
  maxiter <- as.integer(maxiter)
  n <- length(y)
  rho <- lambda

  if (is.null(mu)) mu <- get_mu(n, k, lambda)

  mod <- admm(
    M = maxiter, y = y, x = x, n = n, theta = double(n), z = init$z,
    u = init$u, lambda = lambda, rho = rho, mu = mu, D = init$Dk,
    tol = tol
  )
  convr <- (mod$iter_num < maxiter)
  # nl <- (lambda >= init$lambda_max)

  structure(list(y = y, x = x, R_rate = as.vector(exp(mod$theta)),
                 convr = convr
                 ), class = "admm_rr")
}
