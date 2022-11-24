#' ADMM initialization
#'
#' @param y observed daily infection count
#' @param x scales of reproduction rate
#' @param k order of divided difference matrix D
#' @param z auxiliary variable
#' @param u dual variable
#'
#' @return initialization of auxiliary and dual variables
#' @export
admm_initializer <- function(y, x, k, z = NULL, u = NULL) {
    n <- length(y)
    k <- as.integer(k)
    if (k < 0) {
        abort("`k` is at least 0. ")
    } else if (k == 0) {
        Dk <- generate_I(n)
    } else {
        Dk <- generate_Dk(n, k - 1)
    }

    # lambda_max

    if (!is.null(z)) {
        if (length(z) != n - k) {
            abort("`z` must have length {n - k}.")
        }
    } else {
        z <- double(n - k)
    }
    if (!is.null(u)) {
        if (length(u) != n - k) {
            abort("`u` must have length {n - k}.")
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
#' @param k order of divided difference matrix D
#' @param lambda Lagrangian multiplier
#' @param mu a hyper parameter depending on order, problem size, and
#' @param tol tolerance of convergence of primal & dual residuals
#' @param init ADMM initializer
#'
#' @return model outputs; model convergence
#' @export
#'
#' @examples
#' admm_solver(
#'     y = c(rev(seq(0.2, 0.6, by = 0.1)), seq(0.2, 0.6, by = 0.1)),
#'     x = rep(1, 10), k = 1, mu = 0.078,
#'     init = admm_initializer(
#'         y = c(rev(seq(0.2, 0.6, by = 0.1)), seq(0.2, 0.6, by = 0.1)),
#'         x = rep(1, 10), k = 1
#'     )
#' )
admm_solver <- function(y, x, k, mu, lambda = 0.01, tol = 1e-3, maxiter = 1e3,
                        init = admm_initializer(y, x, k)) {
    if (!inherits(init, "admm_initializer")) {
        abort("`init` should be created with `admm_initializer()`.")
    }
    maxiter <- as.integer(maxiter)
    n <- length(y)
    rho <- lambda

    mod <- admm(
        M = maxiter, y = y, x = x, n = n, theta = double(n), z = init$z,
        u = init$u, lambda = lambda, rho = rho, mu = mu, D = init$Dk,
        tol = tol
    )
    convr <- (mod$iter_num < maxiter)
    # nl <- (lambda >= init$lambda_max)

    structure(
        list(
            primal_var = mod$theta, aux_var = mod$z, dual_var = mod$u,
            primal_res = mod$prim_res, dual_res = mod$dual_res,
            iter_num = mod$iter_num, convergence = convr
        ),
        class = "admm_rr"
    )
}
