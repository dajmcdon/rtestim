#' Estimate Rt using smoothness-penalized Poisson likelihood
#'
#' @description
#' The Effective Reproduction Number \eqn{R_t} of an infectious
#' disease can be estimated by solving the smoothness penalized Poisson
#' regression (trend filtering) of the form:
#'
#' \deqn{\hat{\theta} = \argmin_{\theta} \frac{1}{n} \sum_{i=1}^n (w_i e^{\theta_i} -
#'   y_i\theta_i) + \lambda\Vert D^{(k+1)}\theta\Vert_1, }
#'
#' where \eqn{R_t = e^{\theta}}, \eqn{y_i} is the observed case count at day
#' \eqn{i}, \eqn{w_i} is the weighted past counts at day \eqn{i}, \eqn{\lambda}
#' is the smoothness penalty, and \eqn{D^{(k+1)}} is the \eqn{(k+1)}-th order
#' difference matrix.
#'
#' @param observed_counts vector of the observed daily infection counts
#' @param korder Integer. Degree of the piecewise polynomial curve to be
#'   estimated. For example, `korder = 0` corresponds to a piecewise constant
#'   curve.
#' @param lambda Vector. A user supplied sequence of tuning parameters which
#'   determines the balance between data fidelity and
#'   smoothness of the estimated Rt; larger `lambda` results in a smoother
#'   estimate. The default, `NULL`
#'   results in an automatic computation based on `nlambda`, the largest value
#'   of `lambda` that would result in a maximally smooth estimate, and `lambda_min_ratio`.
#'   Supplying a value of `lambda` overrides
#'   this behaviour. It is likely better to supply a
#'   decreasing sequence of `lambda` values than a single (small) value. If
#'   supplied, the user-defined `lambda` sequence is automatically sorted in
#'   decreasing order.
#' @param maxiter Integer. Maximum number of iterations for the estimation
#'   algorithm.
#' @param init a list of internal configuration parameters of class
#'   `rt_admm_configuration`.
#' @param dist_gamma Vector of length 2. These are the shape and scale for the
#'   assumed serial interval distribution. Roughly, this distribution describes
#'   the probability of an infectious individual infecting someone else after
#'   some period of time after having become infectious.
#'   As in most literature, we assume that this interval follows a gamma
#'   distribution with some shape and scale.
#' @param x a vector of positions at which the counts have been observed. In an
#'   ideal case, we would observe data at regular intervals (e.g. daily or
#'   weekly) but this may not always be the case. May be numeric or Date.
#' @param nsol Integer. The number of tuning parameters `lambda` at which to
#'   compute Rt.
#' @param delay_distn in the case of a non-gamma delay distribution,
#'   a vector or matrix (or `Matrix::Matrix()`) of delay probabilities may be
#'   passed here. For a vector, these will be coerced
#'   to sum to 1, and padded with 0 in the right tail if necessary. If a
#'   time-varying delay matrix, it must be lower-triangular. Each row will be
#'   silently coerced to sum to 1. See also `vignette("delay-distributions")`.
#' @param delay_distn_periodicity Controls the relationship between the spacing
#'   of the computed delay distribution and the spacing of `x`. In the default
#'   case, `x` would be regular on the sequence `1:length(observed_cases)`,
#'   and this would
#'   be 1. But if `x` is a `Date` object or spaced irregularly, the relationship
#'   becomes more complicated. For example, weekly data when `x` is a date in
#'   the form `YYYY-MM-DD` requires specifying `delay_distn_periodicity = "1 week"`.
#'   Or if `observed_cases` were reported on Monday, Wednesday, and Friday,
#'   then `delay_distn_periodicity = "1 day"` would be most appropriate.
#' @param lambdamin Optional value for the smallest `lambda` to use. This should
#'   be greater than zero.
#' @param lambdamax Optional value for the largest `lambda` to use.
#' @param lambda_min_ratio If neither `lambda` nor `lambdamin` is specified, the
#'   program will generate a lambdamin by lambdamax * lambda_min_ratio.
#'   A multiplicative factor for the minimal lambda in the
#'   `lambda` sequence, where `lambdamin = lambda_min_ratio * lambdamax`.
#'   A very small value will lead to the solution `Rt = log(observed_counts)`.
#'   This argument has no effect if there is a user-defined `lambda` sequence.
#'
#' @return An object with S3 class `poisson_rt`. Among the list components:
#' * `observed_counts` the observed daily infection counts.
#' * `x` a vector of positions at which the counts have been observed.
#' * `weighted_past_counts` the weighted sum of past infection counts.
#' * `Rt` the estimated effective reproduction rate. This is a matrix with
#'     each column corresponding to one value of `lambda`.
#' * `lambda` the values of `lambda` actually used in the algorithm.
#' * `korder` degree of the estimated piecewise polynomial curve.
#' * `dof` degrees of freedom of the estimated trend filtering problem.
#' * `niter` the required number of iterations for each value of `lambda`.
#' * `convergence` if number of iterations for each value of `lambda` is less
#'     than the maximum number of iterations for the estimation algorithm.
#'
#' @export
#'
#' @examples
#' y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
#' out <- estimate_rt(y)
#' plot(out)
#'
#' out0 <- estimate_rt(y, korder = 0L, nsol = 10)
#' plot(out0)
estimate_rt <- function(
    observed_counts,
    korder = 3L,
    dist_gamma = c(2.5, 2.5),
    x = 1:n,
    lambda = NULL,
    nsol = 100L,
    delay_distn = NULL,
    delay_distn_periodicity = NULL,
    lambdamin = NULL,
    lambdamax = NULL,
    lambda_min_ratio = 1e-4,
    maxiter = 1e5,
    init = configure_rt_admm()) {

  assert_int(nsol, lower = 1)
  assert_int(maxiter, lower = 1)
  assert_number(lambda_min_ratio, lower = 0, upper = 1)
  assert_number(lambdamin, lower = 0, null.ok = TRUE)
  assert_number(lambdamax, lower = 0, null.ok = TRUE)
  assert_int(delay_distn_periodicity, lower = 1, null.ok = TRUE)
  assert_numeric(dist_gamma, lower = .Machine$double.eps, finite = TRUE, len = 2)
  assert_class(init, "rt_admm_configuration")
  assert_numeric(observed_counts, lower = 0)

  ymiss <- is.na(observed_counts)
  if (any(ymiss)) {
    cli_warn("Missing values in `observed_counts` will be ignored.")
  }
  observed_counts <- observed_counts[!ymiss]
  n <- length(observed_counts)
  x <- x[!ymiss]

  assert_int(korder, lower = 0, upper = n - 2L)

  xin <- x
  if (inherits(x, "Date")) x <- as.numeric(x)
  assert_numeric(x, len = n, any.missing = FALSE)
  if (is.unsorted(x, strictly = TRUE)) {
    cli_abort("`x` must be sorted in increasing order without duplicates.")
  }

  weighted_past_counts <- delay_calculator(
    observed_counts, x, dist_gamma, delay_distn, delay_distn_periodicity
  )



  # checks on lambda, lambdamin, lambdamax
  if (is.null(lambda)) lambda <- double(nsol) # prep for create_lambda
  if (is.null(lambdamin)) lambdamin <- -1.0
  if (is.null(lambdamax)) lambdamax <- -1.0
  if (length(lambda) == 0) {
    msg <- "If lambda is not specified,"
    if (lambda_min_ratio >= 1) {
      cli_abort("{msg} lambda_min_ratio must be in (0,1).")
    }
    if (lambdamin > 0 && lambdamax > 0 && lambdamin >= lambdamax) {
      cli_abort("{msg} lambdamin must be < lambdamax.")
    }
    lambda <- double(nsol)
  }
  if (length(lambda) != nsol) nsol <- length(lambda)
  lambda <- sort(lambda, decreasing = TRUE)

  mod <- rtestim_path(
    observed_counts,
    x,
    weighted_past_counts,
    korder,
    lambda = lambda,
    lambdamax = lambdamax,
    lambdamin = lambdamin,
    nsol = nsol,
    rho = init$rho,
    maxiter = maxiter,
    maxiter_newton = init$maxiter_newton,
    maxiter_line = init$maxiter_line,
    tolerance = init$tolerance,
    lambda_min_ratio = lambda_min_ratio,
    ls_alpha = init$alpha,
    ls_gamma = init$gamma,
    verbose = init$verbose
  )

  structure(
    enlist(
      observed_counts,
      x = xin,
      weighted_past_counts,
      Rt = drop(mod$Rt),
      lambda = drop(mod$lambda),
      korder = mod$korder,
      dof = drop(mod$nknots) + mod$korder + 1,
      niter = drop(mod$niter),
      convergence = (mod$niter < maxiter),
      call = match.call(),
      alp = drop(mod$alp),
      delay_distn_periodicity,
      tolerance = init$tolerance
    ),
    class = "poisson_rt"
  )
}


#' Rt estimation algorithm configuration
#'
#' @param rho Double. An ADMM parameter; coefficient of augmented term in the
#' Lagrangian function.
#' @param alpha Double. A parameter adjusting upper bound in line search algorithm
#'   in `prox_newton` algorithm.
#' @param gamma Double. A parameter adjusting step size in line search algorithm
#'   in `prox_newton` algorithm.
#' @param tolerance Double. Tolerance of ADMM convergence.
#' @param maxiter_newton Integer. Maximum number of iterations for the outer
#'   Newton iteration.
#' @param maxiter_line Integer. Maximum number of iterations for the linesearch
#'   algorithm in the proximal Newton method.
#' @param verbose Integer.
#' @param ... space for future extensions
#'
#' @return a list of model parameters with class `rt_admm_configuration`
#'
#' @export
configure_rt_admm <- function(
    rho = -1,
    alpha = 0.5,
    gamma = 0.9,
    tolerance = 1e-4,
    maxiter_newton = 50L,
    maxiter_line = 20L,
    verbose = 0,
    ...) {
  rlang::check_dots_empty()
  assert_number(alpha, lower = 0, upper = 1)
  assert_number(gamma, lower = 0, upper = 1)
  assert_int(maxiter_newton, lower = 1L)
  assert_int(maxiter_line, lower = 1L)
  assert_int(verbose, lower = 0L)
  assert_number(tolerance, lower = 0)
  assert_number(rho)
  if (abs(rho + 1) > sqrt(.Machine$double.eps)) assert_number(rho, lower = 0)

  structure(
    enlist(
      rho,
      alpha,
      gamma,
      tolerance,
      maxiter_newton,
      maxiter_line,
      verbose
    ),
    class = "rt_admm_configuration"
  )
}
