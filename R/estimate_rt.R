#' Estimate Rt using smoothness-penalized Poisson likelihood
#'
#' @description
#' The Effective Reproduction Number \eqn{R_t} of an infectious
#' disease can be estimated by solving the smoothness penalized Poisson
#' regression of the form:
#'
#' \eqn{R_t = argmin_{\theta} (\frac{1}{n} \sum_{i=1}^n e^{\theta_i} -
#'   y_i\theta_i) + \lambda||D^{(k+1)}\theta||_1}
#'
#' where \eqn{y_i} is the observed case count at day \eqn{i},
#' \eqn{\theta_i = \sum_{a=1}y_{a}w_{t-a}} is the weighted past counts
#' at day \eqn{i}, \eqn{\lambda} is the smoothness penalty, and \eqn{D^{(k+1)}}
#' is the \eqn{(k+1)}-th order difference matrix.
#'
#' @param observed_counts vector of the observed daily infection counts
#' @param degree Integer. Degree of the piecewise polynomial curve to be
#'   estimated. Ror example, `degree = 0` corresponds to a piecewise constant
#'   curve.
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
#'   weekly) but this may not always be the case.
#' @param nsol Integer. The number of tuning parameters `lambda` at which to
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
#' @param algo the algorithm to be used in computation. `linear_admm`:
#'   linearized ADMM; `irls_admm`: iteratively reweighted least squares with
#'   standard ADMM.
#'
#' @return An object with S3 class `"poison_rt"`. Among the list components:
#' * `observed_counts` the observed daily infection counts
#' * `weighted_past_counts` the weighted sum of past infection counts
#' * `R` the estimated effective reproduction rate. This is a matrix with
#'     each column corresponding to one value of `lambda`.
#' * `lambda` the value of `lambda` actually used in the algorithm.
#' * `niter` the required number of iterations for each value of `lambda`
#' * `convr` if the model converges `convr==TRUE` or not `convr==FALSE`
#'
#' @export
#'
#' @examples
#' y <- c(rev(seq(2, 6, by = 1)), seq(2, 6, by = 1))
#' estimate_rt(
#'   observed_counts = y, degree = 2, lambda = .1,
#'   algo = "linear_admm",
#'   init = rt_admm_configuration(y, degree = 1)
#' )
estimate_rt <- function(observed_counts,
                        degree = 3L,
                        dist_gamma = c(2.5, 2.5),
                        x = NULL,
                        lambda = NULL,
                        nsol = 100L,
                        lambdamin = NULL,
                        lambdamax = NULL,
                        lambda_min_ratio = 1e-4,
                        algo = c("linear_admm", "irls_admm"),
                        maxiter = 1e4,
                        init = NULL) {

  arg_is_nonneg_int(degree)
  arg_is_pos_int(nsol, maxiter)
  arg_is_scalar(degree, nsol, lambda_min_ratio)
  arg_is_scalar(lambdamin, lambdamax, allow_null = TRUE)
  arg_is_positive(lambdamin, lambdamax, allow_null = TRUE)
  arg_is_positive(lambda_min_ratio, dist_gamma)
  arg_is_length(2, dist_gamma)
  algo <- match.arg(algo)

  # create weighted past cases, do setup
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

  n <- length(observed_counts)

  # (1) check that counts are non-negative, integer
  if (any(observed_counts < 0))
    cli::cli_abort("`observed_counts` must be non-negative")

  # (2) checks on lambda, lambdamin, lambdamax
  if (is.null(lambda)) lambda <- double(0) # prep for create_lambda
  if (is.null(lambdamin)) lambdamin <- -1.0
  if (is.null(lambdamax)) lambdamax <- -1.0
  if (length(lambda) == 0) {
    msg <- "If lambda is not specified,"
    if (lambda_min_ratio >= 1)
      cli::cli_abort("{msg} lambda_min_ratio must be in (0,1)")
    if (lambdamin > 0 && lambdamax > 0 && lambdamin >= lambdamax)
      cli::cli_abort("{msg} lambdamin must be < lambdamax.")
  }

  # (3) check that x is a double vector of length 0 or n
  x <- x %||% 1:n
  if (!is.numeric(x)) cli::cli_abort("x must be a numeric vector")
  x <- as.double(x)
  if (!(length(x) == n | length(x) == 0))
    cli::cli_abort("x must be length 0 or n")

  # (4) check algorithm
  algo <- match(algo, c("linear_admm", "irls_admm"))
  algo <- as.integer(algo)

  mod <- rtestim_path(
    algo,
    observed_counts,
    x,
    weighted_past_counts,
    init$degree,
    lambda = lambda,
    lambdamax = lambdamax,
    lambdamin = lambdamin,
    nsol = nsol,
    rho = init$rho,
    maxiter = maxiter,
    tolerance = init$tolerance,
    lambda_min_ratio = lambda_min_ratio,
    ls_alpha = init$alpha,
    ls_gamma = init$gamma,
    verbose = init$verbose)

  structure(
    list(
      observed_counts = observed_counts,
      x = x,
      weighted_past_counts = weighted_past_counts,
      Rt = mod$Rt,
      lambda = mod$lambda,
      degree = mod$degree,
      maxiter = maxiter,
      niter = mod$niter
    ),
    class = "poisson_rt"
  )
}


#' Rt estimation algorithm configuration
#'
#' We should convert this into an S3 method so that we can pass in a
#' vector of counts or an admm_initializer object (and overwrite)
#'
#' @inheritParams estimate_rt
#' @param weighted_past_counts the weighted sum of past infections counts with
#'   corresponding serial interval functions (or its Gamma approximation) as
#'   weights
#' @param primal_var initial values of log(Rt)
#' @param auxi_var auxiliary variable in the ADMM algorithm
#' @param dual_var dual variable in the ADMM algorithm
#' @param alpha Double. A parameter adjusting upper bound in line search algorithm
#'   in `irls_admm` algorithm.
#' @param gamma Double. A parameter adjusting step size in line search algorithm
#'   in `irls_admm` algorithm.
#'
#' @return a list of model parameters with class `rt_admm_configuration`
#'
#' @export
rt_admm_configuration <- function(observed_counts,
                                  degree,
                                  weighted_past_counts = NULL,
                                  primal_var = NULL,
                                  auxi_var = NULL,
                                  dual_var = NULL,
                                  rho = -1,
                                  rho_adjust = -1,
                                  alpha = 0.5,
                                  gamma = 0.9,
                                  tolerance = 1e-4,
                                  verbose = 0) {
  n <- length(observed_counts)
  arg_is_scalar(degree, rho, rho_adjust, alpha, gamma, tolerance, verbose)
  arg_is_positive(alpha, gamma, tolerance)
  arg_is_numeric(rho, rho_adjust, tolerance, verbose)
  arg_is_nonneg_int(degree)
  if (alpha >= 1) cli::cli_abort("alpha must be in (0, 1).")
  if (gamma > 1) cli::cli_abort("gamma must be in (0, 1].")

  if (is.null(primal_var)) {
    if (!is.null(weighted_past_counts))
      primal_var <- log(observed_counts / (n * weighted_past_counts))
  } else {
    arg_is_length(n, primal_var)
  }
  if (is.null(auxi_var)) auxi_var <- double(n - degree)
  else arg_is_length(n - degree, auxi_var)

  if (is.null(dual_var)) dual_var <- double(n - degree)
  else arg_is_length(n - degree, dual_var)

  structure(
    list(
      degree = degree,
      primal_var = primal_var,
      auxi_var = auxi_var,
      dual_var = dual_var,
      rho = rho,
      rho_adjust = rho_adjust,
      tolerance = tolerance,
      alpha = alpha,
      gamma = gamma,
      verbose = as.integer(verbose)
    ),
    class = "rt_admm_configuration"
  )
}
