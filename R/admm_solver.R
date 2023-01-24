#' Estimate Rt using smoothness-penalized Poisson likelihood
#'
#' @description
#' The Effective Reproduction Number \eqn{R_t} of an infectious
#' disease can be estimated by solving the smoothnesss penalized Poisson
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
#' @param nlambda Integer. The number of tuning parameters `lambda` at which to
#'   compute Rt.
#' @param lambdamin Optional value for the smallest `lambda` to use. This should
#'   be greater than zero.
#' @param lambdamax Optional value for the largest `lambda` to use.
#' @param lambdamin_ratio If neither `lambda` nor `lambdamin` is specified, the
#'   program will generate a lambdamin by lambdamax * lambda_min_ratio
#'   A multiplicative factor for the minimal lambda in the
#'   `lambda` sequence, where `lambdamin = lambdamin_ratio * lambdamax`.
#'   A very small value will lead to the solution `Rt = log(observed_counts)`.
#'   This argument has no effect if there is user-defined `lambda` sequence.
#'
#' @return An object with S3 class `"poisson_rt"`. Among the list components:
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
#' admm_solver(
#'   observed_counts = y, weighted_past_counts = rep(1, 10), degree = 1,
#'   init = admm_initializer(observed_counts = y,
#'   weighted_past_counts = rep(1, 10), degree = 1)
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

  # (1) check that counts are non-negative, integer
  if (any(observed_counts < 0)) cli::cli_abort("`observed_counts` must be non-negative")
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
    if (lambdamin < 0)
      cli::cli_abort("{msg} lambdamin must be positive.")
    if (lambdamin >= lambdamax)
      cli::cli_abort("{msg} lambdamin must be < lambdamax.")
  }


  # (3) check that x is a double vector of length 0 or n
  if (!is.double(x)) x = as.double(x)
  if (!is.numeric(x)) cli::cli_abort("x must be a numeric vector")
  if (!(length(x) == n | length(x) == 0))
    cli::cli_abort("x must be length 0 or n")

  # (4) check degree smaller than data length
  if (init$degree < n) cli::cli_abort("degree must be < than length of observed_counts")


  mod <- rtestim_path(
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
    verbose = init$verbose)

  structure(
    list(
      observed_counts = observed_counts,
      x = x %||% 1:n,
      weighted_past_counts = weighted_past_counts,
      Rt = mod$Rt,
      lambda = mod$lambda,
      degree = mod$degree,
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
#' @param observed_counts vector of daily infections
#' @param weighted_past_counts the weighted sum of past infections counts with
#'   corresponding serial interval functions (or its Gamma approximation) as
#'   weights
#' @param degree degree of the piecewise polynomial curves to be fitted,
#'   e.g., degree = 0 corresponds to piecewise constant curves
#' @param primal_var initial values of log(Rt)
#' @param auxi_var auxiliary variable in the ADMM algorithm
#' @param dual_var dual variable in the ADMM algorithm
#'
#' @return a list of model parameters with class `admm_initializer`
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
                                  tolerance = 1e-4,
                                  verbose = 0) {
  n <- length(observed_counts)
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
    verbose = 0L
  }

  if (is.null(primal_var)) {
    if (!is.null(weighted_past_counts)) {
      # should we divide by n?
      # what do we do when observed_counts == 0
      primal_var <- log(observed_counts / (n * weighted_past_counts))
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
    class = "rt_admm_configuration"
  )
}
