% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_rt.R
\name{estimate_rt}
\alias{estimate_rt}
\title{Estimate Rt using smoothness-penalized Poisson likelihood}
\usage{
estimate_rt(
  observed_counts,
  korder = 3L,
  dist_gamma = c(2.5, 2.5),
  x = 1:n,
  lambda = NULL,
  nsol = 50L,
  delay_distn = NULL,
  delay_distn_periodicity = NULL,
  lambdamin = NULL,
  lambdamax = NULL,
  lambda_min_ratio = 1e-04,
  maxiter = 1e+05,
  init = configure_rt_admm()
)
}
\arguments{
\item{observed_counts}{vector of the observed daily infection counts}

\item{korder}{Integer. Degree of the piecewise polynomial curve to be
estimated. For example, \code{korder = 0} corresponds to a piecewise constant
curve.}

\item{dist_gamma}{Vector of length 2. These are the shape and scale for the
assumed serial interval distribution. Roughly, this distribution describes
the probability of an infectious individual infecting someone else after
some period of time after having become infectious.
As in most literature, we assume that this interval follows a gamma
distribution with some shape and scale.}

\item{x}{a vector of positions at which the counts have been observed. In an
ideal case, we would observe data at regular intervals (e.g. daily or
weekly) but this may not always be the case. May be numeric or Date.}

\item{lambda}{Vector. A user supplied sequence of tuning parameters which
determines the balance between data fidelity and
smoothness of the estimated Rt; larger \code{lambda} results in a smoother
estimate. The default, \code{NULL}
results in an automatic computation based on \code{nlambda}, the largest value
of \code{lambda} that would result in a maximally smooth estimate, and \code{lambda_min_ratio}.
Supplying a value of \code{lambda} overrides
this behaviour. It is likely better to supply a
decreasing sequence of \code{lambda} values than a single (small) value. If
supplied, the user-defined \code{lambda} sequence is automatically sorted in
decreasing order.}

\item{nsol}{Integer. The number of tuning parameters \code{lambda} at which to
compute Rt.}

\item{delay_distn}{in the case of a non-gamma delay distribution,
a vector or matrix (or \code{Matrix::Matrix()}) of delay probabilities may be
passed here. For a vector, these will be coerced
to sum to 1, and padded with 0 in the right tail if necessary. If a
time-varying delay matrix, it must be lower-triangular. Each row will be
silently coerced to sum to 1. See also \code{vignette("delay-distributions")}.}

\item{delay_distn_periodicity}{Controls the relationship between the spacing
of the computed delay distribution and the spacing of \code{x}. In the default
case, \code{x} would be regular on the sequence \code{1:length(observed_cases)},
and this would
be 1. But if \code{x} is a \code{Date} object or spaced irregularly, the relationship
becomes more complicated. For example, weekly data when \code{x} is a date in
the form \code{YYYY-MM-DD} requires specifying \code{delay_distn_periodicity = "1 week"}.
Or if \code{observed_cases} were reported on Monday, Wednesday, and Friday,
then \code{delay_distn_periodicity = "1 day"} would be most appropriate.}

\item{lambdamin}{Optional value for the smallest \code{lambda} to use. This should
be greater than zero.}

\item{lambdamax}{Optional value for the largest \code{lambda} to use.}

\item{lambda_min_ratio}{If neither \code{lambda} nor \code{lambdamin} is specified, the
program will generate a lambdamin by lambdamax * lambda_min_ratio.
A multiplicative factor for the minimal lambda in the
\code{lambda} sequence, where \code{lambdamin = lambda_min_ratio * lambdamax}.
A very small value will lead to the solution \code{Rt = log(observed_counts)}.
This argument has no effect if there is a user-defined \code{lambda} sequence.}

\item{maxiter}{Integer. Maximum number of iterations for the estimation
algorithm.}

\item{init}{a list of internal configuration parameters of class
\code{rt_admm_configuration}.}
}
\value{
An object with S3 class \code{poisson_rt}. Among the list components:
\itemize{
\item \code{observed_counts} the observed daily infection counts.
\item \code{x} a vector of positions at which the counts have been observed.
\item \code{weighted_past_counts} the weighted sum of past infection counts.
\item \code{Rt} the estimated effective reproduction rate. This is a matrix with
each column corresponding to one value of \code{lambda}.
\item \code{lambda} the values of \code{lambda} actually used in the algorithm.
\item \code{korder} degree of the estimated piecewise polynomial curve.
\item \code{dof} degrees of freedom of the estimated trend filtering problem.
\item \code{niter} the required number of iterations for each value of \code{lambda}.
\item \code{convergence} if number of iterations for each value of \code{lambda} is less
than the maximum number of iterations for the estimation algorithm.
}
}
\description{
The Effective Reproduction Number \eqn{R_t} of an infectious
disease can be estimated by solving the smoothness penalized Poisson
regression (trend filtering) of the form:

\deqn{\hat{\theta} = \arg\min_{\theta} \frac{1}{n} \sum_{i=1}^n (w_i e^{\theta_i} -
  y_i\theta_i) + \lambda\Vert D^{(k+1)}\theta\Vert_1, }

where \eqn{R_t = e^{\theta}}, \eqn{y_i} is the observed case count at day
\eqn{i}, \eqn{w_i} is the weighted past counts at day \eqn{i}, \eqn{\lambda}
is the smoothness penalty, and \eqn{D^{(k+1)}} is the \eqn{(k+1)}-th order
difference matrix.
}
\examples{
y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
out <- estimate_rt(y)
out
plot(out)

out0 <- estimate_rt(y, korder = 0L, nsol = 40)
out0
plot(out0)
}
