% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-cv.R
\name{predict.cv_poisson_rt}
\alias{predict.cv_poisson_rt}
\title{Predict observed data using estimated Rt}
\usage{
\method{predict}{cv_poisson_rt}(object, which_lambda = c("lambda.min", "lambda.1se"), ...)
}
\arguments{
\item{object}{result of cross validation of type \code{cv_poisson_rt}}

\item{which_lambda}{Select which lambdas from the object to use. If not
provided, all Rt's are returned. Note that new lambdas not originally
used in the estimation procedure may be provided, but the results will be
calculated by linearly interpolating the estimated Rt's.

The strings \code{lambda.min} or \code{lambda.1se} are allowed to choose either
the lambda that minimizes the cross validation score or the largest lambda
whose corresponding cross validation score is within 1 standard error of
the minimal cross validation score.}

\item{...}{not used.}
}
\value{
A vector or matrix of predicted case counts.
}
\description{
Given an object of class \code{poisson_rt} produced with \code{\link[=estimate_rt]{estimate_rt()}},
calculate predicted observed cases for the estimated Rt values.
Note: This function is not intended for "new x" or to produce forecasts, but
rather to examine how Rt relates to observables.
}
\examples{
y <- c(1, rpois(100, dnorm(1:100, 50, 15) * 500 + 1))
cv <- cv_estimate_rt(y, korder = 3, nfold = 3, nsol = 30)
p <- predict(cv)
p <- predict(cv, which_lambda = cv$lambda[1])
p <- predict(cv, which_lambda = "lambda.1se")
p <- predict(cv, which_lambda = NULL)
plot(y)
matlines(p, lty = 2)
}
