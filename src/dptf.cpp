#include <RcppArmadillo.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include "tf_dp.h"
#include "dptf.h"
// [[Rcpp::depends(RcppArmadillo)]]

void tf_dp(int n, double* y, double lam, double* beta);
void tf_dp_past(int n, double* y, double* past, double lam, double* beta);

/**
 * Filter mean (i.e., exponential of natural parameter) trend for 0th-order
 * unweighted Poisson trend filtering
 * @param y signal length
 * @param lam hyperparameter
 * @return filtered mean (exponential of natural parameter) trend
 */
// [[Rcpp::export]]
arma::vec dptf(arma::vec y, double lam) {
  int n = y.size();
  arma::vec beta(n);
  tf_dp(n, y.begin(), lam, beta.begin());  // get pointers to the vectors
  return beta;
}

/**
 * Filter mean (i.e., exponential of natural parameter) trend for 0th-order
 * weighted Poisson trend filtering
 * @param y signal length
 * @param lam hyperparameter
 * @param w signal weights
 * @return filtered mean (exponential of natural parameter) trend
 */
// [[Rcpp::export]]
arma::vec dptf_past(arma::vec y, double lam, arma::vec w) {
  int n = y.size();
  arma::vec beta(n);
  tf_dp_past(n, y.begin(), w.begin(), lam, beta.begin());
  return beta;
}
