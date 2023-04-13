#include <Rcpp.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include "tf_dp.h"
#include "dptf.h"

/**
 * Filter mean (i.e., exponential of natural parameter) trend for 0th-order
 * unweighted Poisson trend filtering
 * @param y signal length
 * @param lam hyperparameter
 * @return filtered mean (exponential of natural parameter) trend
 */
// [[Rcpp::export]]
Rcpp::NumericVector dptf(Rcpp::NumericVector y, double lam) {
  int n = y.size();
  Rcpp::NumericVector beta(n);
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
Rcpp::NumericVector dptf_past(Rcpp::NumericVector y,
                              double lam,
                              Rcpp::NumericVector w) {
  int n = y.size();
  Rcpp::NumericVector beta(n);
  tf_dp_past(n, y.begin(), w.begin(), lam, beta.begin());
  return beta;
}
