#include <RcppArmadillo.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include "tf_dp.h"
#include "dptf.h"
// [[Rcpp::depends(RcppArmadillo)]]

void tf_dp(int n, double* y, double lam, double* beta);

// [[Rcpp::export]]
arma::vec dptf(arma::vec y, double lam) {
  int n = y.size();
  arma::vec beta(n);

  tf_dp(n, y.begin(), lam, beta.begin());  // get pointers to the vectors
  return beta;
}
