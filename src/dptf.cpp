#include <RcppArmadillo.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include "tf_dp.h"
#include "dptf.h"
// [[Rcpp::depends(RcppArmadillo)]]

void tf_dp(int n, double* y, double lam, double* beta);
void tf_dp_weight(int n, double* y, double* w, double lam, double* beta);
void tf_dp_past_weight(int n,
                       double* y,
                       double* past,
                       double* w,
                       double lam,
                       double* beta);

// [[Rcpp::export]]
arma::vec dptf(arma::vec y, double lam) {
  int n = y.size();
  arma::vec beta(n);

  tf_dp(n, y.begin(), lam, beta.begin());  // get pointers to the vectors
  return beta;
}

// [[Rcpp::export]]
arma::vec dptf_weight(arma::vec y, double lam, arma::vec w) {
  int n = y.size();
  arma::vec beta(n);

  tf_dp_weight(n, y.begin(), w.begin(), lam, beta.begin());
  return beta;
}

// [[Rcpp::export]]
arma::vec dptf_past_weight(arma::vec y, double lam, arma::vec x, arma::vec w) {
  int n = y.size();
  arma::vec beta(n);

  tf_dp_past_weight(n, y.begin(), x.begin(), w.begin(), lam, beta.begin());
  return beta;
}
