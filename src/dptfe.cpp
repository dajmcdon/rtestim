#include <Rcpp.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include "tf_dp.h"
#include "dptfe.h"



// [[Rcpp::export]]
Rcpp::NumericVector dptfe(Rcpp::NumericVector y, double lam) {
  int n = y.size();
  Rcpp::NumericVector beta(n);
  tf_dp(n, y.begin(), lam, beta.begin());  // get pointers to the vectors
  return beta;
}

// [[Rcpp::export]]
Rcpp::NumericVector dptfe_past(Rcpp::NumericVector y,
                               double lam,
                               Rcpp::NumericVector w) {
  int n = y.size();
  Rcpp::NumericVector beta(n);
  tf_dp_past(n, y.begin(), w.begin(), lam, beta.begin());
  return beta;
}
