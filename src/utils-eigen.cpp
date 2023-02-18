#include <Eigen/Sparse>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <dspline.h>
#include "utils-eigen.h"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(dspline)]]

using Rcpp::NumericVector;

Eigen::SparseMatrix<double> identity(int n) {
  Eigen::SparseMatrix<double> Id(n, n);
  Id.setIdentity();
  return Id;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> get_Dtil(int k, NumericVector xd) {
  int n = xd.size();
  return dspline::rcpp_b_mat(k, xd, false, Rcpp::seq(0, n-k-1), true);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> get_D(int k, NumericVector xd) {
  k++;
  int n = xd.size();
  return dspline::rcpp_b_mat(k, xd, true, Rcpp::seq(0, n-k-1), true);
}
