#include <Eigen/Sparse>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <dspline.h>
#include "utils-eigen.h"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(dspline)]]

using namespace Rcpp;

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

// [[Rcpp::export]]
NumericVector doDv(NumericVector v, int k, NumericVector xd) {
  return dspline::rcpp_d_mat_mult(v, k, xd, false, false);
}

// [[Rcpp::export]]
NumericVector doDtv(NumericVector v, int k, NumericVector xd) {
  return dspline::rcpp_d_mat_mult(v, k, xd, false, true);
}

// [[Rcpp::export]]
NumericVector doDtDv(NumericVector v, int k, NumericVector xd) {
  NumericVector tmp = doDv(v, k, xd);
  return doDtv(tmp, k, xd);
}

// [[Rcpp::export]]
double one_norm(NumericVector const& z) {
  double out = sum(abs(z));
  return out;
}

// [[Rcpp::export]]
double pois_obj(int ord,
                NumericVector const& y,
                NumericVector const& x,
                NumericVector const& w,
                NumericVector& theta,
                double lambda) {
  NumericVector Dtheta = doDv(theta, ord, x);
  NumericVector v = -y * theta + w * exp(theta);
  double obj = mean(v) + lambda * sum(abs(Dtheta));
  return obj;
}

// [[Rcpp::export]]
NumericVector gaussianized_data(NumericVector const& y,
                                NumericVector const& w,
                                NumericVector& theta) {
  int n = y.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    if (w(i) * exp(theta(i)) > 1e-3) {
      out(i) = y(i) * exp(-theta(i)) / w(i) - 1 + theta(i);
    } else {  // deal with overflow using approximation
      out(i) = y(i) - exp(theta(i)) / w(i) + theta(i);
    }
  }
  return out;
}

// [[Rcpp::export]]
double line_search(double s,
                   double lambda,
                   double alpha,
                   double gamma,
                   NumericVector const& y,
                   NumericVector const& x,
                   NumericVector const& w,
                   int n,
                   int ord,
                   NumericVector& theta,
                   NumericVector& theta_old,
                   int M) {
  NumericVector gradient(n);
  double grades;
  NumericVector dir = theta - theta_old;

  // initialize upper bound
  double bound = mean(dir * (w * exp(theta) - y));

  NumericVector Dth = doDv(theta, ord, x);
  NumericVector Dth_old = doDv(theta_old, ord, x);
  double Dth_old_norm = one_norm(Dth_old);
  bound += lambda * (one_norm(Dth) - Dth_old_norm);
  NumericVector exp_theta_old = exp(theta_old);

  s = 1;
  for (int i = 0; i < M; i++) {
    // compute gradient/ grades
    gradient = -s * dir * y + w * exp(theta_old + s * dir) - w * exp_theta_old;
    if (i > 0) Dth = doDv(dir, ord, x) * s + Dth_old;
    grades = mean(gradient) + lambda * (one_norm(Dth) - Dth_old_norm);

    // adjust upper bound
    bound *= alpha * s;

    // check criteria
    if (grades <= bound) break;
    else s *= gamma;
  }
  return s;
}
