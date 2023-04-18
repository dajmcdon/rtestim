/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 */
#include <testthat.h>
#include <RcppEigen.h>
#include <Eigen/Sparse>
#include "utils.h"
#include "admm.h"

typedef Eigen::COLAMDOrdering<int> Ord;

using Eigen::SparseMatrix;
using Eigen::SparseQR;
using Eigen::VectorXd;


using namespace Rcpp;

context("Check equivalence of two algorithms for k=0,1,2") {

  /* k = 0 */
  test_that("two algo have same estimates when k = 0") {
    int n = 10;
    int M = 1e3;
    int M_inner = 3;
    int Mline = 5;
    int iter = 0;
    int ord;
    double lambda = 0.1;
    double rho = 0.1;
    double alpha = 0.5;
    double gamma = 0.9;
    double tol = 1e-4;
    double mu;
    NumericVector y;
    NumericVector x(n);  // = arma::linspace(1, n, n);
    for (int i = 0; i < n; i++) {
      x(i) = i + 1;
    }
    NumericVector w(n);
    for (int i = 0; i < n; i++) {
      w(i) = 1;
    }
    Eigen::SparseMatrix<double> D;
    Eigen::SparseMatrix<double> DD;
    y = {9, 11, 10, 5, 12, 9, 10, 9, 9, 13};
    Rcout << "Data is " << y << "\n";

    ord = 0;
    D = get_Dtil(ord, x);
    DD = D.transpose() * D;
    mu = 2 * pow(4, ord+1); // or `ord`?
    iter = 0;
    NumericVector theta1(n);
    NumericVector z1(n - ord);
    NumericVector u1(n - ord);
    linear_admm(M, y, x, w, n, ord, theta1, z1, u1, lambda, rho, mu, tol, iter);
    Rcout << "The estimated theta of linear_admm is : " << theta1 << "\n";

    iter = 0;
    NumericVector theta2(n);
    NumericVector z2(n - ord);
    NumericVector u2(n - ord);
    prox_newton(M, M_inner, n, ord, y, x, w, theta2, z2, u2, lambda, rho, alpha,
                gamma, DD, tol, Mline, iter);
    Rcout << "The estimated theta of prox_newton is : " << theta2 << "\n";

    expect_true(mean(abs(exp(theta1) - exp(theta2))) <= .1);
  }

  /* k = 1 */
  test_that("two algo have same estimates when k = 1") {
    int n = 10;
    int M = 1e3;
    int M_inner = 3;
    int Mline = 5;
    int iter = 0;
    int ord;
    double lambda = 0.1;
    double rho = 0.1;
    double alpha = 0.5;
    double gamma = 0.9;
    double tol = 1e-4;
    double mu;
    NumericVector y;
    NumericVector x(n);
    for (int i = 0; i < n; i++) {
      x(i) = i + 1;
    }
    NumericVector w(n);
    for (int i = 0; i < n; i++) {
      w(i) = 1;
    }
    Eigen::SparseMatrix<double> D;
    Eigen::SparseMatrix<double> DD;

    y = {9, 11, 10, 5, 12, 9, 10, 9, 9, 13};

    ord = 1;
    D = get_Dtil(ord, x);
    DD = D.transpose() * D;
    mu = 2 * pow(4, ord);
    iter = 0;
    NumericVector theta1(n);
    NumericVector z1(n - ord);
    NumericVector u1(n - ord);
    linear_admm(M, y, x, w, n, ord, theta1, z1, u1, lambda, rho, mu, tol, iter);
    Rcout << "The estimated theta of linear_admm is : " << theta1 << "\n";

    iter = 0;
    NumericVector theta2(n);
    NumericVector z2(n - ord);
    NumericVector u2(n - ord);
    prox_newton(M, M_inner, n, ord, y, x, w, theta2, z2, u2, lambda, rho, alpha,
                gamma, DD, tol, Mline, iter);
    Rcout << "The estimated theta of prox_newton is : " << theta2 << "\n";

    expect_true(mean(abs(exp(theta1) - exp(theta2))) <= .1);
  }

  /* k = 2 */
  test_that("two algo have same estimates when k = 2") {
    int n = 10;
    int M = 1e3;
    int M_inner = 3;
    int Mline = 5;
    int iter = 0;
    int ord;
    double lambda = 0.1;
    double rho = 0.1;
    double alpha = 0.5;
    double gamma = 0.9;
    double tol = 1e-4;
    double mu;
    NumericVector y;
    NumericVector x(n);  // = arma::linspace(1, n, n);
    for (int i = 0; i < n; i++) {
      x(i) = i + 1;
    }
    NumericVector w(n);
    for (int i = 0; i < n; i++) {
      w(i) = 1;
    }

    Eigen::SparseMatrix<double> D;
    Eigen::SparseMatrix<double> DD;

    y = {9, 11, 10, 5, 12, 9, 10, 9, 9, 13};
    Rcout << "Data is " << y << "\n";
    ord = 2;
    D = get_Dtil(ord, x);
    DD = D.transpose() * D;
    mu = 2 * pow(4, ord);
    iter = 0;
    NumericVector theta1(n);
    NumericVector z1(n - ord);
    NumericVector u1(n - ord);
    linear_admm(M, y, x, w, n, ord, theta1, z1, u1, lambda, rho, mu, tol, iter);
    Rcout << "The estimated theta of linear_admm is : " << theta1 << "\n";

    iter = 0;
    NumericVector theta2(n);
    NumericVector z2(n - ord);
    NumericVector u2(n - ord);
    prox_newton(M, M_inner, n, ord, y, x, w, theta2, z2, u2, lambda, rho, alpha,
                gamma, DD, tol, Mline, iter);
    Rcout << "The estimated theta of prox_newton is : " << theta2 << "\n";

    expect_true(mean(abs(exp(theta1) - exp(theta2))) <= .1);
  }
}
