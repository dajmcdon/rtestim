/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 */
#include <testthat.h>
#include <RcppArmadillo.h>
#include "utils.h"
#include "admm.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

context("Check equivalence of two algorithms for k=0,1,2") {
  int n = 100;
  int M = 1e3;
  int M_inner = 3;
  int iter;
  int ord;
  double lambda = 0.1;
  double rho = 0.1;
  double alpha = 0.5;
  double gamma = 0.9;
  double tol = 1e-4;
  double mu;
  vec y;
  vec x = arma::linspace(1, n);
  vec w = ones(n);
  vec theta1(n);
  vec theta2(n);
  vec z;
  vec u;
  sp_mat D;

  /* k = 0 */
  ord = 0;
  D = buildDx_tilde(n, ord, x);
  mu = 2 * pow(4, ord);
  z = u = zeros(n - ord);
  for (int i = 0; i < 1; i++) {
    y = Rcpp::rpois(n, 10);
    y.brief_print("y:");
    // y.brief_print("y:");
    linear_admm(M, y, x, w, n, ord, theta1, z, u, lambda, rho, mu, tol, iter);
    prox_newton(M, M_inner, n, ord, y, x, w, theta2, z, u, lambda, rho, alpha,
                gamma, D, tol, iter);
    exp(theta1).brief_print("est1:");
    exp(theta2).brief_print("est2:");
    test_that("two algo have same estimates when k = 0") {
      expect_true(mean(abs(exp(theta1) - exp(theta2))) <= 1);
    }
  }

  /* k = 1 */
  ord = 1;
  D = buildDx_tilde(n, ord, x);
  mu = 2 * pow(4, ord);
  for (int i = 0; i < 1; i++) {
    y = Rcpp::rpois(n, 10);
    linear_admm(M, y, x, w, n, ord, theta1, z, u, lambda, rho, mu, tol, iter);
    prox_newton(M, M_inner, n, ord, y, x, w, theta2, z, u, lambda, rho, alpha,
                gamma, D, tol, iter);
    test_that("two algo have same estimates when k = 1") {
      expect_true(mean(abs(exp(theta1) - exp(theta2))) <= 1);
    }
  }

  /* k = 2 */
  ord = 2;
  D = buildDx_tilde(n, ord, x);
  mu = 2 * pow(4, ord);
  for (int i = 0; i < 1; i++) {
    y = Rcpp::rpois(n, 10);
    linear_admm(M, y, x, w, n, ord, theta1, z, u, lambda, rho, mu, tol, iter);
    prox_newton(M, M_inner, n, ord, y, x, w, theta2, z, u, lambda, rho, alpha,
                gamma, D, tol, iter);
    test_that("two algo have same estimates when k = 2") {
      expect_true(mean(abs(exp(theta1) - exp(theta2))) <= 1);
    }
  }
}
