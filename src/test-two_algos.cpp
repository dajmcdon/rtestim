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
  int n = 10;
  int M = 1e3;
  int M_inner = 3;
  int iter = 0;
  int ord;
  double lambda = 0.1;
  double rho = 0.1;
  double alpha = 0.5;
  double gamma = 0.9;
  double tol = 1e-4;
  double mu;
  vec y;
  vec x = arma::linspace(1, n, n);
  vec w = ones(n);
  vec theta1(n);
  vec theta2(n);
  vec z1;
  vec u1;
  vec z2;
  vec u2;
  sp_mat D;

  y = {9, 11, 10, 5, 12, 9, 10, 9, 9, 13};
  y.print("y:");

  /* k = 0 */
  test_that("two algo have same estimates when k = 0") {
    ord = 0;
    D = buildDx_tilde(n, ord, x);
    mu = 2 * pow(4, ord);
    z1 = zeros(n - ord);
    z2 = zeros(n - ord);
    u1 = zeros(n - ord);
    u2 = zeros(n - ord);
    theta1 = zeros(n);
    theta2 = zeros(n);
    iter = 0;

    linear_admm(M, y, x, w, n, ord, theta1, z1, u1, lambda, rho, mu, tol, iter);
    prox_newton(M, M_inner, n, ord, y, x, w, theta2, z2, u2, lambda, rho, alpha,
                gamma, D, tol, iter);
    exp(theta1).brief_print("est1:");
    exp(theta2).brief_print("est2:");
    expect_true(mean(abs(exp(theta1) - exp(theta2))) <= .1);
  }

  /* k = 1 */
  test_that("two algo have same estimates when k = 1") {
    ord = 1;
    D = buildDx_tilde(n, ord, x);
    mu = 2 * pow(4, ord);
    z1 = zeros(n - ord);
    z2 = zeros(n - ord);
    u1 = zeros(n - ord);
    u2 = zeros(n - ord);
    theta1 = zeros(n);
    theta2 = zeros(n);
    iter = 0;
    linear_admm(M, y, x, w, n, ord, theta1, z1, u1, lambda, rho, mu, tol, iter);
    prox_newton(M, M_inner, n, ord, y, x, w, theta2, z2, u2, lambda, rho, alpha,
                gamma, D, tol, iter);
    expect_true(mean(abs(exp(theta1) - exp(theta2))) <= .1);
  }

  /* k = 2 */
  test_that("two algo have same estimates when k = 2") {
    ord = 2;
    D = buildDx_tilde(n, ord, x);
    mu = 2 * pow(4, ord);
    z1 = zeros(n - ord);
    z2 = zeros(n - ord);
    u1 = zeros(n - ord);
    u2 = zeros(n - ord);
    theta1 = zeros(n);
    theta2 = zeros(n);
    iter = 0;
    linear_admm(M, y, x, w, n, ord, theta1, z1, u1, lambda, rho, mu, tol, iter);
    prox_newton(M, M_inner, n, ord, y, x, w, theta2, z2, u2, lambda, rho, alpha,
                gamma, D, tol, iter);
    exp(theta1).brief_print("est1:");
    exp(theta2).brief_print("est2:");
    expect_true(mean(abs(exp(theta1) - exp(theta2))) <= .1);
  }
}
