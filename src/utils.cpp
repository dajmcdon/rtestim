#include <RcppArmadillo.h>
#include <cmath>
#include "utils.h"

using namespace Rcpp;

double chooseC(int n, int k) {
  return Rf_choose(n, k);
}


// [[Rcpp::export]]
arma::sp_mat buildD(int n, int ord) {
  int m = n - ord - 1;
  arma::sp_mat D(m, n);
  int c1 = ord + 1;
  double z;
  for (int i = 0; i <= c1; i++) {
    z = chooseC(c1, i) * std::pow(-1, i + ord + 1);
    D.diag(i) += z;
  }
  return D;
}

// [[Rcpp::export]]
arma::sp_mat buildDx(int k, const arma::vec& x) {
  // x should be in increasing order
  int n = x.size();
  arma::sp_mat Dmat(n - 1, n); // output
  arma::sp_mat D1(n - 1, n);
  D1.diag(0) -= 1;
  D1.diag(1) += 1;
  if (k == 0) return D1; // k=0 is the same as usual
  Dmat += D1;

  for (int tk = 1; tk <= k; tk++) {
    arma::sp_mat D1_cp = arma::resize(D1, n - tk - 1, n - tk);
    arma::sp_mat delx(n - tk, n - tk);
    arma::vec lagdiff = tk / (x.tail(n - tk) - x.head(n - tk));
    delx.diag(0) = lagdiff;
    D1_cp = D1_cp * delx;
    Dmat = D1_cp * Dmat;
  }
  return Dmat;
}

// [[Rcpp::export]]
arma::sp_mat buildDx_tilde(int k, const arma::vec& x) {
  int n = x.size();
  if (k == 0) { // Only useful for k >= 1
    arma::sp_mat Dk(n, n);
    Dk.eye();
    return(Dk);
  }
  arma::sp_mat Dk = buildDx(k - 1, x);
  arma::vec lagdiff = k / (x.tail(n - k) - x.head(n - k));
  arma::sp_mat delx(n - k, n - k);
  delx.diag(0) = lagdiff;
  Dk = delx * Dk;
  return(Dk);
}

void create_lambda(arma::vec& lambda,
                   double lambdamin,
                   double lambdamax,
                   double lambda_min_ratio,
                   int nsol) {
  if (lambda.size() > 0) {
    lambdamin = lambda.min();
    lambdamax = lambda.max();
    nsol = lambda.size();
  } else {
    lambda.set_size(nsol);
    lambdamin = (lambdamin < 0) ? lambda_min_ratio * lambdamax : lambdamin;
    int sm_lam = lambdamin < 0;
    if (sm_lam == 1) {
      double lmpad = 1e-8;
      lambda.tail(nsol - 1) = arma::logspace(log10(lmpad), log10(lambdamax), nsol - 1);
      lambda(0) = 0;
    } else {
      lambda = arma::logspace(log10(lambdamin), log10(lambdamax), nsol);
    }
  }
}


void EntrywiseSoftThreshold(arma::vec& z, double lam) {
  z.transform( [&](double x) {
    return ((x > 0) - (x < 0)) * std::max(0.0, std::fabs(x) - lam);
  } );
}
