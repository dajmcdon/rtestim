#include <RcppArmadillo.h>
#include <boost/math/special_functions/lambert_w.hpp>
#include "dptf.h"
#include "utils.h"
#include "admm.h"


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;
using namespace arma;

double update_pois(double c, double mu, int n) {
  // @elvis, explain why this is the solution
  if (c < 500) { // deal with potential overflow from big exp(c)
    c -= boost::math::lambert_w0(exp(c) / (mu * n));
  } else {
    // See: https://en.wikipedia.org/wiki/Lambert_W_function#Asymptotic_expansions
    // We use the first four terms.
    double la, lb;
    la = c - log(mu * n);
    lb = log(la);
    c -= la - lb + lb / la + (lb * (lb - 2)) / (2 * la * la);
  }
  return c;
}


void admm(int M,
          arma::vec const& y,
          arma::vec const& w,
          int n,
          arma::vec& theta,
          arma::vec& z,
          arma::vec& u,
          double lambda,
          double rho,
          double mu,
          arma::sp_mat const& DD,
          arma::sp_mat const& D,
          double tol) {

  int iter;
  double r_norm, s_norm;
  vec z_old = z;
  sp_mat Dt = D.t();
  double lam_z = lambda / rho;
  vec c;
  vec beta = z;

  // start of iteration:
  for (iter = 0; iter < M; iter++) {
    if (iter % 1000 == 0) Rcpp::checkUserInterrupt();
    // update primal variable:
    c = y / n - rho * DD * theta + rho * Dt * (z - u) + mu * theta;
    c = c / mu + log(w);
    c.transform([&](double c) { return update_pois(c, mu, n); });
    theta = c - log(w);

    // update alternating variable:
    c = D * theta + u;
    //    beta.transform([*](double beta) { return tf_dp(n, c, lam_z, beta); });
    //    z = beta;
    z = dptf(c, lam_z);

    // update dual variable:
    u += D * theta - z;

    // stopping criteria check:
    r_norm = sqrt(mean(square(D * theta - z)));
    // dual residuals:
    s_norm = sqrt(mean(square(z_old - z)));

    if ((r_norm < tol) && (s_norm < tol)) {
      iter++;
      break;
    }
    // auxiliary variables update:
    z_old = z;
  }
}

// This is a wrapper around the void function to use in test_that()
// [[Rcpp::export]]
List admm_testing(int M,
                  arma::vec const& y,
                  arma::vec const& w,
                  int n,
                  arma::vec theta,
                  arma::vec z,
                  arma::vec u,
                  double lambda,
                  double rho,
                  double mu,
                  arma::sp_mat const& DD,
                  arma::sp_mat const& D,
                  double tol) {
  admm(M, y, w, n, theta, z, u, lambda, rho, mu, DD, D, tol);
  List out = List::create(
    Named("y") = y,
    Named("n") = n,
    Named("lambda") = lambda,
    Named("theta") = exp(theta),
    Named("z") = z,
    Named("u") = u,
    Named("niter") = M
  );
  return out;
}

