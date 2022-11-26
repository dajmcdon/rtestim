#include <RcppArmadillo.h>
#include <boost/math/special_functions/lambert_w.hpp>
#include "dp.h"
#include "admm.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;
using namespace arma;

double update_primal(double c, double mu, int n) {
  c -= boost::math::lambert_w0(exp(c) / (mu * n));
  return c;
}

// [[Rcpp::export]]
Rcpp::List admm(int M,
                arma::vec y,
                arma::vec x,
                int n,
                arma::vec theta,
                arma::vec z,
                arma::vec u,
                double lambda,
                double rho,
                double mu,
                arma::sp_mat D,
                double tol) {
  int iter;
  double r_norm, s_norm;
  vec z_old = z;
  sp_mat Dt = D.t();
  sp_mat DD = Dt * D;
  double lam_z = lambda / rho;
  vec c;

  // start of iteration:
  for (iter = 0; iter < M; iter++) {
    // update primal variable:
    c = y / n - rho * DD * theta + rho * Dt * (z - u) + mu * theta;
    c = c / mu + log(x);
    c.transform([&](double c) { return update_primal(c, mu, n); });
    theta = c - log(x);

    // update alternating variable:
    c = D * theta + u;
    z = prox_dp_norm(c, lam_z);

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

  return List::create(Named("theta") = wrap(theta), Named("z") = wrap(z),
                      Named("u") = wrap(u), Named("prim_res") = r_norm,
                      Named("dual_res") = s_norm, Named("iter_num") = iter);
}
