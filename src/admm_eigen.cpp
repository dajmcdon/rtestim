#include <RcppEigen.h>
#include <boost/math/special_functions/lambert_w.hpp>
#include "dptfe.h"
#include "utils-eigen.h"


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;


double update_epois(double c, double mu, int n) {
  // @elvis, explain why this is the solution
  if (c < 500) {  // deal with potential overflow from big exp(c)
    c -= boost::math::lambert_w0(exp(c) / (mu * n));
  } else {
    // See:
    // https://en.wikipedia.org/wiki/Lambert_W_function#Asymptotic_expansions We
    // use the first four terms.
    double la, lb;
    la = c - log(mu * n);
    lb = log(la);
    c -= la - lb + lb / la + (lb * (lb - 2)) / (2 * la * la);
  }
  return c;
}



void admm_eigen(int M,
          NumericVector const& y,
          NumericVector const& x,
          NumericVector const& w,
          int n,
          int ord,
          NumericVector& theta,
          NumericVector& z,
          NumericVector& u,
          double lambda,
          double rho,
          double mu,
          double tol,
          int& iter) {
  double r_norm, s_norm;
  NumericVector z_old = clone(z);
  double lam_z = lambda / rho;
  NumericVector tmp_n(n);
  NumericVector tmp_m(z.size());

  // start of iteration:
  for (iter = 0; iter < M; iter++) {
    if (iter % 500 == 0) Rcpp::checkUserInterrupt();
    // update primal variable:
    tmp_n = doDtDv(theta, ord, x);
    tmp_n = y / n - rho * tmp_n + rho * doDtv(z - u, ord, x) + mu * theta;
    tmp_n = tmp_n / mu + log(w);
    for (int it = 0; it < n; it++) {
      tmp_n(it) = update_epois(tmp_n(it), mu, n);
    }
    theta = tmp_n - log(w);

    // update alternating variable:
    tmp_m = doDv(theta, ord, x);
    tmp_m += u;
    z = dptfe(tmp_m, lam_z);
    // z = dptf(tmp_m, lam_z);

    // update dual variable:
    tmp_m = doDv(theta, ord, x);
    u += tmp_m - z;

    // stopping criteria check:
    r_norm = sqrt(mean(pow(tmp_m - z, 2)));
    // dual residuals:
    s_norm = rho * sqrt(mean(pow(z_old - z, 2)));

    if ((r_norm < tol) && (s_norm < tol)) break;
    // auxiliary variables update:
    z_old = z;
  }
}

// This is a wrapper around the void function to use in test_that()
// [[Rcpp::export]]
List admm_eigen_testing(int M,
                  NumericVector const& y,
                  NumericVector const& x,
                  NumericVector const& w,
                  int n,
                  int ord,
                  NumericVector theta,
                  NumericVector z,
                  NumericVector u,
                  double lambda,
                  double rho,
                  double mu,
                  double tol,
                  int iter) {
  admm_eigen(M, y, x, w, n, ord, theta, z, u, lambda, rho, mu, tol, iter);
  List out =
    List::create(Named("y") = y, Named("n") = n, Named("lambda") = lambda,
                 Named("theta") = exp(theta), Named("z") = z, Named("u") = u,
                 Named("niter") = iter + 1);
  return out;
}
