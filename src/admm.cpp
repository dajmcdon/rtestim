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

void admm(int M,
          arma::vec const& y,
          arma::vec const& x,
          arma::vec const& w,
          int n,
          int ord,
          arma::vec& theta,
          arma::vec& z,
          arma::vec& u,
          double lambda,
          double rho,
          double mu,
          double tol,
          int& iter) {
  double r_norm, s_norm;
  vec z_old = z;
  double lam_z = lambda / rho;
  vec c(n);  // a buffer
  vec c2(n);
  vec c3(n);
  vec c4(z.size());

  // start of iteration:
  for (iter = 0; iter < M; iter++) {
    if (iter % 1000 == 0)
      Rcpp::checkUserInterrupt();
    // update primal variable:
    calcDTDvline(n, ord, x, theta, c2);  // c2 = DD * theta
    z -= u;
    calcDTvline(n, ord, x, z, c3);  // c3 = Dt * (z - u)
    c = y / n - rho * c2 + rho * c3 + mu * theta;
    c = c / mu + log(w);
    c.transform([&](double c) { return update_pois(c, mu, n); });
    theta = c - log(w);

    // update alternating variable:
    calcDvline(n, ord, x, theta, c4);
    c4 += u;  // c4 = D * theta + u;
    z = dptf(c4, lam_z);

    // update dual variable:
    calcDvline(n, ord, x, theta, c4);  // c4 = D * theta
    u += c4 - z;

    // stopping criteria check:
    r_norm = sqrt(mean(square(c4 - z)));
    // dual residuals:
    s_norm = rho * sqrt(mean(square(z_old - z)));

    if ((r_norm < tol) && (s_norm < tol))
      break;
    // auxiliary variables update:
    z_old = z;
  }
}

// This is a wrapper around the void function to use in test_that()
// [[Rcpp::export]]
List admm_testing(int M,
                  arma::vec const& y,
                  arma::vec const& x,
                  arma::vec const& w,
                  int n,
                  int ord,
                  arma::vec theta,
                  arma::vec z,
                  arma::vec u,
                  double lambda,
                  double rho,
                  double mu,
                  double tol,
                  int iter) {
  admm(M, y, w, n, theta, z, u, lambda, rho, mu, DD, D, tol, iter);
  List out =
      List::create(Named("y") = y, Named("n") = n, Named("lambda") = lambda,
                   Named("theta") = exp(theta), Named("z") = z, Named("u") = u,
                   Named("niter") = iter + 1);
  return out;
}

// [[Rcpp::export]]
arma::vec admm_gauss(int M,
                     int n,
                     int ord,
                     arma::vec const& y,
                     arma::vec const& x,
                     arma::vec const& w,
                     arma::vec& theta,
                     arma::vec& z,
                     arma::vec& u,
                     double rho,
                     double lam_z,
                     double r_norm,
                     double s_norm,
                     arma::sp_mat const& DD,
                     double tol) {
  mat dDD(DD);
  dDD *= n * rho;
  mat W(n, n);
  vec c(n);
  vec c2(z.size());
  vec z_old(z);

  for (int iter = 0; iter < M; iter++) {
    if (iter % 1000 == 0)
      Rcpp::checkUserInterrupt();
    // solve for primal variable - theta:
    W = dDD;
    W.diag() += exp(theta);
    z -= u;
    calcDTvline(n, ord, x, z, c);  // c = Dt * (z - u)
    theta = solve(W, exp(theta) % y + n * rho * c);
    // solve for alternating variable - z:
    calcDvline(n, ord, x, theta, c2);
    c2 += u;  // c2 = D * theta + u;
    z = dptf(c2, lam_z);
    // update dual variable - u:
    calcDvline(n, ord, x, theta, c2);
    u += c2 - z;  // u += D * theta - z;

    // primal residuals:
    r_norm = sqrt(mean(square(c2 - z)));
    // dual residuals:
    z_old -= z;
    calcDTvline(n, ord, x, z_old, c);  // c = Dt * (z_old - z);
    s_norm = rho * sqrt(mean(square(c)));
    // stopping criteria check:
    if (r_norm < tol && s_norm < tol)
      break;

    // auxiliary variables update:
    z_old = z;
  }
  return theta;
}

void irls_admm(int M,
               int n,
               int ord,
               arma::vec const& y,
               arma::vec const& x,
               arma::vec const& w,
               arma::vec& theta,
               arma::vec& z,
               arma::vec& u,
               double lambda,
               double rho,
               double mu,
               double alpha,
               double gamma,
               arma::sp_mat const& D,
               double tol,
               int& iter) {
  double s;             // step size
  vec obj_list(M + 1);  // objective list for each iterate
  double obj = 1e4;     // initialize it to be large
  vec theta_old(n);     // a buffer for line search
  double lam_z = lambda / rho;
  sp_mat const DD = D.t() * D;
  int m = z.size();
  vec c1(n);  // a buffer
  vec c2(m);  // a buffer
  double r_norm = 0.0;
  double s_norm = 0.0;

  // initialize objective
  vec Dv(m);
  obj = pois_obj(y, w, theta, lambda, Dv);
  obj_list[0] = obj;
  int iter_best = 0;

  for (iter = 0; iter < M; iter++) {
    if (iter % 1000 == 0)
      Rcpp::checkUserInterrupt();

    theta_old = theta;

    // define new(fake) data for least squares problem
    c1 = fake_data(y, w, theta);
    // solve least squares problem (Gaussian TF)
    theta = admm_gauss(M, n, ord, c1, x, w, theta, z, u, rho, lam_z, r_norm,
                       s_norm, DD, tol);

    // line search for step size
    s = line_search(s, lambda, alpha, gamma, y, x, w, n, ord, theta, theta_old,
                    c1, c2, M);
    if (s < 0)
      break;
    // update theta
    theta *= s;
    theta += (1 - s) * theta_old;

    // compute objective
    Dv = D * theta;
    obj = pois_obj(y, w, theta, lambda, Dv);
    obj_list[iter + 1] = obj;
    // check stopping criteria
    if (obj < obj_list[iter_best]) {
      if (obj_list[iter_best] - obj <= fabs(obj_list[iter_best]) * tol) {
        obj_list[iter_best] = obj;
        iter_best = iter + 1;
        break;
      }
      obj_list[iter_best] = obj;
      iter_best = iter + 1;
    }
    if (iter >= iter_best + 4 && iter_best != 0)  // adjust the iterate steps
      break;
  }
}
