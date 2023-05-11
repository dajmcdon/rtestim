#include <RcppEigen.h>
#include <Eigen/Sparse>
#include <boost/math/special_functions/lambert_w.hpp>
#include "utils.h"
#include "admm.h"
#include "dptf.h"


typedef Eigen::COLAMDOrdering<int> Ord;

using Eigen::SparseMatrix;
using Eigen::SparseQR;
using Eigen::VectorXd;
SparseQR<SparseMatrix<double>, Ord> qradmm;

using namespace Rcpp;

/**
 * Solving Lambert_0 function for primal step of linearized ADMM
 * @param c primal variable to be computed
 * @param mu upper bound for the primal step
 * @param n signal length
 * @return updated primal variable
 */
double update_pois(double c, double mu, int n) {
  // solve x such that exp(x) * x = exp(c) / (mu * n) which is a Lambert_0
  // function; then let theta = c - x
  if (c < 500) {
    double cz = exp(c) / (mu * n);
    c -= boost::math::lambert_w0(cz);
  } else {  // deal with potential overflow from big exp(c)
    // See:
    // https://en.wikipedia.org/wiki/Lambert_W_function#Asymptotic_expansions
    // We use the first four terms.
    double la, lb;
    la = c - log(mu * n);
    lb = log(la);
    c -= la - lb + lb / la + (lb * (lb - 2)) / (2 * la * la);
  }
  return c;
}

void linear_admm(int& M,
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
    if (iter % 100 == 0)
      Rcpp::checkUserInterrupt();
    // update primal variable:
    tmp_n = doDtDv(theta, ord, x);
    tmp_n = y / n - rho * tmp_n + rho * doDtv(z - u, ord, x) + mu * theta;
    tmp_n = tmp_n / mu + log(w);
    for (int it = 0; it < n; it++)
      tmp_n(it) = update_pois(tmp_n(it), mu, n);
    theta = tmp_n - log(w);

    // update alternating variable:
    tmp_m = doDv(theta, ord, x);
    tmp_m += u;
    z = dptf(tmp_m, lam_z);

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

/**
 * This is a wrapper around the void function to use in test_that()
 */
// [[Rcpp::export]]
Rcpp::List linear_admm_testing(int M,
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
  linear_admm(M, y, x, w, n, ord, theta, z, u, lambda, rho, mu, tol, iter);
  List out =
      List::create(Named("y") = y, Named("n") = n, Named("lambda") = lambda,
                   Named("theta") = exp(theta), Named("z") = z, Named("u") = u,
                   Named("niter") = iter);
  return out;
}

/**
 * ADMM for Gaussian trend filtering
 * @param M maximum iteration of the algos
 * @param n signal length
 * @param ord degree of Poisson trend filtering
 * @param y observed signals
 * @param x signal locations
 * @param w signal weights
 * @param theta primal variable of length `n`
 * @param z auxiliary variable of length `n-ord`
 * @param u dual variable of length `n-ord`
 * @param rho Lagrangian parameter of ADMM
 * @param lam_z hyperparameter of the auxiliary step of ADMM
 * @param r_norm primal residual
 * @param s_norm dual residual
 * @param DD D^T * D
 * @param tol tolerance of stopping criteria
 * @return updated primal variable
 */
// [[Rcpp::export]]
Rcpp::NumericVector admm_gauss(int M,
                               int n,
                               int ord,
                               Rcpp::NumericVector const& y,
                               Rcpp::NumericVector const& x,
                               Rcpp::NumericVector const& w,
                               Rcpp::NumericVector& theta,
                               Rcpp::NumericVector& z,
                               Rcpp::NumericVector& u,
                               double rho,
                               double lam_z,
                               double r_norm,
                               double s_norm,
                               Eigen::SparseMatrix<double> const& DD,
                               double tol) {
  NumericVector z_old = clone(z);
  NumericVector tmp_n(n);
  NumericVector tmp_m(z.size());
  NumericVector Dth(z.size());
  VectorXd tmp_theta(n);
  VectorXd etheta = nvec_to_evec(theta);
  SparseMatrix<double> cDD = DD * n * rho;  // a copy that doesn't change
  // cDD.diagonal() += etheta.exp();
  for (int i = 0; i < n; i++) {
    cDD.diagonal()(i) += exp(etheta(i));
  }
  qradmm.compute(cDD);
  NumericVector W = exp(theta);  // fix the weights

  for (int iter = 0; iter < M; iter++) {
    if (iter % 50 == 0)
      Rcpp::checkUserInterrupt();
    // solve for primal variable - theta:
    tmp_n = doDtv(z - u, ord, x) * n * rho;
    tmp_n += W * y;
    tmp_theta = nvec_to_evec(tmp_n);
    tmp_theta = qradmm.solve(tmp_theta);
    theta = evec_to_nvec(tmp_theta);
    // solve for alternating variable - z:
    Dth = doDv(theta, ord, x);
    tmp_m = Dth + u;
    z = dptf(tmp_m, lam_z);

    // update dual variable - u:
    u += Dth - z;

    // primal residuals:
    r_norm = sqrt(mean(pow(Dth - z, 2)));
    // dual residuals:
    tmp_n = doDtv(z - z_old, ord, x);
    s_norm = rho * sqrt(mean(pow(tmp_n, 2)));
    // stopping criteria check:
    if (r_norm < tol && s_norm < tol)
      break;

    // auxiliary variables update:
    z_old = z;
  }
  return theta;
}

void prox_newton(int& M,
                 int Minner,
                 int n,
                 int ord,
                 Rcpp::NumericVector const& y,
                 Rcpp::NumericVector const& x,
                 Rcpp::NumericVector const& w,
                 Rcpp::NumericVector& theta,
                 Rcpp::NumericVector& z,
                 Rcpp::NumericVector& u,
                 double lambda,
                 double rho,
                 double alpha,
                 double gamma,
                 Eigen::SparseMatrix<double> const& DD,
                 double tol,
                 int Mline,
                 int& iter) {
  double s;                       // step size
  NumericVector obj_list(M + 1);  // objective list for each iterate
  double obj = 1e4;               // initialize it to be large
  NumericVector theta_old(n);     // a buffer for line search
  double lam_z = lambda / rho;
  double r_norm = 0.0;
  double s_norm = 0.0;

  // initialize objective
  obj = pois_obj(ord, y, x, w, theta, lambda);
  obj_list(0) = obj;
  int iter_best = 0;

  for (iter = 0; iter < M; iter++) {
    if (iter % 50 == 0)
      Rcpp::checkUserInterrupt();
    theta_old = theta;

    // define new(Gaussianized) data for least squares problem
    NumericVector std_y = gaussianized_data(y, w, theta);
    // solve least squares problem (Gaussian TF)
    theta = admm_gauss(Minner, n, ord, std_y, x, w, theta, z, u, rho, lam_z,
                       r_norm, s_norm, DD, tol);

    // line search for step size
    s = line_search(s, lambda, alpha, gamma, y, x, w, n, ord, theta, theta_old,
                    Mline);
    if (s < 0)
      break;
    // update theta
    theta = theta * s + (1 - s) * theta_old;

    // compute objective
    obj = pois_obj(ord, y, x, w, theta, lambda);
    obj_list(iter + 1) = obj;
    // check stopping criteria
    if (obj < obj_list(iter_best)) {
      if (obj_list(iter_best) - obj <= abs(obj_list(iter_best)) * tol) {
        obj_list(iter_best) = obj;
        iter_best = iter + 1;
        break;
      }
      obj_list(iter_best) = obj;
      iter_best = iter + 1;
    }
    // adjust the iterate steps
    if (iter >= iter_best + 4 && iter_best != 0)
      break;
  }
}

/**
 * This is a wrapper around the void function to use in test_that()
 */
// [[Rcpp::export]]
Rcpp::List prox_newton_testing(int M,
                               int Minner,
                               int n,
                               int ord,
                               Rcpp::NumericVector const& y,
                               Rcpp::NumericVector const& x,
                               Rcpp::NumericVector const& w,
                               Rcpp::NumericVector& theta,
                               Rcpp::NumericVector& z,
                               Rcpp::NumericVector& u,
                               double lambda,
                               double rho,
                               double alpha,
                               double gamma,
                               Eigen::SparseMatrix<double> const& DD,
                               double tol,
                               int Mline,
                               int iter) {
  prox_newton(M, Minner, n, ord, y, x, w, theta, z, u, lambda, rho, alpha,
              gamma, DD, tol, Mline, iter);
  List out =
      List::create(Named("y") = y, Named("n") = n, Named("lambda") = lambda,
                   Named("theta") = exp(theta), Named("z") = z, Named("u") = u,
                   Named("niter") = iter);
  return out;
}
