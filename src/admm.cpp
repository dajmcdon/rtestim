#include <RcppEigen.h>
#include <Eigen/Sparse>
#include "utils.h"
#include "dptf.h"
#include "admm.h"

typedef Eigen::COLAMDOrdering<int> Ord;

using Eigen::SparseMatrix;
using Eigen::SparseQR;
using Eigen::VectorXd;
SparseQR<SparseMatrix<double>, Ord> qradmm;

using namespace Rcpp;

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
 * @param DD D^T * D
 * @param tol tolerance of stopping criteria
 * @param iter interation index
 */
void admm_gauss(int M,
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
                Eigen::SparseMatrix<double> const& DD,
                double tol,
                int& iter) {
  double r_norm = 0.0;
  double s_norm = 0.0;
  NumericVector z_old = clone(z);
  NumericVector tmp_n(n);
  NumericVector tmp_m(z.size());
  NumericVector Dth(z.size());
  VectorXd tmp_theta(n);
  SparseMatrix<double> cDD = DD * n * rho;  // a copy that doesn't change
  NumericVector W = exp(theta) * w;         // update the weights
  VectorXd eW = nvec_to_evec(W);
  for (int i = 0; i < n; i++) {
    cDD.diagonal()(i) += eW(i);
  }
  qradmm.compute(cDD);

  for (iter = 0; iter < M; iter++) {
    if (iter % 1000 == 0)
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
}

void prox_newton(int M,
                 int& Minner,
                 int Mline,
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
                 int& total_iter) {
  double s;                       // step size
  NumericVector obj_list(M + 1);  // objective list for each iterate
  double obj = 1e4;               // initialize it to be large
  NumericVector theta_old(n);     // a buffer for line search
  double lam_z = lambda / rho;
  int inner_iter;

  // initialize objective
  obj = pois_obj(ord, y, x, w, theta, lambda);
  obj_list(0) = obj;
  int iter_best = 0;
  NumericVector std_y(y);

  for (int iter = 0; iter < M; iter++) {
    if (iter % 50 == 0)
      Rcpp::checkUserInterrupt();
    theta_old = theta;

    // define new(Gaussianized) data for least squares problem
    std_y = gaussianized_data(y, w, theta);
    // solve least squares problem (Gaussian TF)
    admm_gauss(Minner, n, ord, std_y, x, w, theta, z, u, rho, lam_z, DD, tol,
               inner_iter);
    total_iter += inner_iter;
    //  line search for step size
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
        iter_best = iter + 1;
        break;
      }
      iter_best = iter + 1;
    }

    // adjust the iterate steps
    if (iter >= iter_best + 4)
      break;
  }
}

/**
 * This is a wrapper around the void function to use in test_that()
 */
// [[Rcpp::export]]
Rcpp::List prox_newton_testing(int M,
                               int Minner,
                               int Mline,
                               int ord,
                               Rcpp::NumericVector const& y,
                               Rcpp::NumericVector const& x,
                               Rcpp::NumericVector const& w,
                               double lambda,
                               double ls_alpha,
                               double ls_gamma,
                               double tol) {
  Eigen::SparseMatrix<double> Dk;
  Eigen::SparseMatrix<double> DkDk;
  Dk = get_Dtil(ord, x);
  DkDk = Dk.transpose() * Dk;
  int m = Dk.rows();
  int n = y.size();
  NumericVector beta(n);
  NumericVector z(m);
  NumericVector u(m);
  double rho = lambda;
  int iter = 0;
  prox_newton(M, Minner, Mline, n, ord, y, x, w, beta, z, u, lambda, rho,
              ls_alpha, ls_gamma, DkDk, tol, iter);
  List out = List::create(Named("lambda") = lambda, Named("theta") = exp(beta),
                          Named("niter") = iter);
  return out;
}
