#include <RcppEigen.h>
#include <Eigen/Sparse>
#include <boost/math/special_functions/lambert_w.hpp>
#include "dptfe.h"
#include "utils-eigen.h"

typedef Eigen::COLAMDOrdering<int> Ord;

using Eigen::SparseMatrix;
using Eigen::SparseQR;
using Eigen::VectorXd;
SparseQR<SparseMatrix<double>, Ord> qr;



// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;

NumericVector evec_to_nvec(VectorXd evec) {
  NumericVector nvec(wrap(evec));
  return nvec;
}

VectorXd nvec_to_evec(NumericVector nvec) {
  VectorXd evec = as<Eigen::Map<VectorXd> >(nvec);
  return evec;
}


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

// [[Rcpp::export]]
NumericVector admm_gauss(int M,
                         int n,
                         int ord,
                         NumericVector const& y,
                         NumericVector const& x,
                         NumericVector const& w,
                         NumericVector& theta,
                         NumericVector& z,
                         NumericVector& u,
                         double rho,
                         double lam_z,
                         double r_norm,
                         double s_norm,
                         SparseMatrix<double> const& DD,
                         double tol) {

  SparseMatrix<double> cDD = DD * n * rho; // a copy that doesn't change

  NumericVector z_old = clone(z);
  NumericVector tmp_n(n);
  NumericVector tmp_m(z.size());
  NumericVector Dth(z.size());
  VectorXd tmp_theta(n);

  for (int iter = 0; iter < M; iter++) {
    if (iter % 1000 == 0) Rcpp::checkUserInterrupt();
    // solve for primal variable - theta:
    SparseMatrix<double> wDD(cDD); // another copy that we alter
    SparseMatrix<double> W = identity(n);
    VectorXd etheta = nvec_to_evec(theta);
    W = W * etheta.exp();
    tmp_n = doDtv(z - u, ord, x) * n * rho;
    tmp_n += exp(theta) * y;
    tmp_theta = nvec_to_evec(tmp_n);
    qr.compute(W);
    tmp_theta = qr.solve(tmp_theta);
    theta = evec_to_nvec(tmp_theta);
    // solve for alternating variable - z:
    Dth = doDv(theta, ord, x);
    tmp_m = Dth + u;
    z = dptfe(tmp_m, lam_z);
    // update dual variable - u:
    u = Dth - z;

    // primal residuals:
    r_norm = sqrt(mean(pow(Dth - z, 2)));
    // dual residuals:
    tmp_n = doDtv(z - z_old, ord, x);
    s_norm = rho * sqrt(mean(pow(tmp_n, 2)));
    // stopping criteria check:
    if (r_norm < tol && s_norm < tol) break;

    // auxiliary variables update:
    z_old = z;
  }
  return theta;
}

void irls_admm(int M,
               int n,
               int ord,
               NumericVector const& y,
               NumericVector const& x,
               NumericVector const& w,
               NumericVector& theta,
               NumericVector& z,
               NumericVector& u,
               double lambda,
               double rho,
               double mu,
               double alpha,
               double gamma,
               SparseMatrix<double> const& DD,
               double tol,
               int& iter) {
  double s;             // step size
  NumericVector obj_list(M + 1);  // objective list for each iterate
  double obj = 1e4;     // initialize it to be large
  NumericVector theta_old(n);     // a buffer for line search
  double lam_z = lambda / rho;
  int m = z.size();
  double r_norm = 0.0;
  double s_norm = 0.0;

  // initialize objective
  obj = pois_obj(ord, y, x, w, theta, lambda);
  obj_list(0) = obj;
  int iter_best = 0;

  for (iter = 0; iter < M; iter++) {
    if (iter % 500 == 0) Rcpp::checkUserInterrupt();
    theta_old = theta;

    // define new(Gaussianized) data for least squares problem
    NumericVector std_y = gaussianized_data(y, w, theta);
    // solve least squares problem (Gaussian TF)
    theta = admm_gauss(M, n, ord, std_y, x, w, theta, z, u, rho, lam_z, r_norm,
                       s_norm, DD, tol);

    // line search for step size
    s = line_search(s, lambda, alpha, gamma,
                    y, x, w, n, ord, theta, theta_old, M);
    if (s < 0) break;
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
    if (iter >= iter_best + 4 && iter_best != 0)  break;
  }
}

