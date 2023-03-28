#include <Eigen/Sparse>
#include <RcppEigen.h>
#include "dptf.h"
#include "admm.h"
#include "utils.h"

typedef Eigen::COLAMDOrdering<int> Ord;

using Eigen::SparseMatrix;
using Eigen::SparseQR;
using Eigen::VectorXd;
SparseQR<SparseMatrix<double>, Ord> qr;

// [[Rcpp::depends(RcppEigen)]]


using namespace Rcpp;

// [[Rcpp::export]]
List rtestim_path(int algo,
                  NumericVector y,
                  NumericVector x,  // positions
                  NumericVector w,  // weighted past cases
                  int korder,
                  NumericVector lambda,
                  double lambdamax = -1,
                  double lambdamin = -1,
                  int nsol = 100,
                  double rho = -1,
                  int maxiter = 1e5,
                  double tolerance = 1e-3,
                  double lambda_min_ratio = 1e-4,
                  double ls_alpha = 0.5,
                  double ls_gamma = 0.9,
                  int maxiter_inner = 30,
                  int verbose = 0) {
  int n = y.size();

  // Placeholders for solutions
  NumericMatrix theta(n, nsol);
  NumericVector niter(nsol);
  NumericVector dof(nsol);

  // Build D matrices as needed
  Eigen::SparseMatrix<double> D;
  Eigen::SparseMatrix<double> Dk;
  D = get_D(korder, x);
  qr.compute(D.transpose());
  int m = n;
  if (korder > 0) {
    Dk = get_Dtil(korder, x);
    Eigen::SparseMatrix<double> DkDk = Dk.transpose() * Dk;
    m = Dk.rows();
  }
  // Generate lambda sequence if necessary
  if (abs(lambda[n-1]) < tolerance / 100 && lambdamax <= 0) {
    VectorXd b(n - korder);
    VectorXd wy = nvec_to_evec(w * y);
    b = qr.solve(wy); // very approximate;
    NumericVector bp = evec_to_nvec(b);
    lambdamax = max(abs(bp)) / n;
  }
  create_lambda(lambda, lambdamin, lambdamax, lambda_min_ratio, nsol);

  // ADMM parameters
  double _rho;
  double _mu;

  // ADMM variables
  NumericVector beta(n);
  NumericVector alpha(m);
  NumericVector u(m);
  double mu = 2 * pow(4, korder);  // unevenly-spaced version?
  int iters = 0;
  int nsols = nsol;

  // Outer loop to compute solution path
  for (int i = 0; i < nsol; i++) {
    if (verbose > 0) Rcout << ".";
    Rcpp::checkUserInterrupt();

    if (korder == 0) {
      beta = dptf_past(y, lambda[i], w);
      niter[i] = 1;
    } else {
      _rho = (rho < 0) ? lambda[i] : rho;
      _mu = mu * lambda[i];
      switch (algo) {
        case 1:
          admm(maxiter, y, x, w, n, korder, beta, alpha, u, lambda[i], _rho,
               _mu, tolerance, iters);  // add rho_adjust?
          break;
        case 2:
          irls_admm(maxiter, n, korder, y, x, w, beta, alpha, u, lambda[i],
                    _rho, ls_alpha, ls_gamma, Dk, tolerance,
                    maxiter_inner, iters);
          break;
      }
      niter[i] = iters + 1;
      maxiter -= iters;
      if (maxiter < 0) nsols = i + 1;
    }


    // Store solution
    if (korder == int(0)) {
      theta(_, i) = beta;
      dof[i] = sum(abs(diff(beta)) > tolerance);
    } else {
      theta(_, i) = exp(beta);
      dof[i] = sum(abs(diff(alpha)) > tolerance);
    }


    // Verbose handlers
    if (verbose > 1) Rcout << niter(i);
    if (verbose > 2) Rcout << "(" << lambda(i) << ")";
    if (verbose > 0) Rcout << std::endl;
    if (maxiter < 0) break;
  }


  // Return
  List out = List::create(
    Named("Rt") = theta(_, Range(0, nsols-1)),
    Named("lambda") = lambda[Range(0, nsols-1)],
    Named("degree") = korder + 1,
    Named("dof") = dof[Range(0, nsols - 1)],
    Named("niter") = niter[Range(0, nsols-1)]
  );
  return out;
}
