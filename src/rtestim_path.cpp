#include <RcppArmadillo.h>
#include "dptf.h"
#include "admm.h"
#include "utils.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
List rtestim_path(int algo,
                  arma::vec y,
                  arma::vec x,  // positions
                  arma::vec w,  // weighted past cases
                  int korder,
                  arma::vec lambda,
                  double lambdamax = -1,
                  double lambdamin = -1,
                  int nsol = 100,
                  double rho = -1,
                  int maxiter = 1e5,
                  double tolerance = 1e-3,
                  double lambda_min_ratio = 1e-4,
                  double ls_alpha = 0.5,
                  double ls_gamma = 0.9,
                  int verbose = 0) {
  int n = y.n_elem;
  if (lambda.size() > 0)
    nsol = lambda.size();

  // Placeholders for solutions
  arma::mat theta(n, nsol);
  arma::vec niter(nsol);
  arma::vec dof(nsol);

  // Build D matrix
  arma::sp_mat D;
  arma::sp_mat Dk;
  D = buildDx(n, korder - 1, x);
  Dk = buildDx_tilde(n, korder, x);
  arma::sp_mat DkDk = Dk.t() * Dk;

  // Generate lambda sequence if necessary
  if (lambda.size() == 0 && lambdamax <= 0) {
    arma::vec b(n - korder);
    arma::mat matD(D);                // convert to dense mat to avoid spsolve
    arma::solve(b, matD.t(), w % y);  // very approximate; for unevenly space?
    lambdamax = arma::norm(b, "inf");
    lambdamax /= n;
  }
  create_lambda(lambda, lambdamin, lambdamax, lambda_min_ratio, nsol);

  // ADMM parameters
  double _rho = (rho < 0) ? lambda(0) : rho;

  // ADMM variables
  arma::vec beta(n, arma::fill::zeros);
  arma::vec alpha(Dk.n_rows, arma::fill::zeros);
  arma::vec u(Dk.n_rows, arma::fill::zeros);
  double mu = 2 * pow(4, korder);  // unevenly-spaced version?
  int iters = 0;

  // Outer loop to compute solution path
  for (int i = 0; i < nsol; i++) {
    if (verbose > 0)
      Rcout << ".";
    Rcpp::checkUserInterrupt();

    if (korder == 0) {
      beta = dptf_past(y, lambda(i), w);
      niter(i) = 1;
    } else {
      if (i > 0)
        _rho = (rho < 0) ? lambda(i) : rho;
      switch (algo) {
        case 1:
          admm(maxiter, y, x, w, n, korder, beta, alpha, u, lambda(i), _rho,
               mu * lambda(i), tolerance, iters);  // add rho_adjust?
          break;
        case 2:
          irls_admm(maxiter, n, korder, y, x, w, beta, alpha, u, lambda(i),
                    _rho, mu * lambda(i), ls_alpha, ls_gamma, Dk, tolerance,
                    iters);
          break;
      }
      niter(i) = iters + 1;
    }

    // Store solution
    if (korder == int(0)) {
      theta.col(i) = beta;  // the algorithm returns mean (exp of the
                            // natural parameter)
      dof(i) = arma::sum(arma::abs(arma::diff(beta)) > tolerance);
    } else {
      theta.col(i) = exp(beta);
      dof(i) = arma::sum(arma::abs(arma::diff(alpha)) > tolerance);
    }

    // Verbose handlers
    if (verbose > 1)
      Rcout << niter(i);
    if (verbose > 2)
      Rcout << "(" << lambda(i) << ")";
    if (verbose > 0)
      Rcout << std::endl;
  }

  // Return
  List out = List::create(Named("Rt") = theta, Named("lambda") = lambda,
                          Named("degree") = korder + 1, Named("dof") = dof,
                          Named("niter") = niter);
  return out;
}
