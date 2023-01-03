#include <RcppArmadillo.h>
#include "admm.h"
#include "utils.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
List rtestim_path(
    arma::vec y,
    arma::vec x, // positions
    arma::vec w, // weighted past cases
    int korder,
    arma::vec lambda,
    double lambdamax = -1,
    double lambdamin = -1,
    int nsol = 100,
    double rho_adjust = -1,
    double rho = -1,
    int maxiter = 1e5,
    double tolerance = 1e-3,
    double lambda_min_ratio = 1e-4,
    int verbose = 0) {

  int n = y.n_elem;
  if (lambda.size() > 0) nsol = lambda.size();

  // Placeholders for solutions
  arma::mat theta(n, nsol);
  arma::vec niter(nsol);

  // Build D matrix
  arma::sp_mat D;
  arma::sp_mat Dk;
  D = buildDx(n, korder, x);
  Dk = buildDx_tilde(n, korder, x);
  arma::sp_mat DkDk = Dk.t() * Dk;

  // Generate lambda sequence if necessary
  if (lambda.size() == 0 && lambdamax <= 0) {
    arma::vec b(n - korder );
    arma::spsolve(b, D.t(), w % y); // very approximate
    lambdamax = arma::norm(b, "inf");
    lambdamax *= n;
  }
  create_lambda(lambda, lambdamin, lambdamax, lambda_min_ratio, nsol);

  // ADMM parameters
  double tolerance_abs = std::sqrt(n) * tolerance; // adjust for number of parameters
  double _rho = (rho < 0) ? lambda(0) : rho;


  // ADMM variables
  arma::vec beta(n, arma::fill::zeros);
  arma::vec alpha(D.n_rows, arma::fill::zeros);
  arma::vec u(D.n_rows, arma::fill::zeros);
  double mu = 2 * pow(4, korder);

  // Outer loop to compute solution path
  for (int i = 0; i < nsol; i++) {
    if (i > 0 && rho < 0) _rho = lambda(i);
    if (verbose > 0) Rcout << ".";
    Rcpp::checkUserInterrupt();
    int iters = 0;



    admm(maxiter, y, w, n, beta, alpha, u, lambda(i),
         _rho, mu * lambda(i), DkDk, D, tolerance_abs);

      // old version
      // DkDk, Dk,
      // y, w,
      // beta, alpha, u,
      // lambda(i), _rho,
      // rho_adjust,
      // maxiter, tolerance_abs, iters);

    // Store solution
    theta.col(i) = exp(beta);
    niter(i) = iters;

    // Verbose handlers
    if (verbose > 1) Rcout << niter(i);
    if (verbose > 2) Rcout << "(" << lambda(i) << ")";
    if (verbose > 0) Rcout << std::endl;
  }



  // Return
  List out = List::create(
    Named("y") = y,
    Named("n") = n,
    Named("lambda") = lambda,
    Named("theta") = theta,
    Named("niter") = niter
  );
  return out;
}
