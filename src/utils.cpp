#include <RcppArmadillo.h>
#include <cmath>
#include "utils.h"

using namespace Rcpp;

/**
 * Generate a (banded) divided difference matrix of an arbitrary order for
 * equally spaced cases.
 * @param n an integer; signal length.
 * @param ord an integer; ord = -1, 0, 1, ... Order of the divided difference
 * matrix. If `ord == -1`, return an identity matrix (of `ord + 2 = 1` band).
 * @return a sparse matrix D
 */
// [[Rcpp::export]]
arma::sp_mat buildD(int n, int ord) {
  // stop if not: n > ord + 1; ord + 2 > 0 (; n > 0)
  int c1 = ord + 1;
  int m = n - c1;
  arma::sp_mat D(m, n);
  double z;
  for (int i = 0; i <= c1; i++) {
    z = Rf_choose(c1, i) * std::pow(-1, i + c1);
    D.diag(i) += z;
  }
  return D;
}

/**
 * Generate a (banded) divided difference matrix of an arbitrary order for
 * unequally spaced cases.
 * @param n an integer; signal length
 * @param ord an integer; ord = 0, 1, ... Order of the divided difference
 * matrix.
 * @param x a vector of signal locations. If no input, return equally spaced
 * matrix.
 * @return a sparse matrix Dmat
 */
// [[Rcpp::export]]
arma::sp_mat buildDx(int n, int ord, const arma::vec& x) {
  // stop if not: n > ord + 1; ord >= 0 (; n > 1)
  if (x.size() == 0) {
    arma::sp_mat Dmat = buildD(n, ord);
    return Dmat;
  }
  // stop if not: x.size() == n.size()
  // x should be in an increasing order (with no replicated numbers)
  arma::sp_mat D1(n - 1, n);
  D1.diag(0) -= 1;
  D1.diag(1) += 1;
  if (ord == 0)
    return D1;                      // ord = 0 is the same as usual
  arma::sp_mat Dmat = D1;           // output
  arma::sp_mat delx(n - 1, n - 1);  // diagonal matrix adjusting locations

  for (int tk = 1; tk <= ord; tk++) {  // use the way in Eq.12, Sec4, TF14supp
    D1.resize(n - tk - 1, n - tk);     // D^(1) with shrinking dim
    delx.diag(0) = tk / (x.tail(n - tk) - x.head(n - tk));
    Dmat = D1 * delx * Dmat;
    delx.set_size(n - tk - 1, n - tk - 1);
  }
  return Dmat;
}

/**
 * Generate a (banded) divided difference matrix D of an arbitrary order for
 * unequally spaced cases.
 * @param n an integer; signal length
 * @param ord an integer; ord = 0, 1, ... Order of the divided difference matrix
 * minus 1.
 * @param x a vector of signal locations. If no input, return equally spaced
 * matrix.
 * @return a sparse matrix Dmat/Dk
 */
// [[Rcpp::export]]
arma::sp_mat buildDx_tilde(int n, int ord, const arma::vec& x) {
  // stop if not: n > ord + 1; ord >= 0 (; n > 1)
  if (x.size() == 0) {
    arma::sp_mat Dmat = buildD(n, ord - 1);
    return Dmat;
  }
  // stop if not: x.size() == n.size()
  // x should be in an increasing order (with no replicated numbers)
  if (ord == 0) {
    arma::sp_mat Dk(n, n);
    Dk.eye();
    return Dk;
  }
  arma::sp_mat Dk = buildDx(n, ord - 1, x);
  arma::vec lagdiff = ord / (x.tail(n - ord) - x.head(n - ord));
  arma::sp_mat delx(n - ord, n - ord);
  delx.diag(0) = lagdiff;
  Dk = delx * Dk;
  return Dk;
}

void create_lambda(arma::vec& lambda,
                   double lambdamin,
                   double lambdamax,
                   double lambda_min_ratio,
                   int nsol) {
  if (lambda.size() > 0) {
    lambdamin = lambda.min();
    lambdamax = lambda.max();
    nsol = lambda.size();
  } else {
    lambda.set_size(nsol);
    lambdamin = (lambdamin < 0) ? lambda_min_ratio * lambdamax : lambdamin;
    double lmpad = 1e-20;
    if (lambdamin < lmpad) {
      lambda.tail(nsol - 1) =
          arma::logspace(log10(lmpad), log10(lambdamax), nsol - 1);
      lambda(0) = lambdamin;
    } else {
      lambda = arma::logspace(log10(lambdamin), log10(lambdamax), nsol);
    }
  }
}
