#include <RcppArmadillo.h>
#include <cmath>
#include "utils.h"

using namespace Rcpp;
using namespace arma;

/**
 * Generate a (banded) divided difference matrix of an arbitrary order for
 * equally spaced cases.
 * @param n an integer; signal length.
 * @param ord an integer; ord = -1, 0, 1, 2, ...
 * @return a sparse matrix D of order `ord + 1`. If `ord == -1`, return an
 * identity matrix.
 */
// [[Rcpp::export]]
arma::sp_mat buildD(int n, int ord) {
  if (ord + 1 < 0)
    stop("ord+1 must be non-negative.");
  if (n <= ord + 1)
    stop("n must be larger than ord + 1.");

  if (ord == -1) {
    sp_mat D = speye<sp_mat>(n, n);
    return D;
  }
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
 * @param ord an integer; ord = -1, 0, 1, ...
 * @param x a vector of signal locations. If no input, return equally spaced
 * matrix.
 * @return a sparse matrix Dmat. Unequally-spaced divided difference matrix of
 * order `ord+1`.
 */
// [[Rcpp::export]]
arma::sp_mat buildDx(int n, int ord, const arma::vec& x) {
  if (ord + 1 < 0)
    stop("ord + 1 must be non-negative.");
  if (n <= ord + 1)
    stop("n must be larger than ord + 1.");

  if (ord == -1 || x.size() == 0) {
    arma::sp_mat Dmat = buildD(n, ord);
    return Dmat;
  } else if (x.size() == n && !x.is_sorted("strictascend")) {
    stop("`x` must be in an increasing order with no replicated numbers ");
  } else if (x.size() != n)
    stop("x must be NULL or of length n.");

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
  if (ord < 0)
    stop("ord must be non-negative.");
  if (n <= ord + 1)
    stop("n must be larger than ord + 1.");

  if (x.size() == 0) {
    arma::sp_mat Dmat = buildD(n, ord - 1);
    return Dmat;
  } else if (x.size() == n && !x.is_sorted("strictascend")) {
    stop("`x` must be in an increasing order with no replicated numbers ");
  } else if (x.size() != n)
    stop("x must be NULL or of length n.");

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

/**
 * Compute lambda sequence
 * @param lambda hyperparameter
 * @param lambdamax max of lambda sequence
 * @param lambdamin min of lambda sequence
 * @param lambda_min_ratio ratio of lamdbamax/lambdamin
 * @param nsol length of lambda sequence
 */
void create_lambda(arma::vec& lambda,
                   double& lambdamin,
                   double& lambdamax,
                   double& lambda_min_ratio,
                   int& nsol) {
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

/**
 * A wrapper for unit testing
 */
// [[Rcpp::export()]]
arma::vec create_lambda_test(arma::vec lambda,
                             double lambdamin,
                             double lambdamax,
                             double lambda_min_ratio,
                             int nsol) {
  create_lambda(lambda, lambdamin, lambdamax, lambda_min_ratio, nsol);
  return lambda;
}

/**
 * Compute Gaussianized signals for proximal Newton algorithm
 * @param y observed signals
 * @param w weights
 * @param theta primal variable
 * @return gaussianized signals
 */
// [[Rcpp::export]]
arma::vec gaussianized_data(arma::vec const& y,
                            arma::vec const& w,
                            arma::vec& theta) {
  int n = y.size();
  vec c(n);
  for (int i = 0; i < n; i++) {
    if (w[i] * exp(theta[i]) > 1e-3)
      c[i] = y[i] * exp(-theta[i]) / w[i] - 1 + theta[i];
    else  // deal with overflow using approximation
      c[i] = y[i] - exp(theta[i]) / w[i] + theta[i];
  }
  return c;
}

/**
 * Compute Poisson loss
 * @param y observed signals
 * @param w signal weights
 * @param theta primal variable
 * @param lambda hyperparameter
 * @param Dv D*theta
 * @return Poisson loss
 */
double pois_obj(arma::vec const& y,
                arma::vec const& w,
                arma::vec& theta,
                double lambda,
                arma::vec& Dv) {
  vec v = -y % theta + w % exp(theta);
  double obj = mean(v) + lambda * norm(Dv, 1);
  return obj;
}

/**
 * Compute step size using backtracking line search for proximal Newton
 * algorithm
 * @param s step size to be returned
 * @param lambda hyperparameter
 * @param alpha scale in (0,0.5] adjusting upper bound. alpha=0.5 is
 * recommended.
 * @param gamma scale in (0,1) adjusting step size. gamma=0.9 is recommended.
 * @param y Gaussianized signals returned by `gaussianized_data`
 * @param x signal locations
 * @param w signal weights
 * @param n signal length
 * @param ord degree of trend filtering
 * @param theta primal variable updated by `admm_gauss`
 * @param theta_old primal variable before using `admm_gauss`
 * @param c1 a buffer of length `n-ord`
 * @param c2 a buffer of length `n-ord`
 * @param M maximum iteration
 * @return step size
 */
// [[Rcpp::export]]
double line_search(double s,
                   double lambda,
                   double alpha,
                   double gamma,
                   arma::vec const& y,
                   arma::vec const& x,
                   arma::vec const& w,
                   int n,
                   int ord,
                   arma::vec& theta,
                   arma::vec& theta_old,
                   arma::vec& c1,
                   arma::vec& c2,
                   int M) {
  calcDvline(n, ord, x, theta, c1);      // c1 = D * theta
  calcDvline(n, ord, x, theta_old, c2);  // c2 = D * theta_old
  double bound;
  double gradient;
  double grad_h;
  vec grad_g(n);
  vec Dv(n);
  vec dir = theta - theta_old;  // proximal Newton direction
  s = 1;                        // initialize step size

  for (int i = 0; i < M; i++) {
    // compute gradient/ grades
    grad_g = -s * dir % y + w % exp(theta_old + s * dir) - w % exp(theta_old);
    if (i > 0) {
      Dv = c1 - c2;
      Dv *= s;
      Dv += c2;  // if s=1, Dv stays same
    }
    grad_h = lambda * (norm(Dv, 1) - norm(c2, 1));
    gradient = mean(grad_g) + grad_h;

    // compute upper bound
    Dv = dir % (w % exp(theta_old) - y);
    bound = mean(Dv) * s;
    bound += grad_h;
    bound *= alpha;

    // check criteria
    if (gradient <= bound)
      break;
    else
      s *= gamma;
  }
  return s;
}

/**
 * Calculate b = D * v in place
 * @param n signal length
 * @param ord degree of the Poisson trend filtering problem
 * @param x signal locations
 * @param v vector of size n
 * @param b vector of size (n - ord); a buffer storing results
 */
void calcDvline(int n,
                int ord,
                arma::vec const& x,
                arma::vec& v,
                arma::vec& b) {
  b = v;
  int fct = 1;
  for (int i = 0; i < ord; i++) {
    if (i != 0 && x.size() > 0)
      b /= (x.tail(n - i) - x.head(n - i));
    b = b.tail(n - i - 1) - b.head(n - i - 1);
    b.resize(n - i - 1);
  }
  for (int i = 2; i < ord; i++)
    fct *= i;
  b *= fct;
}

/**
 * A wrapper for testing `calcDvline`
 */
// [[Rcpp::export]]
arma::vec calcDvline_test(int n,
                          int ord,
                          arma::vec const& x,
                          arma::vec& v,
                          arma::vec& b) {
  arma::sp_mat D = buildDx_tilde(n, ord, x);
  b = D * v;
  return b;
}

/**
 * Calculate b = D^T * v in place
 * @param v vec of length `n-ord`
 * @param b vec of length n; a buffer storing results
 */
void calcDTvline(int n,
                 int ord,
                 arma::vec const& x,
                 arma::vec& v,
                 arma::vec& b) {
  b.head(n - ord) = v;
  int fct = 1;

  for (int i = ord; i > 0; i--) {
    b(n - i) = b(n - i - 1);
    for (int j = n - i - 1; j > 0; j--)
      b(j) = b(j - 1) - b(j);
    b[0] = -b[0];
    if (i != 1 && x.size() > 0) {
      b.head(n - i + 1) /= (x.tail(n - i + 1) - x.head(n - i + 1));
    }
  }
  for (int i = 2; i < ord; i++)
    fct *= i;
  b *= fct;
}

/*
 * A wrapper for testing `calcDTvline`
 */
// [[Rcpp::export]]
arma::vec calcDTvline_test(int n,
                           int ord,
                           arma::vec const& x,
                           arma::vec& v,
                           arma::vec& b) {
  arma::sp_mat D = buildDx_tilde(n, ord, x);
  b = D.t() * v;
  return b;
}

/**
 * Calculate b = D^T * D * v in place
 * @param v vec of length n
 * @param b vec of length n; a buffer storing results
 */
void calcDTDvline(int n,
                  int ord,
                  arma::vec const& x,
                  arma::vec& v,
                  arma::vec& b) {
  vec c = v;
  int fct = 1;
  for (int i = 0; i < ord; i++) {
    if (i != 0 && x.size() > 0)
      c /= (x.tail(n - i) - x.head(n - i));
    c = c.tail(n - i - 1) - c.head(n - i - 1);
    c.resize(n - i - 1);
  }
  b.head(n - ord) = c;
  for (int i = ord; i > 0; i--) {
    b(n - i) = b(n - i - 1);
    for (int j = n - i - 1; j > 0; j--)
      b(j) = b(j - 1) - b(j);
    b[0] = -b[0];
    if (i != 1 && x.size() > 0) {
      b.head(n - i + 1) /= (x.tail(n - i + 1) - x.head(n - i + 1));
    }
  }
  for (int i = 2; i < ord; i++)
    fct *= i * i;
  b *= fct;
}

/*
 * A wrapper for testing `calcDTDvline`
 */
// [[Rcpp::export]]
arma::vec calcDTDvline_test(int n,
                            int ord,
                            arma::vec const& x,
                            arma::vec& v,
                            arma::vec& b) {
  arma::sp_mat D = buildDx_tilde(n, ord, x);
  b = D.t() * D * v;
  return b;
}
