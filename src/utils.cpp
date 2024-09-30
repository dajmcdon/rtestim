#include <Eigen/Sparse>
#include <RcppEigen.h>
#include <dspline.h>
#include "utils.h"

using namespace Rcpp;
using Eigen::VectorXd;
using Eigen::MatrixXd;

NumericVector evec_to_nvec(VectorXd evec) {
  NumericVector nvec(wrap(evec));
  return nvec;
}

VectorXd nvec_to_evec(NumericVector nvec) {
  VectorXd evec = as<Eigen::Map<VectorXd> >(nvec);
  return evec;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> get_Dtil(int k, NumericVector xd) {
  int n = xd.size();
  return dspline::rcpp_b_mat(k, xd, false, Rcpp::seq(0, n - k - 1), true);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> get_D(int k, NumericVector xd) {
  k++;
  int n = xd.size();
  return dspline::rcpp_b_mat(k, xd, true, Rcpp::seq(0, n - k - 1), true);
}

void create_lambda(NumericVector& lambda,
                   double& lambdamin,
                   double& lambdamax,
                   double& lambda_min_ratio,
                   int& nsol) {
  if (all(lambda == 0).is_false()) {
    lambdamin = min(lambda);
    lambdamax = max(lambda);
    nsol = lambda.size();
  } else {
    lambdamin = (lambdamin < 0) ? lambda_min_ratio * lambdamax : lambdamin;
    double lmpad = 1e-20;
    double ns = static_cast<double>(nsol);
    double p;
    if (lambdamin > lmpad) {
      p = pow(lambdamin / lambdamax, 1 / (ns - 1));
      lambda(0) = lambdamax;
      for (int i = 1; i < nsol; i++) lambda[i] = lambda[i - 1] * p;
    } else {
      p = pow(lmpad / lambdamax, 1 / (ns - 2));
      lambda(0) = lambdamax;
      for (int i = 1; i < nsol - 1; i++) lambda[i] = lambda[i - 1] * p;
      lambda(nsol - 1) = lambdamin;
    }
  }
}

// [[Rcpp::export]]
NumericVector create_lambda_test(NumericVector lambda,
                                 double lambdamin,
                                 double lambdamax,
                                 double lambda_min_ratio,
                                 int nsol) {
  create_lambda(lambda, lambdamin, lambdamax, lambda_min_ratio, nsol);
  return lambda;
}

// [[Rcpp::export]]
NumericVector doDv(NumericVector v, int k, NumericVector xd) {
  return dspline::rcpp_d_mat_mult(v, k, xd, false, false);
}

// [[Rcpp::export]]
NumericVector doDtv(NumericVector v, int k, NumericVector xd) {
  return dspline::rcpp_d_mat_mult(v, k, xd, false, true);
}

// [[Rcpp::export]]
double one_norm(NumericVector const& z) {
  double out = sum(abs(z));
  return out;
}

// [[Rcpp::export]]
double pois_obj(int ord,
                NumericVector const& y,
                NumericVector const& x,
                NumericVector const& w,
                NumericVector& theta,
                double lambda) {
  NumericVector Dtheta = doDv(theta, ord, x);
  NumericVector v = -y * theta + w * exp(theta);
  double obj = mean(v) + lambda * sum(abs(Dtheta));
  return obj;
}

// [[Rcpp::export]]
NumericVector centered_data(NumericVector const& y, NumericVector const& w,
                            NumericVector& theta) {
  int n = y.size();
  NumericVector out = clone(theta);
  for (int i = 0; i < n; i++) {
    if (w(i) * exp(theta(i)) > 1e-3) {
      out(i) += y(i) * exp(-theta(i)) / w(i) - 1;
    } else if (w(i) < 1e-6) {
      out(i) += y(i) - exp(theta(i)) / 1e-6;
    } else { // deal with overflow using approximation
      out(i) += y(i) - exp(theta(i)) / w(i);
    }
  }
  return out;
}

/**
 * Backtracking line search for proximal Newton method
 * @param s step size
 * @param alpha
 * @param gamma contraction factor
 */
// [[Rcpp::export]]
double line_search(double s,
                   double lambda,
                   double alpha,
                   double gamma,
                   NumericVector const& y,
                   NumericVector const& x,
                   NumericVector const& w,
                   int n,
                   int ord,
                   NumericVector& theta,
                   NumericVector& theta_old,
                   int M) {
  NumericVector dir = theta - theta_old;
  double bound;
  double gradient;
  double bound_s = 0.0;
  NumericVector grad_g(n);
  NumericVector Dv(n);

  // initialize upper bound
  bound = mean(dir * (w * exp(theta) - y));
  NumericVector Dth = doDv(theta, ord, x);
  NumericVector Dth_old = doDv(theta_old, ord, x);
  double Dth_old_norm = one_norm(Dth_old);
  bound += lambda * (one_norm(Dth) - Dth_old_norm);
  NumericVector exp_theta_old = exp(theta_old);

  s = 1;  // initial step length
  for (int i = 0; i < M; i++) {
    // compute gradient: f(theta_old + s * dir) - f(theta_old)
    grad_g = -s * dir * y + w * exp(theta_old + s * dir) - w * exp_theta_old;
    Dv = Dth_old + s * (Dth - Dth_old);  // if s=1, Dv stays same
    gradient = mean(grad_g) + lambda * (one_norm(Dv) - Dth_old_norm);

    // adjust upper bound: alpha * s * (f'(theta_old))^T * dir
    bound_s = bound * alpha * s;

    // check the Armijo (sufficient decrease) condition
    if (gradient <= bound_s)
      break;
    else
      s *= gamma;
  }
  return s;
}



// [[Rcpp::depends(BH)]]

#include <boost/integer/common_factor.hpp>

// [[Rcpp::export]]
int compute_gcd(IntegerVector x) {
  std::sort(x.begin(), x.end());
  int g = 0;
  int n = x.length();
  if (n == 0) {
    g = 0;
  } else if (n == 1) {
    g = x[0];
  } else if (n == 2) {
    g = boost::integer::gcd(x[0], x[1]);
  } else {
    g = boost::integer::gcd(x[0], x[1]);
    for (int i = 2; i < n; i++) {
      g = boost::integer::gcd(g, x[i]);
      if (g == 1) break;
    }
  }
  return g;
}

// [[Rcpp::export]]
NumericVector calc_delays(NumericVector x, NumericVector y) {
  int n = x.size();
  double s = 0.0;
  NumericVector out(n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      out[i] += x[i-j] * y[j];
    }
    s += y[i];
    if (s > 1e-16) out[i] /= s;
  }
  return out;
}


// in-place computation of Ptemp given P and A
Eigen::MatrixXd computePtemp(Eigen::MatrixXd A, Eigen::MatrixXd P) {
  int k = P.rows();
  Eigen::MatrixXd temp = A.row(0) * P;
  Eigen::MatrixXd var = A.row(0) * temp.transpose();

  // in-place block replacement starting from the last component:
  P.block(1, 1, k - 1, k - 1).reverse() = P.block(0, 0, k - 1, k - 1).reverse();
  P(0, 0) = var(0, 0);
  temp.conservativeResize(1, k - 1);  // drop the last component
  P.block(0, 1, 1, k - 1) = temp;
  P.block(1, 0, k - 1, 1) = temp.transpose();
  return P;
}

// [[Rcpp::export]]
Eigen::MatrixXd smat_to_mat(const Eigen::SparseMatrix<double>& sparseMat, int k, bool equal_spaced) {
  int rows = sparseMat.rows(); // n-k
  Eigen::MatrixXd denseMat(rows, k + 1);
  std::vector<double> rowNonzeros;
  // Iterate over nonzero coefficients in each row of the sparse matrix
  for (int i = 0; i < sparseMat.outerSize(); ++i) {
    std::vector<double> rowNonzeros;
    for (Eigen::SparseMatrix<double>::InnerIterator it(sparseMat, i); it; ++it) 
      rowNonzeros.push_back(it.value());
    int m = rowNonzeros.size();
    for (int j = 0; j < m; j++) {
      if (i < rows) 
        denseMat(i - j, j) = rowNonzeros[m - 1 - j];
      else  
        denseMat(rows - 1 - j, j + i - rows + 1) = rowNonzeros[m - 1 - j];
    }
    if (equal_spaced && i == k) {
      denseMat.conservativeResize(1, k + 1);
      return denseMat;
    } 
  }
  return denseMat;
}
