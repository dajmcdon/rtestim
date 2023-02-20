#ifndef __ADMM_EIGEN_H
#define __ADMM_EIGEN_H


double update_epois(double c, double mu, int n);

void admm_eigen(int M,
                Rcpp::NumericVector const& y,
                Rcpp::NumericVector const& x,
                Rcpp::NumericVector const& w,
                int n,
                int ord,
                Rcpp::NumericVector& theta,
                Rcpp::NumericVector& z,
                Rcpp::NumericVector& u,
                double lambda,
                double rho,
                double mu,
                double tol,
                int& iter);


Rcpp::List admm_eigen_testing(int M,
                              Rcpp::NumericVector const& y,
                              Rcpp::NumericVector const& x,
                              Rcpp::NumericVector const& w,
                              int n,
                              int ord,
                              Rcpp::NumericVector theta,
                              Rcpp::NumericVector z,
                              Rcpp::NumericVector u,
                              double lambda,
                              double rho,
                              double mu,
                              double tol,
                              int iter);

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
                               double tol);

void irls_admm(int M,
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
               double mu,
               double alpha,
               double gamma,
               Eigen::SparseMatrix<double> const& DD,
               double tol,
               int& iter)

#endif
