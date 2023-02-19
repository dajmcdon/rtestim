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

#endif
