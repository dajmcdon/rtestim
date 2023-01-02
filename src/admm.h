#ifndef __ADMM_H
#define __ADMM_H

#include <Rcpp.h>

double update_pois(double c, double mu, int n);

void admm(int M, arma::vec const& y, arma::vec const& w, int n,
          arma::vec& theta, arma::vec& z, arma::vec& u,
          double lambda, double rho, double mu,
          arma::sp_mat const& DD, arma::sp_mat const& D,
          double tol);

Rcpp::List admm_testing(int M,
                  arma::vec const& y,
                  arma::vec const& w,
                  int n,
                  arma::vec theta,
                  arma::vec z,
                  arma::vec u,
                  double lambda,
                  double rho,
                  double mu,
                  arma::sp_mat const& DD,
                  arma::sp_mat const& D,
                  double tol);


#endif
