#ifndef __ADMM_H
#define __ADMM_H

double update_primal(double c, double mu2);

Rcpp::List admm(int M,
                arma::vec y,
                arma::vec x,
                int n,
                arma::vec theta,
                arma::vec z,
                arma::vec u,
                double lambda,
                double rho,
                double mu,
                arma::sp_mat D,
                double tol);

#endif