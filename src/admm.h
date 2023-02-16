#ifndef __ADMM_H
#define __ADMM_H

double update_pois(double c, double mu, int n);

void admm(int M,
          arma::vec const& y,
          arma::vec const& x,
          arma::vec const& w,
          int n,
          int ord,
          arma::vec& theta,
          arma::vec& z,
          arma::vec& u,
          double lambda,
          double rho,
          double mu,
          double tol,
          int& iter);

Rcpp::List admm_testing(int M,
                        arma::vec const& y,
                        arma::vec const& x,
                        arma::vec const& w,
                        int n,
                        int ord,
                        arma::vec& theta,
                        arma::vec& z,
                        arma::vec& u,
                        double lambda,
                        double rho,
                        double mu,
                        double tol,
                        int iter);

arma::vec admm_gauss(int M,
                     int n,
                     int ord,
                     arma::vec const& y,
                     arma::vec const& x,
                     arma::vec const& w,
                     arma::vec& theta,
                     arma::vec& z,
                     arma::vec& u,
                     double rho,
                     double lam_z,
                     double r_norm,
                     double s_norm,
                     arma::sp_mat const& DD,
                     double tol);
void irls_admm(int M,
               int n,
               int ord,
               arma::vec const& y,
               arma::vec const& x,
               arma::vec const& w,
               arma::vec& theta,
               arma::vec& z,
               arma::vec& u,
               double lambda,
               double rho,
               double mu,
               double alpha,
               double gamma,
               arma::sp_mat const& D,
               double tol,
               int& iter);

#endif
