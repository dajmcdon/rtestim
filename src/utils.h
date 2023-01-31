#ifndef __UTILS_H
#define __UTILS_H

arma::sp_mat buildD(int n, int ord);
arma::sp_mat buildDx(int n, int k, const arma::vec& x);
arma::sp_mat buildDx_tilde(int n, int k, const arma::vec& x);
void create_lambda(arma::vec& lambda,
                   double& lambdamin,
                   double& lambdamax,
                   double& lambda_min_ratio,
                   int& nsol);
arma::vec fake_data(arma::vec const& y, arma::vec const& w, arma::vec& theta);
double pois_obj(arma::vec const& y,
                arma::vec const& w,
                arma::vec& theta,
                double lambda,
                arma::vec& Dv);
double line_search(double s,
                   double lambda,
                   double alpha,
                   double gamma,
                   arma::vec const& y,
                   arma::vec const& w,
                   int n,
                   arma::vec& theta,
                   arma::vec& theta_old,
                   arma::vec& c1,
                   arma::vec& c2,
                   arma::sp_mat const& D,
                   int M);

#endif
