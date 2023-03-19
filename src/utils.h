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
arma::vec create_lambda_test(arma::vec lambda,
                             double lambdamin = -1,
                             double lambdamax = -1,
                             double lambda_min_ratio = 1e-4,
                             int nsol = 50);
arma::vec gaussianized_data(arma::vec const& y,
                            arma::vec const& w,
                            arma::vec& theta);
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
                   arma::vec const& x,
                   arma::vec const& w,
                   int n,
                   int ord,
                   arma::vec& theta,
                   arma::vec& theta_old,
                   arma::vec& c1,
                   arma::vec& c2,
                   int M);
void calcDvline(int n, int ord, arma::vec const& x, arma::vec& v, arma::vec& b);
arma::vec calcDvline_test(int n,
                          int ord,
                          arma::vec const& x,
                          arma::vec& v,
                          arma::vec& b);
void calcDTvline(int n,
                 int ord,
                 arma::vec const& x,
                 arma::vec& v,
                 arma::vec& b);
arma::vec calcDTvline_test(int n,
                           int ord,
                           arma::vec const& x,
                           arma::vec& v,
                           arma::vec& b);
void calcDTDvline(int n,
                  int ord,
                  arma::vec const& x,
                  arma::vec& v,
                  arma::vec& b);
arma::vec calcDTDvline_test(int n,
                            int ord,
                            arma::vec const& x,
                            arma::vec& v,
                            arma::vec& b);

#endif
