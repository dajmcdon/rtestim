#ifndef __UTILS_H
#define __UTILS_H

arma::sp_mat buildD(int n, int ord);
arma::sp_mat buildDx(int n, int k, const arma::vec& x);
arma::sp_mat buildDx_tilde(int n, int k, const arma::vec& x);
void create_lambda(arma::vec& lambda,
                   double lambdamin,
                   double lambdamax,
                   double lambda_min_ratio,
                   int nsol);

#endif
