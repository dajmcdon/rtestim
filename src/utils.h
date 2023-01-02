#ifndef __UTILS_H
#define __UTILS_H

double chooseC(int n, int k);
arma::sp_mat buildD(int m, int n, int ord, int wrap);
arma::sp_mat buildDx(int k, const arma::vec& x);
arma::sp_mat buildDx_tilde(int k, const arma::vec& x);
void EntrywiseSoftThreshold(arma::vec& z, double lam);
void create_lambda(arma::vec& lambda, double lambdamin, double lambdamax,
                   double lambda_min_ratio, int nsol);

#endif
