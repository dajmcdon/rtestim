#ifndef __DPTF_H
#define __DPTF_H

arma::vec dptf(arma::vec y, double lam);
arma::vec dptf_weight(arma::vec y, double lam, arma::vec w);
arma::vec dptf_past_weight(arma::vec y, double lam, arma::vec x, arma::vec w);

#endif
