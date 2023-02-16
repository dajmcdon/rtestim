#ifndef __DPTF_H
#define __DPTF_H

arma::vec dptf(arma::vec y, double lam);
arma::vec dptf_past(arma::vec y, double lam, arma::vec x);

#endif
