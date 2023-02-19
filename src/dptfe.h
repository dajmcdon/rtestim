#ifndef __DPTFE_H
#define __DPTFE_H

Rcpp::NumericVector dptfe(Rcpp::NumericVector y, double lam);
Rcpp::NumericVector dptfe_past(Rcpp::NumericVector y, double lam,
                              Rcpp::NumericVector x);


#endif
