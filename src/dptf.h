#ifndef __DPTF_H
#define __DPTF_H

Rcpp::NumericVector dptf(Rcpp::NumericVector y, double lam);
Rcpp::NumericVector weight_dptf(Rcpp::NumericVector y,
                                double lam,
                                Rcpp::NumericVector x);

#endif
