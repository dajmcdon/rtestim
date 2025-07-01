#ifndef __DPTF_H
#define __DPTF_H

Rcpp::NumericVector rcpp_tvdz(
    Rcpp::NumericVector y,
    Rcpp::NumericVector z,
    double lam
);

/* Unused, but provided */
Rcpp::NumericVector rcpp_wtvdz(
    Rcpp::NumericVector y,
    Rcpp::NumericVector z,
    double lam, Rcpp::NumericVector x
);

#endif
