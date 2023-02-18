#ifndef UTILS_EIGEN_H
#define UTILS_EIGEN_H

Eigen::SparseMatrix<double> identity(int n);
Eigen::SparseMatrix<double> get_Dtil(int k, Rcpp::NumericVector xd);
Eigen::SparseMatrix<double> get_D(int k, Rcpp::NumericVector xd);

#endif
