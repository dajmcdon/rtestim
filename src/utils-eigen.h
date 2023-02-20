#ifndef UTILS_EIGEN_H
#define UTILS_EIGEN_H

Eigen::SparseMatrix<double> identity(int n);
Eigen::SparseMatrix<double> get_Dtil(int k, Rcpp::NumericVector xd);
Eigen::SparseMatrix<double> get_D(int k, Rcpp::NumericVector xd);
Rcpp::NumericVector doDv(Rcpp::NumericVector v, int k, Rcpp::NumericVector xd);
Rcpp::NumericVector doDtv(Rcpp::NumericVector v, int k, Rcpp::NumericVector xd);
Rcpp::NumericVector doDtDv(Rcpp::NumericVector v, int k, Rcpp::NumericVector xd);
double one_norm(Rcpp::NumericVector const& z);
double pois_obj(int ord,
                Rcpp::NumericVector const& y,
                Rcpp::NumericVector const& x,
                Rcpp::NumericVector const& w,
                Rcpp::NumericVector& theta,
                double lambda);

Rcpp::NumericVector gaussianized_data(Rcpp::NumericVector const& y,
                                      Rcpp::NumericVector const& w,
                                      Rcpp::NumericVector& theta);

double line_search(double s,
                   double lambda,
                   double alpha,
                   double gamma,
                   Rcpp::NumericVector const& y,
                   Rcpp::NumericVector const& x,
                   Rcpp::NumericVector const& w,
                   int n,
                   int ord,
                   Rcpp::NumericVector& theta,
                   Rcpp::NumericVector& theta_old,
                   int M);

#endif
