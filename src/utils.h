#ifndef UTILS_H
#define UTILS_H

Rcpp::NumericVector evec_to_nvec(Eigen::VectorXd evec);
Eigen::VectorXd nvec_to_evec(Rcpp::NumericVector nvec);
Eigen::SparseMatrix<double> identity(int n);
Eigen::SparseMatrix<double> get_Dtil(int k, Rcpp::NumericVector xd);
Eigen::SparseMatrix<double> get_D(int k, Rcpp::NumericVector xd);
Rcpp::NumericVector doDv(Rcpp::NumericVector v, int k, Rcpp::NumericVector xd);
Rcpp::NumericVector doDtv(Rcpp::NumericVector v, int k, Rcpp::NumericVector xd);
Rcpp::NumericVector doDtDv(Rcpp::NumericVector v, int k, Rcpp::NumericVector xd);
void create_lambda(Rcpp::NumericVector& lambda,
                   double& lambdamin,
                   double& lambdamax,
                   double& lambda_min_ratio,
                   int& nsol);

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
