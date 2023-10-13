// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// prox_newton_testing
Rcpp::List prox_newton_testing(int M, int Minner, int Mline, int korder, Rcpp::NumericVector const& y, Rcpp::NumericVector const& x, Rcpp::NumericVector const& w, double lambda, double ls_alpha, double ls_gamma, double tol);
RcppExport SEXP _rtestim_prox_newton_testing(SEXP MSEXP, SEXP MinnerSEXP, SEXP MlineSEXP, SEXP korderSEXP, SEXP ySEXP, SEXP xSEXP, SEXP wSEXP, SEXP lambdaSEXP, SEXP ls_alphaSEXP, SEXP ls_gammaSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type Minner(MinnerSEXP);
    Rcpp::traits::input_parameter< int >::type Mline(MlineSEXP);
    Rcpp::traits::input_parameter< int >::type korder(korderSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector const& >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector const& >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type ls_alpha(ls_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type ls_gamma(ls_gammaSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(prox_newton_testing(M, Minner, Mline, korder, y, x, w, lambda, ls_alpha, ls_gamma, tol));
    return rcpp_result_gen;
END_RCPP
}
// dptf
Rcpp::NumericVector dptf(Rcpp::NumericVector y, double lam);
RcppExport SEXP _rtestim_dptf(SEXP ySEXP, SEXP lamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    rcpp_result_gen = Rcpp::wrap(dptf(y, lam));
    return rcpp_result_gen;
END_RCPP
}
// weight_dptf
Rcpp::NumericVector weight_dptf(Rcpp::NumericVector y, double lam, Rcpp::NumericVector w);
RcppExport SEXP _rtestim_weight_dptf(SEXP ySEXP, SEXP lamSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(weight_dptf(y, lam, w));
    return rcpp_result_gen;
END_RCPP
}
// rtestim_path
List rtestim_path(NumericVector y, NumericVector x, NumericVector w, int korder, NumericVector lambda, double lambdamax, double lambdamin, int nsol, double rho, int maxiter, int maxiter_newton, int maxiter_line, double tolerance, double lambda_min_ratio, double ls_alpha, double ls_gamma, int verbose);
RcppExport SEXP _rtestim_rtestim_path(SEXP ySEXP, SEXP xSEXP, SEXP wSEXP, SEXP korderSEXP, SEXP lambdaSEXP, SEXP lambdamaxSEXP, SEXP lambdaminSEXP, SEXP nsolSEXP, SEXP rhoSEXP, SEXP maxiterSEXP, SEXP maxiter_newtonSEXP, SEXP maxiter_lineSEXP, SEXP toleranceSEXP, SEXP lambda_min_ratioSEXP, SEXP ls_alphaSEXP, SEXP ls_gammaSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type korder(korderSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lambdamax(lambdamaxSEXP);
    Rcpp::traits::input_parameter< double >::type lambdamin(lambdaminSEXP);
    Rcpp::traits::input_parameter< int >::type nsol(nsolSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter_newton(maxiter_newtonSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter_line(maxiter_lineSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_min_ratio(lambda_min_ratioSEXP);
    Rcpp::traits::input_parameter< double >::type ls_alpha(ls_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type ls_gamma(ls_gammaSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(rtestim_path(y, x, w, korder, lambda, lambdamax, lambdamin, nsol, rho, maxiter, maxiter_newton, maxiter_line, tolerance, lambda_min_ratio, ls_alpha, ls_gamma, verbose));
    return rcpp_result_gen;
END_RCPP
}
// get_Dtil
Eigen::SparseMatrix<double> get_Dtil(int k, NumericVector xd);
RcppExport SEXP _rtestim_get_Dtil(SEXP kSEXP, SEXP xdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xd(xdSEXP);
    rcpp_result_gen = Rcpp::wrap(get_Dtil(k, xd));
    return rcpp_result_gen;
END_RCPP
}
// get_D
Eigen::SparseMatrix<double> get_D(int k, NumericVector xd);
RcppExport SEXP _rtestim_get_D(SEXP kSEXP, SEXP xdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xd(xdSEXP);
    rcpp_result_gen = Rcpp::wrap(get_D(k, xd));
    return rcpp_result_gen;
END_RCPP
}
// create_lambda_test
NumericVector create_lambda_test(NumericVector lambda, double lambdamin, double lambdamax, double lambda_min_ratio, int nsol);
RcppExport SEXP _rtestim_create_lambda_test(SEXP lambdaSEXP, SEXP lambdaminSEXP, SEXP lambdamaxSEXP, SEXP lambda_min_ratioSEXP, SEXP nsolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lambdamin(lambdaminSEXP);
    Rcpp::traits::input_parameter< double >::type lambdamax(lambdamaxSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_min_ratio(lambda_min_ratioSEXP);
    Rcpp::traits::input_parameter< int >::type nsol(nsolSEXP);
    rcpp_result_gen = Rcpp::wrap(create_lambda_test(lambda, lambdamin, lambdamax, lambda_min_ratio, nsol));
    return rcpp_result_gen;
END_RCPP
}
// doDv
NumericVector doDv(NumericVector v, int k, NumericVector xd);
RcppExport SEXP _rtestim_doDv(SEXP vSEXP, SEXP kSEXP, SEXP xdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xd(xdSEXP);
    rcpp_result_gen = Rcpp::wrap(doDv(v, k, xd));
    return rcpp_result_gen;
END_RCPP
}
// doDtv
NumericVector doDtv(NumericVector v, int k, NumericVector xd);
RcppExport SEXP _rtestim_doDtv(SEXP vSEXP, SEXP kSEXP, SEXP xdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xd(xdSEXP);
    rcpp_result_gen = Rcpp::wrap(doDtv(v, k, xd));
    return rcpp_result_gen;
END_RCPP
}
// one_norm
double one_norm(NumericVector const& z);
RcppExport SEXP _rtestim_one_norm(SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector const& >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(one_norm(z));
    return rcpp_result_gen;
END_RCPP
}
// pois_obj
double pois_obj(int ord, NumericVector const& y, NumericVector const& x, NumericVector const& w, NumericVector& theta, double lambda);
RcppExport SEXP _rtestim_pois_obj(SEXP ordSEXP, SEXP ySEXP, SEXP xSEXP, SEXP wSEXP, SEXP thetaSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ord(ordSEXP);
    Rcpp::traits::input_parameter< NumericVector const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector const& >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector const& >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(pois_obj(ord, y, x, w, theta, lambda));
    return rcpp_result_gen;
END_RCPP
}
// gaussianized_data
NumericVector gaussianized_data(NumericVector const& y, NumericVector const& w, NumericVector& theta);
RcppExport SEXP _rtestim_gaussianized_data(SEXP ySEXP, SEXP wSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector const& >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(gaussianized_data(y, w, theta));
    return rcpp_result_gen;
END_RCPP
}
// line_search
double line_search(double s, double lambda, double alpha, double gamma, NumericVector const& y, NumericVector const& x, NumericVector const& w, int n, int ord, NumericVector& theta, NumericVector& theta_old, int M);
RcppExport SEXP _rtestim_line_search(SEXP sSEXP, SEXP lambdaSEXP, SEXP alphaSEXP, SEXP gammaSEXP, SEXP ySEXP, SEXP xSEXP, SEXP wSEXP, SEXP nSEXP, SEXP ordSEXP, SEXP thetaSEXP, SEXP theta_oldSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< NumericVector const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector const& >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector const& >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type ord(ordSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type theta_old(theta_oldSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(line_search(s, lambda, alpha, gamma, y, x, w, n, ord, theta, theta_old, M));
    return rcpp_result_gen;
END_RCPP
}
// compute_gcd
int compute_gcd(IntegerVector x);
RcppExport SEXP _rtestim_compute_gcd(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_gcd(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rtestim_prox_newton_testing", (DL_FUNC) &_rtestim_prox_newton_testing, 11},
    {"_rtestim_dptf", (DL_FUNC) &_rtestim_dptf, 2},
    {"_rtestim_weight_dptf", (DL_FUNC) &_rtestim_weight_dptf, 3},
    {"_rtestim_rtestim_path", (DL_FUNC) &_rtestim_rtestim_path, 17},
    {"_rtestim_get_Dtil", (DL_FUNC) &_rtestim_get_Dtil, 2},
    {"_rtestim_get_D", (DL_FUNC) &_rtestim_get_D, 2},
    {"_rtestim_create_lambda_test", (DL_FUNC) &_rtestim_create_lambda_test, 5},
    {"_rtestim_doDv", (DL_FUNC) &_rtestim_doDv, 3},
    {"_rtestim_doDtv", (DL_FUNC) &_rtestim_doDtv, 3},
    {"_rtestim_one_norm", (DL_FUNC) &_rtestim_one_norm, 1},
    {"_rtestim_pois_obj", (DL_FUNC) &_rtestim_pois_obj, 6},
    {"_rtestim_gaussianized_data", (DL_FUNC) &_rtestim_gaussianized_data, 3},
    {"_rtestim_line_search", (DL_FUNC) &_rtestim_line_search, 12},
    {"_rtestim_compute_gcd", (DL_FUNC) &_rtestim_compute_gcd, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_rtestim(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
