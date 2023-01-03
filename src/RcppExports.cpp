// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// admm_testing
List admm_testing(int M, arma::vec const& y, arma::vec const& w, int n, arma::vec theta, arma::vec z, arma::vec u, double lambda, double rho, double mu, arma::sp_mat const& DD, arma::sp_mat const& D, double tol);
RcppExport SEXP _rtestim_admm_testing(SEXP MSEXP, SEXP ySEXP, SEXP wSEXP, SEXP nSEXP, SEXP thetaSEXP, SEXP zSEXP, SEXP uSEXP, SEXP lambdaSEXP, SEXP rhoSEXP, SEXP muSEXP, SEXP DDSEXP, SEXP DSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat const& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat const& >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(admm_testing(M, y, w, n, theta, z, u, lambda, rho, mu, DD, D, tol));
    return rcpp_result_gen;
END_RCPP
}
// dptf
arma::vec dptf(arma::vec y, double lam);
RcppExport SEXP _rtestim_dptf(SEXP ySEXP, SEXP lamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    rcpp_result_gen = Rcpp::wrap(dptf(y, lam));
    return rcpp_result_gen;
END_RCPP
}
// rtestim_path
List rtestim_path(arma::vec y, arma::vec x, arma::vec w, int korder, arma::vec lambda, double lambdamax, double lambdamin, int nsol, double rho_adjust, double rho, int maxiter, double tolerance, double lambda_min_ratio, int verbose);
RcppExport SEXP _rtestim_rtestim_path(SEXP ySEXP, SEXP xSEXP, SEXP wSEXP, SEXP korderSEXP, SEXP lambdaSEXP, SEXP lambdamaxSEXP, SEXP lambdaminSEXP, SEXP nsolSEXP, SEXP rho_adjustSEXP, SEXP rhoSEXP, SEXP maxiterSEXP, SEXP toleranceSEXP, SEXP lambda_min_ratioSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type korder(korderSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lambdamax(lambdamaxSEXP);
    Rcpp::traits::input_parameter< double >::type lambdamin(lambdaminSEXP);
    Rcpp::traits::input_parameter< int >::type nsol(nsolSEXP);
    Rcpp::traits::input_parameter< double >::type rho_adjust(rho_adjustSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_min_ratio(lambda_min_ratioSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(rtestim_path(y, x, w, korder, lambda, lambdamax, lambdamin, nsol, rho_adjust, rho, maxiter, tolerance, lambda_min_ratio, verbose));
    return rcpp_result_gen;
END_RCPP
}
// buildD
arma::sp_mat buildD(int n, int ord);
RcppExport SEXP _rtestim_buildD(SEXP nSEXP, SEXP ordSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type ord(ordSEXP);
    rcpp_result_gen = Rcpp::wrap(buildD(n, ord));
    return rcpp_result_gen;
END_RCPP
}
// buildDx
arma::sp_mat buildDx(int k, const arma::vec& x);
RcppExport SEXP _rtestim_buildDx(SEXP kSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(buildDx(k, x));
    return rcpp_result_gen;
END_RCPP
}
// buildDx_tilde
arma::sp_mat buildDx_tilde(int k, const arma::vec& x);
RcppExport SEXP _rtestim_buildDx_tilde(SEXP kSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(buildDx_tilde(k, x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rtestim_admm_testing", (DL_FUNC) &_rtestim_admm_testing, 13},
    {"_rtestim_dptf", (DL_FUNC) &_rtestim_dptf, 2},
    {"_rtestim_rtestim_path", (DL_FUNC) &_rtestim_rtestim_path, 14},
    {"_rtestim_buildD", (DL_FUNC) &_rtestim_buildD, 2},
    {"_rtestim_buildDx", (DL_FUNC) &_rtestim_buildDx, 2},
    {"_rtestim_buildDx_tilde", (DL_FUNC) &_rtestim_buildDx_tilde, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_rtestim(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
