// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// inorm
double inorm(double x);
RcppExport SEXP _GCGK_inorm(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(inorm(x));
    return rcpp_result_gen;
END_RCPP
}
// V
double V(double a, double b, double c, double d, double s);
RcppExport SEXP _GCGK_V(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(V(a, b, c, d, s));
    return rcpp_result_gen;
END_RCPP
}
// GCGK
NumericVector GCGK(NumericMatrix X, NumericVector parameters);
RcppExport SEXP _GCGK_GCGK(SEXP XSEXP, SEXP parametersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type parameters(parametersSEXP);
    rcpp_result_gen = Rcpp::wrap(GCGK(X, parameters));
    return rcpp_result_gen;
END_RCPP
}
// GCGKTEST
List GCGKTEST(NumericMatrix X, NumericVector parameters, int iteration);
RcppExport SEXP _GCGK_GCGKTEST(SEXP XSEXP, SEXP parametersSEXP, SEXP iterationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< int >::type iteration(iterationSEXP);
    rcpp_result_gen = Rcpp::wrap(GCGKTEST(X, parameters, iteration));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GCGK_inorm", (DL_FUNC) &_GCGK_inorm, 1},
    {"_GCGK_V", (DL_FUNC) &_GCGK_V, 5},
    {"_GCGK_GCGK", (DL_FUNC) &_GCGK_GCGK, 2},
    {"_GCGK_GCGKTEST", (DL_FUNC) &_GCGK_GCGKTEST, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_GCGK(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
