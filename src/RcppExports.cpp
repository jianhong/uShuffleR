// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rushuffle
CharacterVector rushuffle(CharacterVector x, IntegerVector k, IntegerVector n);
RcppExport SEXP _uShuffleR_rushuffle(SEXP xSEXP, SEXP kSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(rushuffle(x, k, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_uShuffleR_rushuffle", (DL_FUNC) &_uShuffleR_rushuffle, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_uShuffleR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
