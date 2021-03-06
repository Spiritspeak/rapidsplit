// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// applyItersplits
List applyItersplits(int iters, List splits);
RcppExport SEXP _rapidsplit_applyItersplits(SEXP itersSEXP, SEXP splitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type iters(itersSEXP);
    Rcpp::traits::input_parameter< List >::type splits(splitsSEXP);
    rcpp_result_gen = Rcpp::wrap(applyItersplits(iters, splits));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rapidsplit_applyItersplits", (DL_FUNC) &_rapidsplit_applyItersplits, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_rapidsplit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
