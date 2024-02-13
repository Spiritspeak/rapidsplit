// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// colMedians
NumericVector colMedians(NumericMatrix mat);
RcppExport SEXP _rapidsplit_colMedians(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(colMedians(mat));
    return rcpp_result_gen;
END_RCPP
}
// colMedians_mask
NumericVector colMedians_mask(NumericMatrix mat, LogicalMatrix mask);
RcppExport SEXP _rapidsplit_colMedians_mask(SEXP matSEXP, SEXP maskSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type mask(maskSEXP);
    rcpp_result_gen = Rcpp::wrap(colMedians_mask(mat, mask));
    return rcpp_result_gen;
END_RCPP
}
// mediansByMask
NumericVector mediansByMask(NumericVector values, LogicalMatrix mask);
RcppExport SEXP _rapidsplit_mediansByMask(SEXP valuesSEXP, SEXP maskSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type values(valuesSEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type mask(maskSEXP);
    rcpp_result_gen = Rcpp::wrap(mediansByMask(values, mask));
    return rcpp_result_gen;
END_RCPP
}
// colMeans_mask
NumericVector colMeans_mask(NumericMatrix mat, LogicalMatrix mask);
RcppExport SEXP _rapidsplit_colMeans_mask(SEXP matSEXP, SEXP maskSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type mask(maskSEXP);
    rcpp_result_gen = Rcpp::wrap(colMeans_mask(mat, mask));
    return rcpp_result_gen;
END_RCPP
}
// meansByMask
NumericVector meansByMask(NumericVector values, LogicalMatrix mask);
RcppExport SEXP _rapidsplit_meansByMask(SEXP valuesSEXP, SEXP maskSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type values(valuesSEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type mask(maskSEXP);
    rcpp_result_gen = Rcpp::wrap(meansByMask(values, mask));
    return rcpp_result_gen;
END_RCPP
}
// colSds
NumericVector colSds(NumericMatrix mat);
RcppExport SEXP _rapidsplit_colSds(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(colSds(mat));
    return rcpp_result_gen;
END_RCPP
}
// colSds_mask
NumericVector colSds_mask(NumericMatrix mat, LogicalMatrix mask);
RcppExport SEXP _rapidsplit_colSds_mask(SEXP matSEXP, SEXP maskSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type mask(maskSEXP);
    rcpp_result_gen = Rcpp::wrap(colSds_mask(mat, mask));
    return rcpp_result_gen;
END_RCPP
}
// sdsByMask
NumericVector sdsByMask(NumericVector values, LogicalMatrix mask);
RcppExport SEXP _rapidsplit_sdsByMask(SEXP valuesSEXP, SEXP maskSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type values(valuesSEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type mask(maskSEXP);
    rcpp_result_gen = Rcpp::wrap(sdsByMask(values, mask));
    return rcpp_result_gen;
END_RCPP
}
// ExcludeSDOutliers
LogicalMatrix ExcludeSDOutliers(NumericVector rtvec, LogicalMatrix mask, double sdlim);
RcppExport SEXP _rapidsplit_ExcludeSDOutliers(SEXP rtvecSEXP, SEXP maskSEXP, SEXP sdlimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rtvec(rtvecSEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< double >::type sdlim(sdlimSEXP);
    rcpp_result_gen = Rcpp::wrap(ExcludeSDOutliers(rtvec, mask, sdlim));
    return rcpp_result_gen;
END_RCPP
}
// ExcludeSDOutliers_nomask
LogicalMatrix ExcludeSDOutliers_nomask(NumericMatrix mat, double sdlim);
RcppExport SEXP _rapidsplit_ExcludeSDOutliers_nomask(SEXP matSEXP, SEXP sdlimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< double >::type sdlim(sdlimSEXP);
    rcpp_result_gen = Rcpp::wrap(ExcludeSDOutliers_nomask(mat, sdlim));
    return rcpp_result_gen;
END_RCPP
}
// applyItersplits
List applyItersplits(int iters, List splits, bool replace);
RcppExport SEXP _rapidsplit_applyItersplits(SEXP itersSEXP, SEXP splitsSEXP, SEXP replaceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type iters(itersSEXP);
    Rcpp::traits::input_parameter< List >::type splits(splitsSEXP);
    Rcpp::traits::input_parameter< bool >::type replace(replaceSEXP);
    rcpp_result_gen = Rcpp::wrap(applyItersplits(iters, splits, replace));
    return rcpp_result_gen;
END_RCPP
}
// stratified_itersplits
IntegerMatrix stratified_itersplits(int itercount, IntegerVector groupsizes);
RcppExport SEXP _rapidsplit_stratified_itersplits(SEXP itercountSEXP, SEXP groupsizesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type itercount(itercountSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type groupsizes(groupsizesSEXP);
    rcpp_result_gen = Rcpp::wrap(stratified_itersplits(itercount, groupsizes));
    return rcpp_result_gen;
END_RCPP
}
// corByColumns
NumericVector corByColumns(NumericMatrix mat1, NumericMatrix mat2);
RcppExport SEXP _rapidsplit_corByColumns(SEXP mat1SEXP, SEXP mat2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat1(mat1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mat2(mat2SEXP);
    rcpp_result_gen = Rcpp::wrap(corByColumns(mat1, mat2));
    return rcpp_result_gen;
END_RCPP
}
// corByColumns_mask
NumericVector corByColumns_mask(NumericMatrix mat1, NumericMatrix mat2, LogicalMatrix mask);
RcppExport SEXP _rapidsplit_corByColumns_mask(SEXP mat1SEXP, SEXP mat2SEXP, SEXP maskSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat1(mat1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mat2(mat2SEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type mask(maskSEXP);
    rcpp_result_gen = Rcpp::wrap(corByColumns_mask(mat1, mat2, mask));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rapidsplit_colMedians", (DL_FUNC) &_rapidsplit_colMedians, 1},
    {"_rapidsplit_colMedians_mask", (DL_FUNC) &_rapidsplit_colMedians_mask, 2},
    {"_rapidsplit_mediansByMask", (DL_FUNC) &_rapidsplit_mediansByMask, 2},
    {"_rapidsplit_colMeans_mask", (DL_FUNC) &_rapidsplit_colMeans_mask, 2},
    {"_rapidsplit_meansByMask", (DL_FUNC) &_rapidsplit_meansByMask, 2},
    {"_rapidsplit_colSds", (DL_FUNC) &_rapidsplit_colSds, 1},
    {"_rapidsplit_colSds_mask", (DL_FUNC) &_rapidsplit_colSds_mask, 2},
    {"_rapidsplit_sdsByMask", (DL_FUNC) &_rapidsplit_sdsByMask, 2},
    {"_rapidsplit_ExcludeSDOutliers", (DL_FUNC) &_rapidsplit_ExcludeSDOutliers, 3},
    {"_rapidsplit_ExcludeSDOutliers_nomask", (DL_FUNC) &_rapidsplit_ExcludeSDOutliers_nomask, 2},
    {"_rapidsplit_applyItersplits", (DL_FUNC) &_rapidsplit_applyItersplits, 3},
    {"_rapidsplit_stratified_itersplits", (DL_FUNC) &_rapidsplit_stratified_itersplits, 2},
    {"_rapidsplit_corByColumns", (DL_FUNC) &_rapidsplit_corByColumns, 2},
    {"_rapidsplit_corByColumns_mask", (DL_FUNC) &_rapidsplit_corByColumns_mask, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_rapidsplit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}