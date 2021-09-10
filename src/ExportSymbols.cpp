#include <Rcpp.h>

RcppExport SEXP _DFM_KalmanFilter(SEXP ySEXP, SEXP HSEXP, SEXP QSEXP, SEXP RSEXP, SEXP FsEXP, SEXP F0SEXP, SEXP P0SEXP);
RcppExport SEXP _DFM_KalmanSmoother(SEXP FsEXP, SEXP HSEXP, SEXP RSEXP, SEXP FfTSEXP, SEXP FpTSEXP, SEXP PfT_vSEXP, SEXP PpT_vSEXP);
RcppExport SEXP _DFM_KalmanFilterSmoother(SEXP ySEXP, SEXP HSEXP, SEXP QSEXP, SEXP RSEXP, SEXP FsEXP, SEXP F0SEXP, SEXP P0SEXP);
RcppExport SEXP _DFM_Estep(SEXP ySEXP, SEXP HSEXP, SEXP QSEXP, SEXP RSEXP, SEXP FsEXP, SEXP F0SEXP, SEXP P0SEXP);
RcppExport SEXP _DFM_ainv(SEXP FsEXP);
RcppExport SEXP _DFM_apinv(SEXP FsEXP);

static const R_CallMethodDef CallEntries[] = {
  {"Cpp_KalmanFilter",   (DL_FUNC) &_DFM_KalmanFilter,   7},
  {"Cpp_KalmanSmoother", (DL_FUNC) &_DFM_KalmanSmoother, 7},
  {"Cpp_KalmanFilterSmoother", (DL_FUNC) &_DFM_KalmanFilterSmoother, 7},
  {"Cpp_Estep",          (DL_FUNC) &_DFM_Estep,          7},
  {"Cpp_ainv",        (DL_FUNC) &_DFM_ainv,        1},
  {"Cpp_apinv",       (DL_FUNC) &_DFM_apinv,       1},
  {NULL, NULL, 0}
};

RcppExport void R_init_DFM(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
