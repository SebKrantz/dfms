#include <Rcpp.h>

RcppExport SEXP _DFM_KalmanFilter(SEXP XSEXP, SEXP ASEXP, SEXP CSEXP, SEXP QSEXP, SEXP RSEXP, SEXP F0SEXP, SEXP P0SEXP);
RcppExport SEXP _DFM_KalmanSmoother(SEXP ASEXP, SEXP CSEXP, SEXP RSEXP, SEXP ZTfSEXP, SEXP ZTpSEXP, SEXP VTf_vSEXP, SEXP VTp_vSEXP);
RcppExport SEXP _DFM_KalmanFilterSmoother(SEXP XSEXP, SEXP ASEXP, SEXP CSEXP, SEXP QSEXP, SEXP RSEXP, SEXP F0SEXP, SEXP P0SEXP);
RcppExport SEXP _DFM_Estep(SEXP XSEXP, SEXP CSEXP, SEXP QSEXP, SEXP RSEXP, SEXP ASEXP, SEXP F0SEXP, SEXP P0SEXP);
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
