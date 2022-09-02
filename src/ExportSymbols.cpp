#include <Rcpp.h>

RcppExport SEXP _DFM_SKF(SEXP XSEXP, SEXP ASEXP, SEXP CSEXP, SEXP QSEXP, SEXP RSEXP, SEXP F_0SEXP, SEXP P_0SEXP, SEXP retLLSEXP);
RcppExport SEXP _DFM_FIS(SEXP ASEXP, SEXP ZTfSEXP, SEXP ZTpSEXP, SEXP VTf_vSEXP, SEXP VTp_vSEXP, SEXP F_0SEXP, SEXP P_0SEXP);
RcppExport SEXP _DFM_SKFS(SEXP XSEXP, SEXP ASEXP, SEXP CSEXP, SEXP QSEXP, SEXP RSEXP, SEXP F_0SEXP, SEXP P_0SEXP, SEXP retLLSEXP);
RcppExport SEXP _DFM_Estep(SEXP XSEXP, SEXP ASEXP, SEXP CSEXP, SEXP QSEXP, SEXP RSEXP, SEXP F_0SEXP, SEXP P_0SEXP);
RcppExport SEXP _DFM_ainv(SEXP FsEXP);
RcppExport SEXP _DFM_apinv(SEXP FsEXP);

static const R_CallMethodDef CallEntries[] = {
  {"Cpp_SKF", (DL_FUNC) &_DFM_SKF, 8},
  {"Cpp_FIS", (DL_FUNC) &_DFM_FIS, 7},
  {"Cpp_SKFS", (DL_FUNC) &_DFM_SKFS, 8},
  {"Cpp_Estep", (DL_FUNC) &_DFM_Estep, 7},
  {"Cpp_ainv", (DL_FUNC) &_DFM_ainv, 1},
  {"Cpp_apinv", (DL_FUNC) &_DFM_apinv, 1},
  {NULL, NULL, 0}
};

RcppExport void R_init_DFM(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
