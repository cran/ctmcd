#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP ctmcd_RcppExport_registerCCallable(void);
extern SEXP ctmcd_rNijTRiT_ModRej(SEXP, SEXP, SEXP);
extern SEXP ctmcd_rNijTRiT_Unif(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"ctmcd_RcppExport_registerCCallable", (DL_FUNC) &ctmcd_RcppExport_registerCCallable, 0},
  {"ctmcd_rNijTRiT_ModRej",              (DL_FUNC) &ctmcd_rNijTRiT_ModRej,              3},
  {"ctmcd_rNijTRiT_Unif",                (DL_FUNC) &ctmcd_rNijTRiT_Unif,                4},
  {NULL, NULL, 0}
};

void R_init_ctmcd(DllInfo *dll)
{
  R_RegisterCCallable("ctmcd", "ctmcd_rNijTRiT_ModRej", (DL_FUNC) &ctmcd_rNijTRiT_ModRej);
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
