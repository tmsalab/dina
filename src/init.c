#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP dina_DINA_Gibbs(SEXP, SEXP, SEXP, SEXP);
extern SEXP dina_DINAsim(SEXP, SEXP, SEXP, SEXP);
extern SEXP dina_rDirichlet(SEXP);
extern SEXP dina_rmultinomial(SEXP);
extern SEXP dina_update_alpha(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dina_update_sg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"dina_DINA_Gibbs",   (DL_FUNC) &dina_DINA_Gibbs,   4},
    {"dina_DINAsim",      (DL_FUNC) &dina_DINAsim,      4},
    {"dina_rDirichlet",   (DL_FUNC) &dina_rDirichlet,   1},
    {"dina_rmultinomial", (DL_FUNC) &dina_rmultinomial, 1},
    {"dina_update_alpha", (DL_FUNC) &dina_update_alpha, 8},
    {"dina_update_sg",    (DL_FUNC) &dina_update_sg,    8},
    {NULL, NULL, 0}
};

void R_init_dina(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}