#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void binlep(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gaulep(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"binlep", (DL_FUNC) &binlep, 9},
    {"gaulep", (DL_FUNC) &gaulep, 9},
    {NULL, NULL, 0}
};

void R_init_glmlep(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
