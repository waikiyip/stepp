#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Fortran calls */
extern void F77_NAME(cinc)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(crstm)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"cinc",  (DL_FUNC) &F77_NAME(cinc),   7},
    {"crstm", (DL_FUNC) &F77_NAME(crstm), 18},
    {NULL, NULL, 0}
};

void R_init_stepp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
