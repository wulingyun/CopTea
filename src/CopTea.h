#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
using namespace std;

/* Interfaces to R */

extern "C" {
  /* DLL Init */
  void R_init_CopTea(DllInfo *info);

	/* Utils */
  SEXP NE_GetDepths(SEXP _Edges, SEXP _Index, SEXP _Core, SEXP _MaxDepth);
  SEXP NE_CountDepths(SEXP _Depth, SEXP _MaxDepth);
  
  SEXP NE_ColSums(SEXP _Matrix, SEXP _RowSel);
}


/* initialize the list */
#define setValues(r, c, v) for (int i = 0; i < length(r); i++) c[i] = v

/* get/set the list element */
SEXP getListElement(SEXP list, const char *tag);
void setListElement(SEXP list, int i, const char *tag, SEXP value);

/* set dim of array */
void setDim2(SEXP array, int x1, int x2);
void setDim3(SEXP array, int x1, int x2, int x3);
void setDim4(SEXP array, int x1, int x2, int x3, int x4);

/* sampling method */
int *sampleWithoutReplace(int n, int k, int *result = NULL, int *buffer = NULL);
