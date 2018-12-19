#include "CopTea.h"

SEXP NE_ColSums(SEXP _Matrix, SEXP _RowSel)
{
  PROTECT(_Matrix = AS_NUMERIC(_Matrix));
  double *Matrix = NUMERIC_POINTER(_Matrix);
  PROTECT(_RowSel = AS_INTEGER(_RowSel));
  int *RowSel = INTEGER_POINTER(_RowSel);
  int nRowSel = length(_RowSel);

  SEXP _Dim;
  PROTECT(_Dim = GET_DIM(_Matrix));
  int *Dim = INTEGER_POINTER(_Dim);
  int nRow = Dim[0];
  int nCol = Dim[1];

  SEXP _Sums;
  PROTECT(_Sums = NEW_NUMERIC(nCol));
  double *Sums = NUMERIC_POINTER(_Sums);

  for (int i = 0; i < nCol; i++)
  {
    Sums[i] = 0;
    for (int j = 0; j < nRowSel; j++)
    {
      Sums[i] += Matrix[i*nRow + RowSel[j] - 1];
    }
  }

  UNPROTECT(4);
  return (_Sums);
}
