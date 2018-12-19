#include "CopTea.h"

const R_CallMethodDef callMethods[] = {
  {"NE_GetDepths", (DL_FUNC) &NE_GetDepths, 4},
  {"NE_CountDepths", (DL_FUNC) &NE_CountDepths, 2},
  {"NE_ColSums", (DL_FUNC) &NE_ColSums, 2},
  {NULL, NULL, 0}
};

void R_init_CopTea(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}
