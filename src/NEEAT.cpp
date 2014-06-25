#include "NEEAT.h"

const R_CallMethodDef callMethods[] = {
  {"RMultiHyper", (DL_FUNC) &RMultiHyper, 3},
  {"PMultiHyper", (DL_FUNC) &PMultiHyper, 4},
  {"PMultiNom", (DL_FUNC) &PMultiNom, 4},
  {"NE_GetDepths", (DL_FUNC) &NE_GetDepths, 4},
  {"NE_CountDepths", (DL_FUNC) &NE_CountDepths, 2},
  {NULL, NULL, 0}
};

void R_init_NEEAT(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
