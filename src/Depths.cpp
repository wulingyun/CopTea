#include "CopTea.h"

inline int *GetDepth(int *edges, int *index, int nEdges, int *seed, int nGene, int maxDepth, int *depth)
{
	int *queue = (int *) R_alloc(nGene, sizeof(int));
	int queue_head, queue_tail;

	queue_head = queue_tail = 0;

	for (int i = 0; i < nGene; i++)
	{
  	if (seed[i])
    {
      queue[queue_tail++] = i;
  	  depth[i] = 0;
		}
    else
      depth[i] = -1;
	}

  while (queue_head != queue_tail)
  {
		int k = queue[queue_head++];
		int d = depth[k] + 1;
    
    if (d > maxDepth) break;

		for (int i = index[k]; i < index[k+1]; i++)
		{
			int t = edges[nEdges + i] - 1;
			if (depth[t] < 0)
			{
				queue[queue_tail++] = t;
				depth[t] = d;
			}
		}
	}

	return (depth);
}


SEXP NE_GetDepths(SEXP _Edges, SEXP _Index, SEXP _Seed, SEXP _MaxDepth)
{
  PROTECT(_Edges = AS_INTEGER(_Edges));
  int *Edges = INTEGER_POINTER(_Edges);
  PROTECT(_Index = AS_INTEGER(_Index));
  int *Index = INTEGER_POINTER(_Index);
  PROTECT(_Seed = AS_LOGICAL(_Seed));
  int *Seed = LOGICAL_POINTER(_Seed);
  PROTECT(_MaxDepth = AS_INTEGER(_MaxDepth));
  int MaxDepth = INTEGER_POINTER(_MaxDepth)[0];

  SEXP _nEdges;
	PROTECT(_nEdges = GET_DIM(_Edges));
	int nEdges = INTEGER_POINTER(_nEdges)[0];
  int nGene = length(_Seed);

  SEXP _Depth;
	PROTECT(_Depth = NEW_INTEGER(nGene));
	int *Depth = INTEGER_POINTER(_Depth);

  GetDepth(Edges, Index, nEdges, Seed, nGene, MaxDepth, Depth);
  
  UNPROTECT(5);
  return (_Depth);
}


SEXP NE_CountDepths(SEXP _Depth, SEXP _MaxDepth)
{
  PROTECT(_Depth = AS_INTEGER(_Depth));
	int *Depth = INTEGER_POINTER(_Depth);
	PROTECT(_MaxDepth = AS_INTEGER(_MaxDepth));
  int nCounter = INTEGER_POINTER(_MaxDepth)[0] + 2;
    
  SEXP _Counter;
  PROTECT(_Counter = NEW_INTEGER(nCounter));
	int *Counter = INTEGER_POINTER(_Counter);
  
  for (int i = 0; i < nCounter; i++)
    Counter[i] = 0;
  for (int i = 0; i < length(_Depth); i++)
    Counter[Depth[i]+1]++;
    
  UNPROTECT(2);
  return (_Counter);
}
