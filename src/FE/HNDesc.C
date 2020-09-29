// =======================================================================
// @(#)HNDesc.C        1.1 10/30/98
//
// Class:       THNDesc
// Purpose:     describe fixed information for a hanging node
//
// Author:      Gunar Matthies  18.11.97
//
// =======================================================================

#include <HNDesc.h>

/** constructor, filling all data */
THNDesc::THNDesc(int n_nodes, double *coeff)
{
  N_Nodes=n_nodes;

  Coeff=coeff;
}

