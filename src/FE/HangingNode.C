// =======================================================================
// @(#)HangingNode.C        1.1 10/30/98
//
// Class:       THangingNode
// Purpose:     represent a hanging node
//
// Author:      Gunar Matthies  18.11.97
//
// =======================================================================

#include <HangingNode.h>

/** constructor, filling all data */
THangingNode::THangingNode(HNDesc type, int n_, int *all_dofs, int *coupling)
{
  static int i, N_;

  N_ = n_;
  Type=type;

  DOF=new int[N_];

  for(i=0;i<N_;i++)
    DOF[i]=all_dofs[coupling[i]];
}

/** destructor */
THangingNode::~THangingNode()
{
  delete DOF;
}

