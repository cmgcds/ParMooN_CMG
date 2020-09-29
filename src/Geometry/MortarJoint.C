// =======================================================================
// @(#)MortarJoint.C        1.1 10/30/98
// 
// Class:       TMortarJoint
// Purpose:     indicates a mortar joint
//
// Author:      Volker Behns  19.03.98
//
// =======================================================================

#ifndef __MORTAR__
#define __MORTAR__
#endif

#include <MortarJoint.h>

// Constructors
TMortarJoint::TMortarJoint()
{
  ID = MortarJoint;
}

// Methods
int TMortarJoint::CheckMatchingRef(TBaseCell *Me, int J_i,
                    struct StoreGeom &Tmp)
{
  Tmp.Filled = FALSE;
  return 0;
}

#ifdef __MORTAR__
int TMortarJoint::CheckMatchingRef(TBaseCell *Me, int J_i,
                    StoreGeomMortar &Tmp)
{
  Tmp.Filled = FALSE;
  return 0;
}

#endif

