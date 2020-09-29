// =======================================================================
// @(#)MortarBaseJoint.C        1.1 10/30/98
// 
// Class:       TMortarBaseJoint
// Purpose:     mortar joint on macro grid level
//
// Author:      Volker Behns  19.03.98
//
// =======================================================================

#ifndef __MORTAR__
#define __MORTAR__
#endif

#include <MortarBaseJoint.h>

// Constructors
TMortarBaseJoint::TMortarBaseJoint(TBaseCell *neighb0,
                    TBaseCell *neighb1) : TJointEqN(neighb0, neighb1)
{
  ID = MortarBaseJoint;
}

// Methods
int TMortarBaseJoint::CheckMatchingRef(TBaseCell *Me, int J_i,
                    struct StoreGeom &Tmp)
{
  Tmp.Filled = FALSE;
  return 0;
}

#ifdef __MORTAR__
int TMortarBaseJoint::CheckMatchingRef(TBaseCell *Me, int J_i,
                    StoreGeomMortar &Tmp)
{
  Tmp.Filled = FALSE;
  return 0;
}

#endif

