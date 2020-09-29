// =======================================================================
// @(#)PeriodicJoint.C        1.1 04/21/99
// 
// Class:       TPeriodicJoint
// Purpose:     connects two cells with periodic boundary conditions
//
// Author:      Volker Behns  16.04.99
//
// =======================================================================

#include <PeriodicJoint.h>
#include <BaseCell.h>

// Constructors
TPeriodicJoint::TPeriodicJoint(TBaseCell *neighb0) :
                TJointEqN(neighb0)
{
  ID = PeriodicJoint;
}

TPeriodicJoint::TPeriodicJoint(TBaseCell *neighb0, TBaseCell *neighb1) :
                TJointEqN(neighb0, neighb1)
{
  ID = PeriodicJoint;
}

// Methods
TJoint *TPeriodicJoint::NewInst(double newT_0, double newT_1, TBaseCell *Me)
{
  return new TPeriodicJoint(Me);
}

TJoint *TPeriodicJoint::NewInst()
{
  return new TPeriodicJoint(NULL);
}

