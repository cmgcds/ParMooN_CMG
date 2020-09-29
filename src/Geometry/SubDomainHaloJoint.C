// =======================================================================
// @(#)SubDomainHaloJoint.C        
// 
// Class:       SubDomainHaloJoint
// Purpose:     superclass for boundary edges in Halo calls
//
// Author:      Sashikumaar Ganesan  10.09.09
//
// History:     
//
// =======================================================================

#include <JointEqN.h>
#include <BaseCell.h>
#include <Constants.h>
#include <Database.h>
#include <SubDomainHaloJoint.h>

// Constructors
TSubDomainHaloJoint::TSubDomainHaloJoint(TBaseCell *neighb0) :
                                    TJointEqN(neighb0)
{
  ID = SubDomainHaloJoint;
}

TSubDomainHaloJoint::TSubDomainHaloJoint(TBaseCell *neighb0, TBaseCell *neighb1) : TJointEqN(neighb0, neighb1)
{
  ID = SubDomainHaloJoint;
}


