 // =======================================================================
// @(#)SubDomainJoint.C       
// 
// Class:       TSubDomainJoint
// Purpose:     superclass for edges in two/more subdomains
//
// Author:      Sashikumaar Ganesan  07.09.09
//
// History:     
//
// =======================================================================

#include <JointEqN.h>
#include <BaseCell.h>
#include <Constants.h>
#include <Database.h>
#include <RefDesc.h>
#include <SubDomainJoint.h>

// Constructors
TSubDomainJoint::TSubDomainJoint(TBaseCell *neighb0) :
                                    TJointEqN(neighb0)
{
  ID = SubDomainJoint;
  N_NeibSubDomains = 0;
  NeibSubDomainRank = -1;
  NeibSubDomainGlobalCellNo = -1;
  NeibSubDomainLocalJointNo = -1;
}

TSubDomainJoint::TSubDomainJoint(TBaseCell *neighb0, TBaseCell *neighb1) : TJointEqN(neighb0, neighb1)
{
  ID = SubDomainJoint;
  N_NeibSubDomains = 0;
  NeibSubDomainRank = -1;
  NeibSubDomainGlobalCellNo = -1;
  NeibSubDomainLocalJointNo = -1;
}

TSubDomainJoint::TSubDomainJoint(TBaseCell *neighb0, TBaseCell *neighb1, int neib_ID, 
                                 int neibSubDomainGlobalCellNo,  int neibSubDomainLocalJointNo) : TJointEqN(neighb0, neighb1)
{
  ID = SubDomainJoint;
  N_NeibSubDomains = 1;

  NeibSubDomainRank = neib_ID;
  NeibSubDomainGlobalCellNo = neibSubDomainGlobalCellNo;
  NeibSubDomainLocalJointNo = neibSubDomainLocalJointNo;
}






