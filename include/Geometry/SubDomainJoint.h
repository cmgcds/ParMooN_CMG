// =======================================================================
// @(#)SubDomainJoint.h        
// 
// Class:       TSubDomainJoint
// Purpose:     superclass for edges in two/more subdomains
//
// Author:      Sashikumaar Ganesan  07.09.09
//
// History:     
//
// =======================================================================

#ifndef __SUBDOMAINJOINT__
#define __SUBDOMAINJOINT__

#include <stdlib.h>
#include <JointEqN.h>

/** edges which are in interface between two different sub domains */
class TSubDomainJoint : public TJointEqN
{
  protected:
    /** number of subdomains contain this joint */
    int N_NeibSubDomains;

    /** subdomain number of the neib cell, which contains this joint */
    int NeibSubDomainRank;

    /** global cell number of the neibs' cell, which contains this joint */
    int NeibSubDomainGlobalCellNo;

    /** local joint number of this joint in the neib cell */
    int NeibSubDomainLocalJointNo;

  public:
    // Constructors
    TSubDomainJoint(TBaseCell *neighb0);

    /** constructor with one initial neighbour */
    TSubDomainJoint(TBaseCell *neighb0, TBaseCell *neighb1);

     /** constructor with neighbour info */
    TSubDomainJoint(TBaseCell *neighb0, TBaseCell *neighb1, int neibID, 
                    int neibSubDomainGlobalCellNo,  int neibSubDomainLocalJointNo);

    // Methods
    /** return whether this is an interior joint */
    virtual bool InnerJoint()
    { return false; }

    /** return the subdomain number of the neib cell joint */
    int GetNeibRank()
    { return NeibSubDomainRank; }

    /** return the subdomain number of the neib cell joint */
    int GetNeibGlobalCellNo()
    { return NeibSubDomainGlobalCellNo; }

    /** return the subdomain number of the neib cell joint */
    int GetNeibLocalJointNo()
    { return NeibSubDomainLocalJointNo; }


  };
#endif
