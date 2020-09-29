// =======================================================================
// @(#)SubDomainHaloJoint.h        
// 
// Class:       SubDomainHaloJoint
// Purpose:     superclass for boundary edges in Halo calls
//
// Author:      Sashikumaar Ganesan  10.09.09
//
// History:     
//
// =======================================================================

#ifndef __SUBDOMAINHALOJOINT__
#define __SUBDOMAINHALOJOINT__

#include <stdlib.h>
#include <JointEqN.h>

/** boundary edges which are in Halo calls */
class TSubDomainHaloJoint : public TJointEqN
{

  public:
    // Constructors
    TSubDomainHaloJoint(TBaseCell *neighb0);

    /** constructor with one initial neighbour */
    TSubDomainHaloJoint(TBaseCell *neighb0, TBaseCell *neighb1);

};
#endif
