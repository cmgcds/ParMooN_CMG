// =======================================================================
// @(#)Joint.h        
// 
// Class:       TSubDomainFace
// Purpose:     superclass for faces between two subdomains
//
// Author:      Sashikumaar Ganesan  07.09.09
//
// History:     
//
// =======================================================================

#ifdef __3D__

#ifndef __SUBDOMAINFACE__
#define __SUBDOMAINFACE__

#include <stdlib.h>
#include <JointEqN.h>

/** for setting the neib cell information while refinement */
struct SubDomainMap
{
    int ParentLocalCell_No;
    int NeibParentLocalCell_No;
    int ParentLocalEdge_No;
    int NeibParentLocalEdge_No;
    int NewEdge;
    int NeibNewEdge;
    bool MapFilled;
};

/** connects two cells, which are in two/(more in 3D) different sub domains */
class TSubDomainFace : public TJointEqN
{

#endif
#endif