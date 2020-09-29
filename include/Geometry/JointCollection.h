// =======================================================================
// @(#)JointCollection.h        4.1
//
// Class:       TJointCollection
// Purpose:     store Joint in an array
//              used by DG matrices assembler
//
// Author:      Sashikumaar Ganesan  03.11.09
//
// History:     03.11.09 Starting implementation
//
// ======================================================================= 

#ifndef __JOINTCOLLECTION__
#define __JOINTCOLLECTION__

#include <Joint.h>

/** store joints in an array, used by DG matrices assembler */
class TJointCollection
{
  protected:
    /** number of joints stored */
    int N_Joints;

    /** array containing the pointers to the cells */
    TJoint **Joints;

  public:
    /** constructor */
    TJointCollection(int n_joints, TJoint **joints);

    /** return number of joints */
    int GetN_Joints()
    { return N_Joints; }

    /** return joint with index i in joint-array */
    TJoint *GetJoint(int i)
    { return Joints[i]; }

    /** destructor: delete arrays */
    ~TJointCollection();
};

#endif
