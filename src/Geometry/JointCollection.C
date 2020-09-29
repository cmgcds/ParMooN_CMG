// =======================================================================
// @(#)JointCollection.C        4.1
//
// Class:       TJointCollection
// Purpose:     store Joint in an array
//              used by DG matrices assembler
//
// Author:      Sashikumaa Ganesan  03.11.09
//
// History:     03.11.09 Starting implementation
//
// =======================================================================  

#include <JointCollection.h>
#include <Joint.h>

/** constructor */
TJointCollection::TJointCollection(int n_joints, TJoint **joints)
{
  N_Joints = n_joints;
  Joints = joints;
}


/** destructor: delete arrays */
TJointCollection::~TJointCollection()
{
 int i;

 if(N_Joints)
 {
  for(i=0; i<N_Joints; i++)
   {
    delete Joints[i];
   }
   delete [] Joints;
 }


}

