// =======================================================================
// @(#)FEDesc2D.C        1.1 10/30/98
//
// Class:       TFEDesc2D
// Purpose:     store a finite element descriptor for a 2D element
//
// Author:      Gunar Matthies  23.07.98
//
// =======================================================================

#include <FEDesc2D.h>
#include <string.h> // in order to use strdup

/** constructor, setting all data with dof on cell boundary */
TFEDesc2D::TFEDesc2D(char *description, int n_dof, int n_jointdof,
                     int **jointdof, int n_innerdof, int *innerdof,
                     int n_outerdof, int *outerdof)
{
  Description = strdup(description);
  N_DOF = n_dof;
  N_JointDOF = n_jointdof;
  JointDOF = jointdof;
  N_InnerDOF = n_innerdof;
  InnerDOF = innerdof;
  N_OuterDOF = n_outerdof;
  OuterDOF = outerdof;
}

/** constructor, setting all data with dof on cell boundary */
TFEDesc2D::TFEDesc2D(char *description, int n_dof, int n_jointdof,
                     int **jointdof, int n_innerdof, int *innerdof)
{
  Description = strdup(description);
  N_DOF = n_dof;
  N_JointDOF = n_jointdof;
  JointDOF = jointdof;
  N_InnerDOF = n_innerdof;
  InnerDOF = innerdof;
  N_OuterDOF = 0;
  OuterDOF = NULL;
}

/** return joint on which the i-th local degree of freedom is   
If i is not a dof on an edge, return -1

If i is a dof on two edges (i.e. on a vertex), one of these two edges is 
returned. Don't use this function in this case.

*/
int TFEDesc2D::GetJointOfThisDOF(int localDOF) const
{
  int i,j;
  bool is_DOF_on_edge=false;
 
  
  for (i=0;i<N_OuterDOF;i++)
  {
    if(OuterDOF[i]==localDOF)
    {
      is_DOF_on_edge=true;
      break;
    }
  }
  if(!is_DOF_on_edge)
    return -1;
  //else // continue to find the edge
  i=0;
  while (true)
  {
    // this must terminate, since we already know localDOF is a dof on an edge
    for (j=0;j<N_JointDOF;j++)
    {
      if(JointDOF[i][j] == localDOF)
        return i;
    }
    i++;
  }
}
