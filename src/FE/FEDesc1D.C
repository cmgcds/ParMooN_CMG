// =======================================================================
// @(#)FEDesc1D.C        1.1 10/30/98
//
// Class:       TFEDesc1D
// Purpose:     store a finite element descriptor for a 1D element
//
// Author:      Gunar Matthies  23.07.98
//
// =======================================================================

#include <FEDesc1D.h>
#include <string.h>

/** constructor, setting all data */
TFEDesc1D::TFEDesc1D(char *description, int n_dof, int n_jointdof,
                     int **jointdof, int n_innerdof, int *innerdof)
{
  Description = strdup(description);
  N_DOF = n_dof;
  N_JointDOF = n_jointdof;
  JointDOF = jointdof;
  N_InnerDOF = n_innerdof;
  InnerDOF = innerdof;
}
