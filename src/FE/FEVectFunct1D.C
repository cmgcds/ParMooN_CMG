// =======================================================================
// @(#)FEVectFunct2D.C
// 
// Class:       TFEVectFunct1D
// Purpose:     a function from a finite element space in 1D
//
// Author:      Sashikumaar Ganesan (26.09.09)
//
// History:     start of implementation 26.09.09 (Sashikumaar Ganesan)
//
//
// =======================================================================

#include <FEVectFunct1D.h>

/** constructor with vector initialization */
TFEVectFunct1D::TFEVectFunct1D(TFESpace1D *fespace1D, char *name, 
                             char *description, double *values, 
                             int length, int n_components)
  : TFEFunction1D(fespace1D, name, description, values, length)
{
  N_Components = n_components;
}

