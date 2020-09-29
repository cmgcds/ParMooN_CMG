// =======================================================================
// @(#)FE1D.C        1.1 10/30/98
//
// Class:       TFE1D
// Purpose:     store all information for one finite element class
//
// Author:      Gunar Matthies  08.10.98
//
// =======================================================================

#include <FE1D.h>
#include <FEDatabase2D.h>

/** constructor */
TFE1D::TFE1D()
{
}

/** constructor with data */
TFE1D::TFE1D(BaseFunct1D basefunct_id, NodalFunctional1D nodalfunctional_id,
         RefTrans1D reftransid, FEDesc1D fedesc_id, int n_info)
{
  BaseFunct_ID = basefunct_id;
  BaseFunct = TFEDatabase2D::GetBaseFunct1D(BaseFunct_ID);

  NodalFunctional_ID = nodalfunctional_id;
  NodalFunctional  = TFEDatabase2D::GetNodalFunctional1D(NodalFunctional_ID);

  RefTransID = reftransid;

  FEDesc_ID = fedesc_id;
  FEDesc = TFEDatabase2D::GetFEDesc1D(FEDesc_ID);

  N_Info = n_info;
  N_DOF = BaseFunct->GetDimension();

  Size = N_Info + N_DOF;
}

