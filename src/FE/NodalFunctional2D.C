// =======================================================================
// @(#)NodalFunctional2D.C        1.1 10/30/98
// 
// Class:       TNodalFunctional2D
// Purpose:     realize nodal functionals in 2D
//
// Author:      Gunar Matthies (02.09.98)
//
// History:     start of implementation 02.09.98 (Gunar Matthies)
//
// =======================================================================

#include <NodalFunctional2D.h>

/** constructor */
TNodalFunctional2D::TNodalFunctional2D(NodalFunctional2D id,
                int n_allfunctionals, 
                int n_edgefunctionals,
                int n_pointsall, int n_pointsedge,
                double *xi, double *eta, double *t,
                EvalAllNF *evalall,
                EvalJointNF *evaledge)
{
  ID = id;
  N_AllFunctionals = n_allfunctionals;
  N_EdgeFunctionals = n_edgefunctionals;
  N_PointsAll = n_pointsall;
  N_PointsEdge = n_pointsedge;
  Xi = xi;
  Eta = eta;
  T = t;
  EvalAll = evalall;
  EvalEdge = evaledge;
}
