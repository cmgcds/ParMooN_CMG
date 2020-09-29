// =======================================================================
// @(#)NodalFunctional1D.C        1.1 11/20/98
// 
// Class:       TNodalFunctional1D
// Purpose:     realize nodal functionals in 1D
//
// Author:      Gunar Matthies (08.10.98)
//
// History:     start of implementation 08.10.98 (Gunar Matthies)
//
// =======================================================================

#include <NodalFunctional1D.h>

/** constructor */
TNodalFunctional1D::TNodalFunctional1D(NodalFunctional1D id,
                int n_allfunctionals, 
                int n_edgefunctionals,
                int n_pointsall, int n_pointsedge,
                double *xi, double *eta, double *t,
                DoubleFunctVect *evalall,
                DoubleFunctVect *evaledge)
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
