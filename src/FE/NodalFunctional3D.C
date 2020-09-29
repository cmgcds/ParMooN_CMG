// =======================================================================
// %W% %G% 
// 
// Class:       TNodalFunctional3D
// Purpose:     realize nodal functionals in 3D
//
// Author:      Gunar Matthies (30.05.2000)
//
// History:     start of implementation 30.05.00 (Gunar Matthies)
//
// =======================================================================

#include <NodalFunctional3D.h>

/** constructor */
TNodalFunctional3D::TNodalFunctional3D(NodalFunctional3D id,
                   int n_allfunctionals, int *n_facefunctionals,
                   int n_pointsall, int *n_pointsface,
                   double *xi, double *eta, double *zeta,
                   double **xiarray, double **etaarray,
                   double **zetaarray,
                   double *t, double *s,
                   EvalAllNF *evalall,
                   EvalJointNF *evalface)
{
  ID = id;

  N_AllFunctionals = n_allfunctionals;
  N_FaceFunctionals = n_facefunctionals;
  N_PointsAll = n_pointsall;
  N_PointsFace = n_pointsface;

  Xi = xi;
  Eta = eta;
  Zeta = zeta;

  XiArray = xiarray;
  EtaArray = etaarray;
  ZetaArray = zetaarray;

  T = t;
  S = s;

  EvalAll = evalall;
  EvalFace = evalface;
}

/** return information for points for all functionals */
void TNodalFunctional3D::GetPointsForAll(int &n_points, double* &xi,
                        double* &eta, double* &zeta)
{ 
  n_points = N_PointsAll;
  xi = Xi; 
  eta = Eta;
  zeta = Zeta;
}

/** return information for points for face functionals 
    on joint j */
void TNodalFunctional3D::GetPointsForFace(int j, int &n_points,
                        double* &xi, double* &eta, double* &zeta)
{ 
  n_points = N_PointsFace[j];
  xi   = XiArray[j];
  eta  = EtaArray[j];
  zeta = ZetaArray[j];
}

/** return information for points for face functionals */
void TNodalFunctional3D::GetPointsForFace(int &n_points, double* &t,
                                          double* &s)
{
  n_points = N_PointsFace[0],
  t = T;
  s = S;
}

