// =======================================================================
// @(#)BdWall.C        1.2 07/16/99
//
// Class:       TBdWall
// Purpose:     a sandwich side wall as a component of a boundary part
//
// Author:      Volker Behns  05.07.99
//
// =======================================================================

#include <BdWall.h>
#include <MooNMD_Io.h>

// Constructor
TBdWall::TBdWall(int id, TBoundComp2D *bdcomp2d) : TBoundComp3D(id)
{
  Type = Wall;

  BdComp2D = bdcomp2d;
}

// Methods
int TBdWall::SetParams(double drx, double dry, double drz)
{
  DriftX = drx;
  DriftY = dry;
  DriftZ = drz;

  return 0;
}

int TBdWall::GetXYZofTS(double T, double S,
                        double &X, double &Y, double &Z)
{
  double x,y,z;

  BdComp2D->GetXYofT(T, x, y);

  X = x + S * DriftX;
  Y = y + S * DriftY;
  Z = S * DriftZ;

  return 0;
}

/** return parameters and coordinates of a given linear
    combination of vertices */
int TBdWall::GetXYZandTS(int N_Points, double *LinComb,
                         double *xp, double *yp, double *zp,
                         double *tp, double *sp,
                         double &X, double &Y, double &Z,
                         double &T, double &S)
{
  int i;
  double t, s, v;
  double x,y,z;

  t = s = 0;
  for(i=0;i<N_Points;i++)
  {
    v = LinComb[i];
    t += v*tp[i];
    s += v*sp[i];
  }
  T = t;
  S = s;

  BdComp2D->GetXYofT(T, x, y);

  X = x + S * DriftX;
  Y = y + S * DriftY;
  Z = S * DriftZ;

  return 0;
}

int TBdWall::GetTSofXYZ(double X, double Y, double Z,
                        double &T, double &S)
{
  double x,y,t;

  S = Z / DriftZ;

  x = X - S * DriftX;
  y = Y - S * DriftY;

  BdComp2D->GetTofXY(x, y, t);

  T = t;
  
  return 0;
}

int TBdWall::ReadIn(std::ifstream &dat)
{
  return BdComp2D->ReadIn(dat);
}
