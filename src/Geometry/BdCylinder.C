// =======================================================================
// @(#)BdSphere.h        1.1 07/16/99
//
// Class:       TBdSphere
// Superclass:  TBoundComp3D
// Purpose:     a Cylinder as a component of a boundary part
//
// Author:      Andreas Hahn 23.04.2010
//
// =======================================================================

#include <BdCylinder.h>

// CTOR
TBdCylinder::TBdCylinder (int id) : TBoundComp3D(id)
{  
  Type = Cylinder;
}

/** Methods **/

int TBdCylinder::ReadIn(std::ifstream &dat)
{
  return 0;
}

/** return the coordinates {X, Y, Z} of parameter values T and S*/
int TBdCylinder::GetXYZofTS(double T, double S, double &X, double &Y,
                            double &Z)
{
  double alpha = mRadius*cos(T);
  double beta  = mRadius*sin(T);
  
  
  X = mPx + alpha*mNx + beta*mBx + S*mAx;
  Y = mPy + alpha*mNy + beta*mBy + S*mAy;
  Z = mPz + alpha*mNz + beta*mBz + S*mAz;
  
  return 0;
}

/** return the parameter values T and S of coordinates (X, Y, Z) */
int TBdCylinder::GetTSofXYZ(double X, double Y, double Z, double &T,
                            double &S)
{
  double x, y, z, a, b;
  
  x = X - mPx;
  y = Y - mPy;
  z = Z - mPz;
  
  S = x*mAx + y*mAy + z*mAz;
  
  a = x*mNx + y*mNy + z*mNz;
  b = x*mBx + y*mBy + z*mBz;
  
  T = atan2(b, a);
  
  return 0;
}

/** return parameters and coordinates of a given linear
        combination of vertices */
int TBdCylinder::GetXYZandTS(int N_Points, double *LinComb,
                            double *xp, double *yp, double *zp,
                            double *tp, double *sp,
                            double &X, double &Y, double &Z,
                            double &T, double &S)
{
  double tmp, x, y, z;
  
  X = Y = Z = 0;
  
  for (int i=0;i<N_Points;++i)
  {
    tmp = LinComb[i];
    
    X += xp[i]*tmp;
    Y += yp[i]*tmp;
    Z += zp[i]*tmp;
  }
  
  X -= mPx; Y -= mPy; Z -= mPz;
  
  x = X*mNx + Y*mNy + Z*mNz;
  y = X*mBx + Y*mBy + Z*mBz;
  z = X*mAx + Y*mAy + Z*mAz;
  
  tmp = sqrt(x*x + y*y)/mRadius;
  
  x /= tmp;
  y /= tmp;
  
  X = x*mNx + y*mBx + z*mAx + mPx;
  Y = x*mNy + y*mBy + z*mAy + mPy;
  Z = x*mNz + y*mBz + z*mAz + mPz;
  
  GetTSofXYZ ( X, Y, Z, T, S );
  
  return 0;
}

void TBdCylinder::SetParams(double r, double px, double py, double pz,
		   double ax, double ay, double az, double nx, double ny, double nz)
{
  double norm;
  
  mRadius = r;
  
  mPx = px; mPy = py; mPz = pz;
  mAx = ax; mAy = ay; mAz = az;
  mNx = nx; mNy = ny; mNz = nz;
  
  mBx = mAy*mNz - mAz*mNy;
  mBy = mAz*mNx - mAx*mNz;
  mBz = mAx*mNy - mAy*mNx;
  
  norm = sqrt(mAx*mAx + mAy*mAy + mAz*mAz);
  mAx /= norm; mAy /= norm; mAz /= norm;
  
  norm = sqrt(mNx*mNx + mNy*mNy + mNz*mNz);
  mNx /= norm; mNy /= norm; mNz /= norm;
  
  norm = sqrt(mBx*mBx + mBy*mBy + mBz*mBz);
  mBx /= norm; mBy /= norm; mBz /= norm;
}