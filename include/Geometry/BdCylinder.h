// =======================================================================
// @(#)BdSphere.h        1.1 07/16/99
//
// Class:       TBdSphere
// Superclass:  TBoundComp3D
// Purpose:     a Cylinder as a component of a boundary part
//
// Author:      Andreas Hahn 16.04.2010
//
// =======================================================================

#ifndef __BDCYLINDER__
#define __BDCYLINDER__

#include <BoundComp3D.h>

class TBdCylinder : public TBoundComp3D
{
  protected:
    /** radius **/
    double mRadius;
    
    /** coordinates of one point on axis **/
    double mPx, mPy, mPz;
    
    /** direction of axis **/
    double mAx, mAy, mAz;
    
    /** **/
    double mNx, mNy, mNz;
    
    /** **/
    double mBx, mBy, mBz;
    
  protected:
    
  public:
    // CTOR
    TBdCylinder (int id);
    
    virtual ~TBdCylinder () {};
    
    // Methods
    /** return the coordinates {X, Y, Z} of parameter values T and S */
    virtual int GetXYZofTS(double T, double S, double &X, double &Y,
                          double &Z);
    
     /** return the parameter values T and S of coordinates (X, Y, Z) */
    virtual int GetTSofXYZ(double X, double Y, double Z, double &T,
                           double &S);
			   
    /** return parameters and coordinates of a given linear
        combination of vertices */
    virtual int GetXYZandTS(int N_Points, double *LinComb,
                            double *xp, double *yp, double *zp,
                            double *tp, double *sp,
                            double &X, double &Y, double &Z,
                            double &T, double &S);
			  
    virtual int ReadIn(std::ifstream &dat);
    
    void SetParams(double r, double px, double py, double pz,
		   double ax, double ay, double az, double nx, double ny, double nz);
};

#endif