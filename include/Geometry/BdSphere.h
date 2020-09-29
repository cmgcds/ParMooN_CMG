// =======================================================================
// @(#)BdSphere.h        1.1 07/16/99
//
// Class:       TBdSphere
// Superclass:  TBoundComp3D
// Purpose:     a Sphere as a component of a boundary part
//
// Author:      Gunar Matthies 2000/12/04
//
// =======================================================================

#ifndef __BDSPHERE__
#define __BDSPHERE__

#include <BoundComp3D.h>

/** a Sphere as a component of a boundary part */
class TBdSphere : public TBoundComp3D
{
  protected:
      /** coordinates of mid point of the Sphere */
      double M_x, M_y, M_z;
      /** Radius */
      double R;

  public:
    // Constructor
    TBdSphere(int id);

    // Methods
    /** set all parameters to the given values */
    void SetParams (double m_x, double m_y, double m_z,
                    double r);

    /** return the coordinates of parameter value T, S */
    virtual int GetXYZofTS(double T, double S,
                           double &X, double &Y, double &Z);

    /** return the parameter value T, S of coordinates */
    virtual int GetTSofXYZ(double X, double Y, double Z,
                           double &T, double &S);

    /** return parameters and coordinates of a given linear
        combination of vertices */
    virtual int GetXYZandTS(int N_Points, double *LinComb,
                            double *xp, double *yp, double *zp,
                            double *tp, double *sp,
                            double &X, double &Y, double &Z,
                            double &T, double &S);

    /** read parameter from input file */
    virtual int ReadIn(std::ifstream &dat);
};

#endif
