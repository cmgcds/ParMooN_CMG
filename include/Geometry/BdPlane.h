// =======================================================================
// @(#)BdPlane.h        1.1 07/16/99
//
// Class:       TBdPlane
// Superclass:  TBoundComp
// Purpose:     a plane as a component of a boundary part
//
// Author:      Volker Behns  05.07.99
//
// =======================================================================

#ifndef __BDPLANE__
#define __BDPLANE__

#include <BoundComp3D.h>

/** a plane as a component of a boundary part */
class TBdPlane : public TBoundComp3D
{
  protected:
      /** coordinates of one point in the plane */
      double P_x, P_y, P_z;
      /** orthogonal vectors in the plain (used for parametization) */
      double A_x, A_y, A_z, B_x, B_y, B_z;
      /** outer normal Vector */
      double N_x, N_y, N_z;

  public:
    // Constructor
    TBdPlane(int id);
    
    virtual ~TBdPlane () {};

    // Methods
    /** set all parameters to the given values */
    void SetParams (double p_x, double p_y, double p_z,
                    double a_x, double a_y, double a_z,
                    double n_x, double n_y, double n_z);

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
    
    void GetParams (double &p_x, double &p_y, double &p_z,
                    double &a_x, double &a_y, double &a_z,
                    double &n_x, double &n_y, double &n_z);
};

#endif
