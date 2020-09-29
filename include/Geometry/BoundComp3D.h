// =======================================================================
// @(#)BoundComp3D.h        1.5 08/12/99
//
// Class:       TBoundComp3D
// Purpose:     components of boundary faces
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#ifndef __BOUNDCOMP3D__
#define __BOUNDCOMP3D__

#include <BoundComp.h>

/** components of boundary faces */
class TBoundComp3D : public TBoundComp
{
  public:
    // Constructor
    TBoundComp3D(int id);

    // Methods
    /** return the coordinates {X, Y, Z} of parameter values T and S*/
    virtual int GetXYZofTS(double T, double S, double &X, double &Y,
                          double &Z) = 0;
    /** return the parameter values T and S of coordinates (X, Y, Z) */
    virtual int GetTSofXYZ(double X, double Y, double Z, double &T,
                           double &S) = 0;

    /** return parameters and coordinates of a given linear
        combination of vertices */
    virtual int GetXYZandTS(int N_Points, double *LinComb,
                            double *xp, double *yp, double *zp,
                            double *tp, double *sp,
                            double &X, double &Y, double &Z,
                            double &T, double &S) = 0;

};

#endif
