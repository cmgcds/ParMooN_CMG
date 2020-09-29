// =======================================================================
// @(#)BoundComp2D.h        1.5 08/12/99
//
// Class:       TBoundComp2D
// Purpose:     Components of boundary faces
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#ifndef __BOUNDCOMP2D__
#define __BOUNDCOMP2D__

#include <BoundComp.h>

/** Components of boundary faces */
class TBoundComp2D : public TBoundComp
{
  public:
    // Constructor
    TBoundComp2D(int id);

    // Methods
    /** return the coordinates {X,Y} of parameter value T */
    virtual int GetXYofT(double T, double &X, double &Y) = 0;
    /** return the parameter value T of coordinates (X, Y) */
    virtual int GetTofXY(double X, double Y, double &T) = 0;

    /** get number of initial vertices on a Comp2Donent */
    virtual int GetN_InitVerts() = 0;
    virtual int GenInitVerts(double *&points, int I_points,
                             int *&edges, int I_edges) = 0;

};

#endif
