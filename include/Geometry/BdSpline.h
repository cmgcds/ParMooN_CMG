// =======================================================================
// @(#)BdSpline.h        1.2 07/16/99
//
// Class:       TBdSpline
// Superclass:  TBoundComp
// Purpose:     splie function as a component of a boundary part
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#ifndef __BDSPLINE__
#define __BDSPLINE__

#include <BoundComp2D.h>

/** splie function as a component of a boundary part */
class TBdSpline : public TBoundComp2D
{
  protected:
    /** number of subsplines */
    int N_Splines;
    /** array for all parameters */
    double *Params;

  public:
    // Constuctor
    /** constructor initializes the parameter array */
    TBdSpline (int id, int N_Splines);

    // Methods
    /** set all parameters */
    void SetParams (double *params);
    /** get number of splines */
    int GetN_Splines();
    /** return new inner points for a good boundary approximation */
    //double *GetBoundPoints(QuadFormula *formula, int &N_NewPt);

    /** return the coordinates of parameter value T */
    virtual int GetXYofT(double T, double &X, double &Y);

    /** return the parameter value T of coordinates (X, Y) */
    virtual int GetTofXY(double X, double Y, double &T);

    /** read parameter from input file */
    virtual int ReadIn(std::ifstream &dat);

    /** get number of initial vertices on this component */
    virtual int GetN_InitVerts()
    { return 4; }
    virtual int GenInitVerts(double *&points, int I_points,
                             int *&edges, int I_edges)
    { return -1; }
};

#endif
