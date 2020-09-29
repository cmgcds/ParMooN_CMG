// =======================================================================
// @(#)BdLine.h        1.2 07/16/99
//
// Class:       TBdLine
// Superclass:  TBoundComp
// Purpose:     a line as a component of a boundary part
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#ifndef __BDLINE__
#define __BDLINE__

#include <BoundComp2D.h>

/** a line as a component of a boundary part */
class TBdLine : public TBoundComp2D
{
  protected:
      /** x coordinate of the begin of line */
      double Xstart;
      /** y coordinate of the begin of line */
      double Ystart;
      /** x progress of line */
      double delX;
      /** y progress of line */
      double delY;

  public:
    // Constructor
    TBdLine(int id);

    // Methods
    /** set all parameters to the given values */
    void SetParams (double xstart, double ystart, double delx, double dely);

    /** return the coordinates of parameter value T */
    virtual int GetXYofT(double T, double &X, double &Y);

    /** return the parameter value T of coordinates (X, Y) */
    virtual int GetTofXY(double X, double Y, double &T);

    /** read parameter from input file */
    virtual int ReadIn(std::ifstream &dat);

    /** get number of initial vertices on this component */
    virtual int GetN_InitVerts()
    { return 2; }
    virtual int GenInitVerts(double *&points, int I_points,
                             int *&edges, int I_edges)
    { return 0; }
};

#endif
