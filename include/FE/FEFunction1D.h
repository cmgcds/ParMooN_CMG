// =======================================================================
// @(#)FEFunction1D.h
// 
// Class:       TFEFunction1D
// Purpose:     a function from a finite element space in 1D
// 
// Author:      Sashikumaar Ganesan (17.05.2007)
//
// History:     start of implementation 17.05.2007
//
// =======================================================================

#ifndef __FEFUNCTION1D__
#define __FEFUNCTION1D__

#include <AllClasses.h>
#include <FESpace1D.h>
#include <Constants.h>

/** a function from a finite element space */
class TFEFunction1D
{
  protected:
    /** name of the function */
    char *Name;

    /** some more words describing the function */
    char *Description;

    /** space to which this function belongs to */
    TFESpace1D *FESpace1D;

    /** double vector according to FE isomorphism */
    double *Values;

    /** length of vector */
    int Length;

  public:
    /** constructor with vector initialization */
    TFEFunction1D(TFESpace1D *fespace1D, char *name, char *description,
                  double *values, int length);

    /** destructor */
    ~TFEFunction1D();

    /** return name */
    char *GetName()
    { return Name; }

    /** return description */
    char *GetDescription()
    { return Description; }

    /** return fe space */
    TFESpace1D *GetFESpace1D()
    { return FESpace1D; }

    /** return length */
    int GetLength()
    { return Length; }

    /** return vector of data */
    double *GetValues()
    { return Values; }


    /** calculate the interpolation of an exact function */
    void Interpolate(DoubleFunct2D *Exact);

    /** calculate the interpolation of an exact function */
    void Interpolate(int ConstCoord, double x, DoubleFunct2D *Exact);

   /** calculate for 1d function for given 2D , PBE*/
    void Interpolate(double x, double y, DoubleFunct3D *Exact);

    void InterpolateNodalPts(int N_Coord, double *Coords, DoubleFunctND *Exact, double *val);

    /** convert current grid to vector-values FE function */
    void GridToData();

};

#endif
