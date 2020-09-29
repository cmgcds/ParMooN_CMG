// =======================================================================
// @(#)FEVectFunct2D.h        
// 
// Class:       TFEVectFunct1D
// Purpose:     a function from a finite element space in 1D
//
// Author:      Sashikumaar Ganesan (26.09.09)
//
// History:     start of implementation 26.09.09 (Sashikumaar Ganesan)
//
//
// =======================================================================

#ifndef __FEVECTFUNCT1D__
#define __FEVECTFUNCT1D__

#include <FEFunction1D.h>

/** a function from a finite element space */
class TFEVectFunct1D : public TFEFunction1D
{
  protected:
    /** number of components */
    int N_Components;

  public:
    /** constructor with vector initialization */
    TFEVectFunct1D(TFESpace1D *fespace1D, char *name, char *description,
                  double *values, int length, int n_components);

    /** return number of components */
    int GetN_Components()
    { return N_Components; }

    /** return i-th component as FEFunction2D */
    TFEFunction1D *GetComponent(int i)
    {
      return new TFEFunction1D(FESpace1D, Name, Description,
                               Values+i*Length, Length);
    }

    /** convert current grid to vector-values FE function */
//     void GridToData();

    /** use current data for grid replacement */
//     void DataToGrid();



};

#endif
