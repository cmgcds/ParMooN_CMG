// =======================================================================
// @(#)RefTrans1D.h        
//
// Class:      TRefTrans1D
//
// Purpose:    reference transformations for 1D geometric objects
//
// Author:     Sashikumaar Ganesan
//
// History:    17.05.2007 start implementation
// 
// =======================================================================

#ifndef __REFTRANS1D__
#define __REFTRANS1D__

#include <Constants.h>
#include <Enumerations.h>
#include <BaseCell.h>

/** reference transformations for 1D geometric objects */
class TRefTrans1D
{
  protected:
    TBaseCell *Cell;

  public:
    /** constuctor */
    TRefTrans1D() {};

    /** transfer form reference element to original element */
    void GetOrigFromRef(double xi, double &x);


    /** transfer form reference element to original element */
    void GetOrigFromRef(double *ref, double *orig);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double x, double &eta);

    /** transfer from original element to reference element */
    void GetRefFromOrig(double *orig, double *ref);

    /** calculate functions and derivatives from reference element
        to original element */
    void GetOrigValues(TBaseCell *cell);

    /** set original element to cell */
    virtual void SetCell(TBaseCell *cell)
    {  Cell = cell; }

    static RefTrans1D FindRefTrans1D
        (int N_LocalUsedElements, FE1D *LocalUsedElements);


    /** return volume of cell according to reference transformation */
    double GetVolume();
};

#endif

