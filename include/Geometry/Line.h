// =======================================================================
// @(#)Line.h        1.1 10/30/98
//
// Class:       TLine
// Purpose:     shape descriptor of a line
//
// Author:      Volker Behns  16.07.97
//
// =======================================================================

#ifndef __LINEH__
#define __LINEH__

#include <ShapeDesc.h>

#define LINEMAXN_EpV  1

/** shape descriptor of a line */
class TLine : public TShapeDesc
{
  protected:

  public:
    // Constructor
    /** initializes the data structure */
    TLine();

    // Methods
    /** return diameter of a cell */
    virtual double GetDiameter(TVertex **Verts)
    {
      return 0;
    }

    /** return shortest edge of a cell */
    virtual double GetShortestEdge(TVertex **Verts)
    {
      return 0;
    }  
    /** return the length of the cell defined with the reference map */
    virtual double GetLengthWithReferenceMap(TVertex **Verts)
    {
      return 0;
    }  

    /** return measure of a cell */
    virtual double GetMeasure(TVertex **Verts);
};

#endif
