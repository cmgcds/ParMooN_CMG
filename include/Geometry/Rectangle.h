// =======================================================================
// @(#)Rectangle.h        1.1 10/30/98
//
// Class:       TRectangle
// Purpose:     shape descriptor of a quadrangle, especially a rectangle
//
// Author:      Volker Behns  27.03.98
//
// =======================================================================

#ifndef __RECTANGLE__
#define __RECTANGLE__

#include <Quadrangle.h>

/** shape descriptor of a quadrangle, especially a rectangle */
class TRectangle : public TQuadrangle
{
  public:
    // Constructor
    /** initialize all description paramters */
    TRectangle();

    // Methods
    /** return diameter of a cell */
    virtual double GetDiameter(TVertex **Verts);

    /** return shortest edge of a cell */
    virtual double GetShortestEdge(TVertex **Verts);

    /** return the length of the cell defined with the reference map */
    virtual double GetLengthWithReferenceMap(TVertex **Verts);

    /** return measure of a cell */
    virtual double GetMeasure(TVertex **Verts);

};

#endif
