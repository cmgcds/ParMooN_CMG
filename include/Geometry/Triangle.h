// =======================================================================
// @(#)Triangle.h        1.1 10/30/98
//
// Class:       TTriangle
// Purpose:     shape descriptor of a triangle
//
// Author:      Volker Behns  16.07.97
//
// =======================================================================

#ifndef __TRIANGLE__
#define __TRIANGLE__

#include <ShapeDesc.h>

#define TRIMAXN_EpV  2

/** shape descriptor of a triangle */
class TTriangle : public TShapeDesc
{
  public:
    // Constructor
    /** build the shape descriptor for a triangle */
    TTriangle();

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
