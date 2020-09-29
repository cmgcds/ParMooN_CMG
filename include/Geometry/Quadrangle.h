// =======================================================================
// @(#)Quadrangle.h        1.1 10/30/98
//
// Class:       TQuadrangle
// Purpose:     shape descriptor of a quadrangle
//
// Author:      Volker Behns  16.07.97
//
// =======================================================================

#ifndef __QUADRANGLE__
#define __QUADRANGLE__

#include <ShapeDesc.h>

#define QUADMAXN_EpV  2

/** shape descriptor of a quadrangle */
class TQuadrangle : public TShapeDesc
{
  public:
    // Constructor
    /** initialize all description paramters */
    TQuadrangle();

    // Methods
    /** return diameter of a cell */
    virtual double GetDiameter(TVertex **Verts);

    /** return shortest edge of a cell */
    virtual double GetShortestEdge(TVertex **Verts);

    /** return the length of the cell defined with the reference map */
    virtual double GetLengthWithReferenceMap(TVertex **Verts);

    /** return measure of a cell */
    virtual double GetMeasure(TVertex **Verts);

    /** check a special quadrangle whether it is a parallelogram
        or even a square */
    Shapes CheckQuad(TVertex **Vertices);
};

#endif
