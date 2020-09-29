// =======================================================================
// @(#)Hexahedron.h        1.2 10/18/99
//
// Class:       THexahedron
// Purpose:     shape descriptor of a hexahedron
//
// Author:      Volker Behns  16.07.97
//
// =======================================================================

#ifndef __HEXAHEDRON__
#define __HEXAHEDRON__

#include <ShapeDesc.h>

#define HEXAMAXN_EpV   3
#define HEXAMAXN_VpF   4
#define HEXAMAXN_FpV   3
#define HEXAMAXN_EpF   4
#define HEXAMAXN_FpE   2

/** shape descriptor of a hexahedron */
class THexahedron : public TShapeDesc
{
  public:
    // Constructor
    THexahedron();

    // Methods
    /** return diameter of a cell */
    virtual double GetDiameter(TVertex **Verts);

    /** return shortest edge of a cell */
    virtual double GetShortestEdge(TVertex **Verts);

    /** return the length of the cell defined with the reference map */
    virtual double GetLengthWithReferenceMap(TVertex **Verts);

    /** return measure of a cell */
    virtual double GetMeasure(TVertex **Verts);

    /** check whether this hexahedron is a brick */
    Shapes CheckHexa(TVertex **Vertices);
};

#endif
