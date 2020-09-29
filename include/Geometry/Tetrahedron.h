// =======================================================================
// @(#)Tetrahedron.h        1.2 10/18/99
//
// Class:       TTetrahedron
// Purpose:     shape descriptor of a tetrahedron
//
// Author:      Volker Behns  16.07.97
//
// =======================================================================

#ifndef __TETRAHEDRON__
#define __TETRAHEDRON__

#include <ShapeDesc.h>

#define TETRAMAXN_EpV   3
#define TETRAMAXN_VpF   3
#define TETRAMAXN_FpV   3
#define TETRAMAXN_EpF   3
#define TETRAMAXN_FpE   2

/** shape descriptor of a tetrahedron */
class TTetrahedron : public TShapeDesc
{
  public:
    // Constructor
    /** build the shape descriptor for a tetrahedron */
    TTetrahedron();

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
