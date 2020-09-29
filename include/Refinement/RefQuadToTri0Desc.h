// =======================================================================
// @(#)RefQuadToTri0Desc.h        1.1 10/30/98
//
// Class:       TRefQuadToTri0Desc
// Purpose:     refinement descriptor for refinement of a quadrangle
//              into two triangles
//
// Author:      Volker Behns  16.03.98
//
// =======================================================================

#ifndef __REFQUADTOTRI0DESC__
#define __REFQUADTOTRI0DESC__

#include <RefDesc.h>

#ifndef __REFQUADTOTRI1DESC__
  #define QUADTTN_E         3
  #define QUADTTMAXN_VpC    3
  #define QUADTTMAXN_CpV    2
  #define QUADTTMAXN_EpC    3
  #define QUADTTMAXN_CpE    2
  #define QUADTTMAXN_EpV    3
  #define QUADTTMAXN_nVpoE  2
  #define QUADTTMAXN_nEpoE  1
#endif

/** refinement descriptor for refinement of a quadrangle into
    two triangles */
class TRefQuadToTri0Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor */
    TRefQuadToTri0Desc(TShapeDesc *shape);

    // Methods
};

#endif
