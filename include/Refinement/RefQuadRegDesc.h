// =======================================================================
// @(#)RefQuadRegDesc.h        1.1 10/30/98
//
// Class:       TRefQuadRegDesc
// Purpose:     refinement descriptor for regular refinement of a quadrangle
//
// Author:      Volker Behns  18.07.97
//
// =======================================================================

#ifndef __REFQUADREGDESC__
#define __REFQUADREGDESC__

#include <RefDesc.h>

#define QUADRRN_E         4
#define QUADRRMAXN_VpC    4
#define QUADRRMAXN_CpV    4
#define QUADRRMAXN_EpC    4
#define QUADRRMAXN_CpE    2
#define QUADRRMAXN_EpV    4
#define QUADN_V           4
#define QUADRRMAXN_iVpE   1
#define QUADRRMAXN_nVpoE  3
#define QUADRRMAXN_nEpoE  2

/** refinement descriptor for regular refinement of a quadrangle */
class TRefQuadRegDesc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for regular refinement of a quadrangle */
    TRefQuadRegDesc(TShapeDesc *shape);

    // Methods
};

#endif
