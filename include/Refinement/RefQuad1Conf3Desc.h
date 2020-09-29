// =======================================================================
// @(#)RefQuad1Conf3Desc.h        1.1 08/26/99
//
// Class:       TRefQuad1Conf3Desc
// Purpose:     refinement descriptor for conforming closure of a
//              quadrangle with 1 hanging node
//
// Author:      Matthias Ebeling  26.08.99
//
// =======================================================================

#ifndef __REFQUAD1CONF3DESC__
#define __REFQUAD1CONF3DESC__

#include <RefDesc.h>

#define QUADConfN_E          4
#define QUAD1ConfMAXN_VpC    3
#define QUAD1ConfMAXN_CpV    3
#define QUAD1ConfMAXN_EpC    3
#define QUAD1ConfMAXN_CpE    2
#define QUAD1ConfMAXN_EpV    4
#define QUADN_V              4
#define QUAD1ConfMAXN_iVpE   1
#define QUAD1ConfMAXN_nVpoE  3
#define QUAD1ConfMAXN_nEpoE  2

/** refinement descriptor for conforming closure of a quadrangle */
class TRefQuad1Conf3Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for conforming closure of a quadrangle */
    TRefQuad1Conf3Desc(TShapeDesc *shape);

    // Methods
};

#endif
