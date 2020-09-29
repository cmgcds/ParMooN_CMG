// =======================================================================
// @(#)RefQuad2Conf0Desc.h        1.1 08/26/99
//
// Class:       TRefQuad2Conf0Desc
// Purpose:     refinement descriptor for conforming closure of a
//              quadrangle with 2 hanging node
//
// Author:      Matthias Ebeling  20.08.99
//
// =======================================================================

#ifndef __REFQUAD2CONF0DESC__
#define __REFQUAD2CONF0DESC__

#include <RefDesc.h>

#define QUADConfN_E          4
#define QUAD2ConfMAXN_VpC    4
#define QUAD2ConfMAXN_CpV    3
#define QUAD2ConfMAXN_EpC    4
#define QUAD2ConfMAXN_CpE    2
#define QUAD2ConfMAXN_EpV    3
#define QUADN_V              4
#define QUAD2ConfMAXN_iVpE   1
#define QUAD2ConfMAXN_nVpoE  3
#define QUAD2ConfMAXN_nEpoE  2

/** refinement descriptor for conforming closure of a quadrangle */
class TRefQuad2Conf0Desc : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for conforming closure of a quadrangle */
    TRefQuad2Conf0Desc(TShapeDesc *shape);

    // Methods
};

#endif
