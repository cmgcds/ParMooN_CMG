// =======================================================================
// @(#)RefMortarLineDesc.h        1.1 10/30/98
//
// Class:       TRefMortarLineDesc
// Purpose:     refinement descriptor for a mortar line
//
// Author:      Volker Behns  16.01.98
//
// =======================================================================

#ifndef __REFMORTARLINEDESC__
#define __REFMORTARLINEDESC__

#include <RefDesc.h>

#define MLINEMAXN_VpC    2
#define MLINEMAXN_CpV    2

/** refinement descriptor for a mortar line */
class TRefMortarLineDesc : public TRefDesc
{
  protected:

  public:
    // Constructor
    /** build a descriptor for refinement of a mortar line */
    TRefMortarLineDesc(TShapeDesc *shape, int N);

    // Methods
};

#endif
