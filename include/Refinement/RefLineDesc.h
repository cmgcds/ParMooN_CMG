// =======================================================================
// @(#)RefLineDesc.h        1.1 10/30/98
//
// Class:       TRefLineDesc
// Purpose:     refinement descriptor for line
//
// Author:      Volker Behns  23.07.97
//
// =======================================================================

#ifndef __REFLINEDESC__
#define __REFLINEDESC__

#include <RefDesc.h>

#define LINEMAXN_VpC    2
#define LINEMAXN_CpV    2
#define LINEMAXN_V      2

/** refinement descriptor for line */
class TRefLineDesc : public TRefDesc
{
  protected:

  public:
    // Constructor
    /** build a descriptor for regular refinement of a line */
    TRefLineDesc(TShapeDesc *shape);

    // Methods
};

#endif
