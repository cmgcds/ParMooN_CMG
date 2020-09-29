// =======================================================================
// @(#)RefNoRef.h        1.1 10/30/98
//
// Class:       TRefNoRef
// Purpose:     no refinement - only path to shape descriptor
//
// Author:      Volker Behns  21.07.97
//
// =======================================================================

#ifndef __REFNOREF__
#define __REFNOREF__

#include <RefDesc.h>

/** no refinement - only path to shape descriptor */
class TRefNoRef : public TRefDesc
{
  public:
    // Constructor
    /** build a descriptor for a nonrefined object */
    TRefNoRef(TShapeDesc *shape);

    // Methods
    /** return TRUE because object is to refine */
    virtual int IsToRefine()
    { return FALSE; }

};

#endif
