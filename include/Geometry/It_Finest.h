// =======================================================================
// @(#)It_Finest.h        1.1 10/30/98
// 
// Class:       TIt_Finest
// Purpose:     iterator to produce a series of cells which lie
//              on top of cell tree
//
// Author:      Volker Behns  05.08.97
//
// =======================================================================

#ifndef __IT_Finest__
#define __IT_Finest__

#include <It_LE.h>

/** iterator to produce a series of cells which lie on
    top of cell tree */
class TIt_Finest : public TIt_LE
{
  protected:

  public:
    // Constructors

    // Methods
    /** initialize the iterator */
    virtual int Init(int level);

    /** return the maximum level */
    virtual int GetMaxLevel();
};

#endif
