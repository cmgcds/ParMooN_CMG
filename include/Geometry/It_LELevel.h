// =======================================================================
// @(#)It_LELevel.h        1.1 10/30/98
// 
// Class:       TIt_LELevel
// Purpose:     iterator to produce a series of cells which lie
//              on a special refinement level or on top of cell tree
//
// Author:      Volker Behns  01.09.97
//
// =======================================================================

#ifndef __IT_LELEVEL__
#define __IT_LELEVEL__

#include <It_Search.h>

/** iterator to produce a series of cells which lie on a special
    refinement level or on top of cell tree */
class TIt_LELevel : public TIt_Search
{
  protected:

  public:
    // Constructors

    // Methods
    /** return the next cell */
    virtual TBaseCell *Next(int &info);
    /** return the previous cell */
    virtual TBaseCell *Prev();
};

#endif
