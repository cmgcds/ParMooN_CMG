// =======================================================================
// @(#)It_LE.h        1.1 10/30/98
// 
// Class:       TIt_LE
// Purpose:     iterator to produce a series of cells which lie
//              on a special level or on top of cell tree (with a
//              lower level)
//
// Author:      Volker Behns  05.08.97
//
// =======================================================================

#ifndef __IT_LE__
#define __IT_LE__

#include <It_Search.h>

/** iterator to produce a series of cells which lie on a special
    level or on top of cell tree (with a lower level) */
class TIt_LE : public TIt_Search
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
