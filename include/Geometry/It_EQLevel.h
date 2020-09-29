// =======================================================================
// @(#)It_EQLevel.h        1.1 10/30/98
// 
// Class:       TIt_EQ
// Purpose:     iterator to produce a series of cells which lie
//              exactly on a special refinement level
//
// Author:      Volker Behns  01.09.97
//
// =======================================================================

#ifndef __IT_EQLEVEL__
#define __IT_EQLEVEL__

#include <It_Search.h>

/** iterator to produce a series of cells which lie
    exactly on a special refinement level */
class TIt_EQLevel : public TIt_Search
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
