// =======================================================================
// @(#)It_EQ.h        1.1 10/30/98
// 
// Class:       TIt_EQ
// Purpose:     iterator to produce a series of cells which lie
//              exactly on a special level
//
// Author:      Volker Behns  05.08.97
//
// =======================================================================

#ifndef __IT_EQ__
#define __IT_EQ__

#include <It_Search.h>

/** iterator to produce a series of cells which lie
    exactly on a special level */
class TIt_EQ : public TIt_Search
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
