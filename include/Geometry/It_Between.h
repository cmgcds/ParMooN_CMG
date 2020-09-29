// =======================================================================
// @(#)It_Between.h        1.1 10/30/98
// 
// Class:       TIt_Between
// Purpose:     iterator to produce a series of cells which can be
//              refined or derefined
//
// Author:      Volker Behns  30.09.98
//
// =======================================================================

#ifndef __IT_BETWEEN__
#define __IT_BETWEEN__

#include <It_Search.h>

/** iterator to produce a series of cells which can be
    refined or derefined */
class TIt_Between : public TIt_Search
{
  protected:
    /** bottom level - do not return cells under this level */
    int Level_bottom;

  public:
    // Constructors

    // Methods
    /** initialize the iterator */
    virtual int Init(int level);

    /** return the next cell */
    virtual TBaseCell *Next(int &info);
    /** return the previous cell */
    virtual TBaseCell *Prev();
};

#endif
