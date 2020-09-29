// =======================================================================
// @(#)It_Search.h        1.1 10/30/98
// 
// Class:       TIt_Search
// Purpose:     class of iterators which search the cell tree by 
//              graph-theoretic methods
//
// Author:      Volker Behns  01.09.97
//
// =======================================================================

#ifndef __IT_SEARCH__
#define __IT_SEARCH__

#include <Iterator.h>

/** class of iterators which search the cell tree by
    graph-theoretic methods */
class TIt_Search : public TIterator
{
  protected:
    /** status vector */
    struct {int N_Children, CurrentChild;} Status[MAX_ItLevel];

    /** active level in status vector */
    int ActiveLevel;
    /** active root cell */
    int ActiveRootCell;

    /** active cell */
    TBaseCell *ActiveCell;

  public:
    // Constructors

    // Methods
    /** Initialize on level */
    virtual int Init(int level);

    /** return the maximum level */
    virtual int GetMaxLevel();
};

#endif
