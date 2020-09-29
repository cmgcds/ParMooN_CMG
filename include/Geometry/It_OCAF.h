// =======================================================================
// @(#)It_OCAF.h        1.1 10/30/98
// 
// Class:       TIt_OCAF OneCoarserAsFinest
// Purpose:     iterator to produce a series of cells which lie
//              one level under the finest one
//
// Author:      Volker Behns  30.09.98
//
// =======================================================================

#ifndef __IT_OCAF
#define __IT_OCAF

#include <It_Search.h>

/** iterator to produce a series of cells which lie
    one level under the finest one */
class TIt_OCAF : public TIt_Search
{
  protected:

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
