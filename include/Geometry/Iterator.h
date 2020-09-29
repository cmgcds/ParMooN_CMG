// =======================================================================
// @(#)Iterator.h        1.2 02/08/99
// 
// Class:       TIterator
// Purpose:     iterator to produce a series of cells which lye
//              exactly on a special level
//
// Author:      Volker Behns  04.08.97
//
// =======================================================================

#ifndef __ITERATOR__
#define __ITERATOR__

#define MAX_ItLevel  50
#define N_ITERATORS   9

enum Iterators {It_EQ, It_LE, It_Finest, It_EQLevel, It_LELevel,
                It_Between, It_OCAF, It_Mortar1, It_Mortar2};

#include <Domain.h>

/** iterator to produce a series of cells with some
    special properties */
class TIterator
{
  protected:
    /** current level */
    int Level;

    /** domain on which the iterator works */
    TDomain *Domain;
    /** a copy of tree of cells */
    TBaseCell **CellTree;
    /** number of cells on trees root */
    int N_RootCells;

  public:
    // Constructors

    // Methods
    /** set all parameters to the given values */
    int SetParam(TDomain *domain);

    /** return the next cell */
    virtual TBaseCell *Next(int &info) = 0;
    /** return the previous cell */
    virtual TBaseCell *Prev() = 0;

    /** Initialize at level "level" */
    virtual int Init(int level) = 0;

    /** return the maximum level */
    virtual int GetMaxLevel() = 0;
};

#endif
