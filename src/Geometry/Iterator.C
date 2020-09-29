// =======================================================================
// @(#)Iterator.C        1.1 10/30/98
// 
// Class:       TIterator
// Purpose:     iterator to produce a series of cells with some
//              special properties
//
// Author:      Volker Behns  04.08.97
//
// =======================================================================

#include <Iterator.h>

// Constructors

// Methods
int TIterator::SetParam(TDomain *domain)
{
  Domain = domain;

  Domain->GetTreeInfo(CellTree, N_RootCells);

  return 0;
}
