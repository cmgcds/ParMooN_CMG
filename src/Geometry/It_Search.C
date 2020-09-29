// =======================================================================
// @(#)It_Search.C        1.1 10/30/98
// 
// Class:       TIt_Search
// Purpose:     class of iterators which search the cell tree by
//              graph-theoretic methods
//
// Author:      Volker Behns  01.09.97
//
// =======================================================================

#include <It_Search.h>

// Methods
int TIt_Search::Init(int level)
{
  int i;

  Level = level;
  ActiveLevel = 0;
  ActiveRootCell = 0;

  for (i=0;i<MAX_ItLevel;i++)
    Status[i].N_Children = Status[i].CurrentChild = 0;

  return 0;
}

int TIt_Search::GetMaxLevel()
{
  return -1;
}
