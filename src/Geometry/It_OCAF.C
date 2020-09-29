// =======================================================================
// @(#)It_OCAF.C        1.1 10/30/98
// 
// Class:       TIt_OCAF
// Purpose:     iterator to produce a series of cells which lie
//              one level under finest one

//
// Author:      Volker Behns  01.10.98
//
// =======================================================================

#include <It_OCAF.h>

// Methods
int TIt_OCAF::Init(int level)
{
  TIt_Search::Init(level);

  Level = MAX_ItLevel;

  return 0;
}

TBaseCell *TIt_OCAF::Next(int &info)
{
  bool OK = false;

  while (!OK)
  {
    if (ActiveLevel)
      do
      {
        ActiveCell = ActiveCell->GetParent();
        Status[ActiveLevel].N_Children = 0;
        Status[ActiveLevel--].CurrentChild = 0;
      } while (ActiveLevel &&
          Status[ActiveLevel].N_Children == Status[ActiveLevel].CurrentChild);

    if (!ActiveLevel &&
         Status[ActiveLevel].N_Children == Status[ActiveLevel].CurrentChild)
    {
      if (ActiveRootCell < N_RootCells)
      {
        ActiveCell = CellTree[ActiveRootCell++];
        Status[ActiveLevel].N_Children = ActiveCell->GetN_Children();
        Status[ActiveLevel].CurrentChild = 0;
      }
      else
        return NULL;
    }

    while (Status[ActiveLevel].N_Children)
    {
      ActiveCell = ActiveCell->GetChild(Status[ActiveLevel].CurrentChild++);
      Status[++ActiveLevel].N_Children = ActiveCell->GetN_Children();
    }
  
    if (ActiveLevel)
    {
      ActiveCell = ActiveCell->GetParent();
      Status[ActiveLevel].N_Children = 0;
      Status[ActiveLevel].CurrentChild = 0;
      Status[--ActiveLevel].CurrentChild = ActiveCell->GetN_Children();
      OK = true;
    }
  }

  info = ActiveLevel;
  return ActiveCell;
}

TBaseCell *TIt_OCAF::Prev()
{
  return 0;
}
