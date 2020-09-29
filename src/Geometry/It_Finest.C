// =======================================================================
// @(#)It_Finest.C        1.1 10/30/98
// 
// Class:       TIt_Finest
// Purpose:     iterator to produce a series of cells which lie
//              on top of cell tree

//
// Author:      Volker Behns  08.08.97
//
// =======================================================================

#include <It_Finest.h>

// Methods
int TIt_Finest::Init(int level)
{
  TIt_Search::Init(level);

  Level = MAX_ItLevel;

  return 0;
}

int TIt_Finest::GetMaxLevel()
{
  int MaxLevel = 0;

  do
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
        break;
    }

    while (Status[ActiveLevel].N_Children)
    {
      ActiveCell = ActiveCell->GetChild(Status[ActiveLevel].CurrentChild++);
      Status[++ActiveLevel].N_Children = ActiveCell->GetN_Children();

      if (ActiveLevel > MaxLevel) MaxLevel = ActiveLevel;
    }

  } while (TRUE);

  return MaxLevel;
}
