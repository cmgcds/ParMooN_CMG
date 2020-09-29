// =======================================================================
// @(#)It_EQLevel.C        1.1 10/30/98
// 
// Class:       TIt_EQ
// Purpose:     iterator to produce a series of cells which lie
//              exactly on a special refinement level
//
// Author:      Volker Behns  01.09.97
//
// =======================================================================

#include <It_EQLevel.h>

// Methods
TBaseCell *TIt_EQLevel::Next(int &info)
{
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

    if (!Level  || (!ActiveLevel && 
          Status[ActiveLevel].N_Children == Status[ActiveLevel].CurrentChild))
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

    while (Status[ActiveLevel].N_Children && ActiveLevel != Level)
    {
      ActiveCell = ActiveCell->GetChild(Status[ActiveLevel].CurrentChild++);
      Status[++ActiveLevel].N_Children = ActiveCell->GetN_Children();
    }
  
  } while (ActiveLevel != Level);

  info = ActiveLevel;
  return ActiveCell;
}

TBaseCell *TIt_EQLevel::Prev()
{
  return 0;
}
