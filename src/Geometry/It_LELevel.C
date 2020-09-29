// =======================================================================
// @(#)It_LELevel.C        1.1 10/30/98
// 
// Class:       TIt_LELevel
// Purpose:     iterator to produce a series of cells which lie
//              on a special refinement level or on top of cell tree
//
// Author:      Volker Behns  01.09.97
//
// =======================================================================

#include <It_LELevel.h>

// Methods
TBaseCell *TIt_LELevel::Next(int &info)
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
  
  info = ActiveLevel;
  return ActiveCell;
}

TBaseCell *TIt_LELevel::Prev()
{
  return 0;
}
