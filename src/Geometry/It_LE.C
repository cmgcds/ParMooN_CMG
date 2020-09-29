// =======================================================================
// @(#)It_LE.C        1.1 10/30/98
// 
// Class:       TIt_LE
// Purpose:     iterator to produce a series of cells which lie
//              on a special level or on top of cell tree (with a
//              lower level)
//
// Author:      Volker Behns  04.08.97
//
// =======================================================================

#include <It_LE.h>

// Methods
TBaseCell *TIt_LE::Next(int &info)
{
  TBaseCell *ret;

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

  if (!(ActiveCell->ExistChildren()))
    Status[ActiveLevel].N_Children = Status[ActiveLevel].CurrentChild = 0;

  while (Status[ActiveLevel].N_Children && ActiveLevel != Level)
  {
    ActiveCell = ActiveCell->GetChild(Status[ActiveLevel].CurrentChild++);
    Status[++ActiveLevel].N_Children = ActiveCell->GetN_Children();
  }
  
  info = ActiveLevel;
  ret = ActiveCell;

  if (ActiveLevel)
    do
    {
      ActiveCell = ActiveCell->GetParent();
      Status[ActiveLevel].N_Children = 0;
      Status[ActiveLevel--].CurrentChild = 0;
    } while (ActiveLevel &&
        Status[ActiveLevel].N_Children == Status[ActiveLevel].CurrentChild);

  return ret;
}

TBaseCell *TIt_LE::Prev()
{
  return 0;
}
