// =======================================================================
// @(#)It_Between.C        1.1 10/30/98
// 
// Class:       TIt_Between
// Purpose:     iterator to produce a series of cells whose level lies
//              between UNIFORM_STEPS and LEVELS

//
// Author:      Volker Behns  01.10.98
//
// =======================================================================

#include <It_Between.h>
#include <Database.h>

// Methods
int TIt_Between::Init(int level)
{
  TIt_Search::Init(level);

  Level = TDatabase::ParamDB->LEVELS;
  Level_bottom = TDatabase::ParamDB->UNIFORM_STEPS;

  return 0;
}

TBaseCell *TIt_Between::Next(int &info)
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
    { return NULL;}
  }
  
  while (Status[ActiveLevel].N_Children && ActiveLevel != Level)
  {
    ActiveCell = ActiveCell->GetChild(Status[ActiveLevel].CurrentChild++);
    Status[++ActiveLevel].N_Children = ActiveCell->GetN_Children();
  }

  info = ActiveLevel == Level ? 1 : ActiveLevel <= Level_bottom ? -1 : 0;
  return ActiveCell;
}

TBaseCell *TIt_Between::Prev()
{
  return 0;
}
