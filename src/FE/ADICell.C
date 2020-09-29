// =======================================================================
// @(#)ADICell.C        4.1 07.11.09
// 
// Class:       src for TADICell
// Purpose:     general super class for all ADICells
//              special spaces are implemented in subclasses
//
// Author:      Sashikumaar Ganesan (07.11.09)
//
// History:     start of implementation 07.11.09 (Sashikumaar Ganesan)
//
// =======================================================================

#include <ADICell.h>
#include <Database.h>

TADICell::TADICell(TBaseCell *cell, int N_quadPts)
{
  Cell = cell;
  N_QuadPts = N_quadPts;

  sol = NULL;
  solT = NULL;
  rhs = NULL;
  rhsT = NULL;
}

TADICell::~TADICell()
{

}