// =======================================================================
// @(#)ADISystem.C        4.1 13.11.09
// 
// Class:       TADISYSTEM
// Purpose:     general super class for all ADI/operator-splitting System
//              special spaces are implemented in subclasses
//
// Author:      Sashikumaar Ganesan (13.11.09)
//
// History:     start of implementation 13.11.09 (Sashikumaar Ganesan)
//
// =======================================================================

#include <ADISystem.h>

TADISystem::TADISystem(double *Sol, double *OldSol, double *b, double *Defect, DoubleFunctND *growthAndNuc)
{
  sol = Sol;
  oldsol = OldSol;
  rhs = NULL;
  B = b;
  defect = Defect;
  GetGrowthAndNuc = growthAndNuc;  
}

TADISystem::~TADISystem()
{

}