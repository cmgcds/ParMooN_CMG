// =======================================================================
// @(#)ItMethod.C        1.24 06/27/00
//
// Class:       TItMethod
// Purpose:     iteration methods
//
// Author:      Volker John 24.10.2000
//
// History:     24.10.2000 start of implementation
//
// =======================================================================
#include <MooNMD_Io.h>
#include <Database.h>
#include <LinAlg.h>
#include <NSE_MultiGrid.h>
#include <ItMethod.h>

#include <stdlib.h>
/** constructor with initialization */
TItMethod::TItMethod(MatVecProc *MatVec, DefectProc *Defect, 
                     TItMethod *Prec,
                     int n_aux, int N_Unknowns)
{
}

TItMethod::~TItMethod()
{
  if (N_Aux>0)
    delete AuxArray[0];
}

