// =======================================================================
// @(#)Parallelogram.C        1.1 10/30/98
//
// Class:       TParallelogram
// Purpose:     shape descriptor of a quadrangle, especially a
//              parallelogram
//
// Author:      Volker Behns  06.02.98
//
// =======================================================================

#include <Parallelogram.h>

// Constructor
TParallelogram::TParallelogram() : TQuadrangle()
{
  Type = Parallelogram;
}

// Methods
