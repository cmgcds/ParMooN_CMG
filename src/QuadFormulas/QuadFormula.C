// =======================================================================
// @(#)QuadFormula.C        1.2 05/04/99
//
// Class:    TQuadFormula
//
// Purpose:  base class for quadrature formulas
// Author:   Gunar Matthies
//
// History:  29.08.1997 start implementation
//           07.04.1999 add Accuracy
// 
// =======================================================================

#include <QuadFormula.h>

TQuadFormula::TQuadFormula()
{
  N_QuadPoints=0;
  Weights=NULL;

  Accuracy = -1;
}

double *TQuadFormula::GetCoords(int i)
{
  return NULL;
}

std::ostream& operator << (std::ostream &s, TQuadFormula *qf)
{
  int i,N_;

  s << endl;
  s << "instance of TQuadFormula" << endl;
  N_=qf->N_QuadPoints;
  s << "N_QuadPoints: " << N_ << endl;
  s << setw(3) << "No";
  s << setw(11) << "weights";
  s << endl;
  s.setf(std::ios::fixed);
  for(i=0;i<N_;i++)
  {
    s << setw(3) << i;
    s << setprecision(6) << setw(11) << qf->Weights[i];
    s << endl;
  }

  return s << endl;
}
