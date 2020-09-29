// =======================================================================
// @(#)QuadFormula3D.C        1.2 05/04/99
//
// Class:      TQuadFormula3D
// Superclass: TQuadFormula
//
// Purpose:    quadrature formula for a 3D integral
// Author:     Gunar Matthies
//
// History:    30.08.1997 start implementation
// 
// =======================================================================

#include <QuadFormula3D.h>
#include <FEDatabase3D.h>

TQuadFormula3D::TQuadFormula3D() : TQuadFormula()
{
  Xi=NULL;
  Eta=NULL;
  Zeta=NULL;
}

TQuadFormula3D::TQuadFormula3D(int n_points, double* weights, 
                     double* xi, double* eta, double* zeta, int acc)
{
  InitObject(n_points, weights, xi, eta, zeta, acc);
}

void TQuadFormula3D::InitObject(int n, double* w, 
                      double* xi, double* eta, double* zeta, int acc)
{
  int i;

  N_QuadPoints=n;

  Weights=new double[n];
  Xi=new double[n];
  Eta=new double[n];
  Zeta=new double[n];

  for(i=0;i<n;i++)
  {
    Weights[i]=w[i];
    Xi[i]=xi[i];
    Eta[i]=eta[i];
    Zeta[i]=zeta[i];
  }

  Accuracy = acc;
}

double *TQuadFormula3D::GetCoords(int i)
{
  double *ret=NULL;

  switch(i)
  {
    case 0: ret=Xi;   break;
    case 1: ret=Eta;  break;
    case 2: ret=Zeta; break;
  }
    
  return ret;
}

void TQuadFormula3D::
  GetFormulaData(int &n_points, double* &weights, 
                 double* &xi, double* &eta, double* &zeta)
{
  n_points=N_QuadPoints;
  weights=Weights;
  xi=Xi;
  eta=Eta;
  zeta=Zeta;
}

#ifdef __3D__
void TQuadFormula3D::FindLocalQuadFormula3D
        (int N_LocalUsedElements, FE3D *LocalUsedElements,
         QuadFormula2D &qf1, QuadFormula3D &qf2)
{
  int i,j, MaxPolynomialDegree, PolynomialDegree;
  BF3DRefElements RefElement;
  int *PolynomialDegreeFromFE3D;
  RefTrans3D RefTrans, *RefTransArray, CurrentRefTrans;
  TRefTrans3D *rt;
  BaseFunct3D BaseFuncts[N_FEs3D];

  // find adequate quadrature formula for all elements
  // and find needed reference transformation
  RefTransArray = TFEDatabase3D::GetRefTrans3D_IDFromFE3D();
  RefTrans = RefTransArray[LocalUsedElements[0]];
  PolynomialDegreeFromFE3D = 
        TFEDatabase3D::GetPolynomialDegreeFromFE3D();
  MaxPolynomialDegree = 0;
  RefElement = 
        TFEDatabase3D::GetRefElementFromFE3D(LocalUsedElements[0]);
  for(i=0;i<N_LocalUsedElements;i++)
  {
    BaseFuncts[i] = 
        TFEDatabase3D::GetBaseFunct3D_IDFromFE3D(LocalUsedElements[i]);
    PolynomialDegree = PolynomialDegreeFromFE3D[LocalUsedElements[i]];
    if(PolynomialDegree > MaxPolynomialDegree) 
      MaxPolynomialDegree = PolynomialDegree;

    CurrentRefTrans = RefTransArray[LocalUsedElements[i]];
    if(CurrentRefTrans > RefTrans)
      RefTrans = CurrentRefTrans;
  }

  switch(RefElement)
  {
    case BFUnitHexahedron:
      qf2 = TFEDatabase3D::GetQFHexaFromDegree
                                        (2*MaxPolynomialDegree);
    break;

    case BFUnitTetrahedron:
      qf2 = TFEDatabase3D::GetQFTetraFromDegree
                                        (2*MaxPolynomialDegree);
    break;
  } // endswitch

  // qf1 = TFEDatabase3D::GetQFLineFromDegree
  //                                       (2*MaxPolynomialDegree);
}
#endif // __3D__

std::ostream& operator << (std::ostream &s, TQuadFormula3D *qf)
{
  int i,N_;

  s << endl;
  s << "instance of TQuadFormula3D" << endl;
  N_=qf->N_QuadPoints;
  s << "N_QuadPoints: " << N_ << endl;
  s << setw(3) << "No";
  s << setw(11) << "weights";
  s << setw(11) << "xi";
  s << setw(11) << "eta";
  s << setw(11) << "zeta";
  s << endl;
  s.setf(std::ios::fixed);
  for(i=0;i<N_;i++)
  {
    s << setw(3) << i;
    s << setprecision(6) << setw(11) << qf->Weights[i];
    s << setprecision(6) << setw(11) << qf->Xi[i];
    s << setprecision(6) << setw(11) << qf->Eta[i];
    s << setprecision(6) << setw(11) << qf->Zeta[i];
    s << endl;
  }

  return s << endl;
}
