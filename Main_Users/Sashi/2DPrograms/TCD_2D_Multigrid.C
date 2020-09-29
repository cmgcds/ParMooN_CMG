// =======================================================================
//
// Purpose:     main program with multigrid solver 
//              more than once scalar equations and with a given velocity
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 29.08.2013
//        :     operator splitting first and second order on 29.08.2013
//        :     MPI implementation started on  
// ======================================================================= 
#include <Domain.h>
#include <Database.h>
#include <SquareStructure1D.h>
#include <SquareMatrix1D.h>
#include <FEFunction1D.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <AuxParam2D.h>
#include <QuadAffin.h>
#include <Solver.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <DiscreteForm2D.h>
#include <LinAlg.h>
#include <Collection.h>
#include <JointCollection.h>
// #include <TCD2D.h>
// #include <TimeConvDiff2D.h>
#include <ConvDiff2D_Routines.h>
#include <LocalProjection.h>
#include <BaseFunct1D.h>
#include <NodalFunctional1D.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TriaIsoparametric.h>

#include <BoundEdge.h>
#include <IsoBoundEdge.h>
#include <InterfaceJoint.h>
#include <ADISystem.h>
#include <ADISystem1D.h>
#include <LineAffin.h>

#include <ItMethod.h>
#include <FixedPointIte.h>
#include <MultiGridScaIte.h>
#include <FgmresIte.h>
#include <MultiGrid2D.h>
#include <MGLevel2D.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

double bound = 0;

#include <MainUtilities.h>
#include <TimeUtilities.h>
#include <FEM_TVD_FCT.h>

#ifdef _MPI
#include <MeshPartition.h>
#include <ParFECommunicator2D.h>
#include <MumpsSolver.h>
#include <ParVector.h>
#endif

#ifdef _OMP
#include <ParDirectSolver.h>
#endif

#define AMG 0
#define GMG 1
#define DIRECT 2

// =======================================================================
// include current example
// =======================================================================

// #include "../Examples/TCD_2D/SimPaTurS.h"
// #include "../Examples/TCD_2D/2D_1D_ADI_SimPaTurS.h"
// #include "../Examples/TCD_2D/2D_1D_ADI_SimPaTurS_Aggrigation.h"
// #include "../Examples/TCD_2D/2D_1D_ConstT_UnitSqr.h"
#include "../Examples/TCD_2D/2D_1D_Smmoth.h"
// #include "../Examples/TCD_2D/2D_1D_ConstT_UnitSqr_IntlOnly.h"
 
// ======================================================================
// utilities for main program
// ======================================================================


#ifdef __PBS__

void Sort(double *XArray, double *YArray, int length)
{
  int n=0, l=0, r=length-1, m;
  int i, j, k, *rr, len, s;
  double Mid, Temp;
  double lend = length;

  len=(int)(2*log(lend)/log((double) 2.0)+2);
  rr= new int[len];

  do
  {
    do
    {
      i=l;
      j=r;

      m=(l+r)/2;
      Mid=XArray[m];

      do
      {
        while( XArray[i] > Mid) i++;

        while(XArray[j] < Mid) j--;

        if (i<=j)
        {
          Temp=XArray[i];
          XArray[i]=XArray[j];
          XArray[j]=Temp;

          Temp=YArray[i];
          YArray[i]=YArray[j];
          YArray[j]=Temp;

          i++; j--;
        }
      } while (i<=j);

      if (l<j)
      {
        rr[++n]=r;
        r=j;
      }
    } while (l<j);

    if (n>0) r=rr[n--];

    if (i<r) l=i;

  } while (i<r);

  delete [] rr;

}


int GetIndex(double *x, double *y, int Length, double x0, double y0)
{
  int l=0, r=Length, m=(r+l)/2;
  double MidX, MidY;
  bool update;

  MidX=x[m];
  MidY=y[m];

  //   cout<<endl;
  //   cout << " x " << x0<< " y " << y0<<endl;
  //   cout << " MidX " << MidX<< " MidY " << MidY<<endl;
  //  if((fabs(MidX-x0)<1e-8 &&  fabs(MidY-y0)<1e-8))
  //     cout << "TRUE " <<endl;

  while(!( fabs(MidX-x0)<1.e-4 && fabs(MidY-y0)<1.e-4 ) )
  {

    if(fabs(MidX-x0) <1.e-4 )
    {
      if(MidY> y0)
        { l=m;}
        else
          { r=m; }
    }
    else if(MidX > x0)
      { l=m; }
      else
        { r=m; }

        m=(r+l)/2;
    MidX=x[m];
    MidY=y[m];
    //cout << m << " MidX " << MidX<< " MidY " << MidY<< " x0  " <<  x0 << "  y0 " <<  y0 <<endl;

  }
  //  cout << " outer m " << m<< endl;
  //   cout << r << " MidX " << MidX<< " MidY " << MidY<< " y[m-1] " <<  y[m-1]<< "  y[m+1] " <<  y[m+1]<<endl;
  return m;
}

void GetNodalPtsADI_System(int &N_IntlPts, TFESpace1D *FeSpace_Intl,
TFESpace2D *FESpace2D, int &MaxN_PtsForNodal, double *&IntlX, double *&IntlY)
{
  int i,j,k,l, N_AllLocalPoints, N_V;
  int N_Cells, PolynomialDegree;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points, *RootPtIndex;
  int *DOF, N_Edges, ApproxOrder, N_XCommonPts, N_RootPts, disp;

  double *xi, *eta, x0, y0;
  double *x_loc, *y_loc, *x_loc_origOrder, *y_loc_origOrder;
  double X[MaxN_PointsForNodal2D], Y[MaxN_PointsForNodal2D];
  double AbsDetjk[MaxN_PointsForNodal2D], lastx, lasty;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_PointsForNodal2D];

  boolean IsIsoparametric, update;

  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId;
  TFE2D *Element;
  TFE2D *FE_Obj;
  TNodalFunctional2D *nf;
  TBaseFunct2D *bf;
  RefTrans2D RefTrans, *RefTransArray;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;
  BF2DRefElements RefElement;
  TRefTrans2D *rt;
  RefTrans2D F_K;
  QuadFormula2D QuadFormula;
  TVertex *Vertices;

  
  Coll = FESpace2D->GetCollection();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();

  //   #ifdef _MPI
  //   N_Cells = Coll->GetN_OwnCells();
  //   #else
  N_Cells = Coll->GetN_Cells();
  //   #endif

  RefTransArray = TFEDatabase2D::GetRefTrans2D_IDFromFE2D();

  //first find the nimber of AllLocal points
  N_AllLocalPoints = 0;
  MaxN_PtsForNodal = -1;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    nf = Element->GetNodalFunctional2D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_AllLocalPoints +=N_Points;
    if(MaxN_PtsForNodal<N_Points) MaxN_PtsForNodal=N_Points;
  }

  x_loc = new double[N_AllLocalPoints];
  y_loc = new double[N_AllLocalPoints];
  x_loc_origOrder = new double[N_AllLocalPoints];
  y_loc_origOrder = new double[N_AllLocalPoints];

  N_AllLocalPoints = 0;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    nf = Element->GetNodalFunctional2D();
    nf->GetPointsForAll(N_Points, xi, eta);

    PolynomialDegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    ApproxOrder = TFEDatabase2D::GetAccuracyFromFE2D(FEId);
    RefElement = Element->GetBaseFunct2D()->GetRefElement();

    switch(RefElement)
    {
      case BFUnitSquare:
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*PolynomialDegree);
        N_Edges = 4;
        break;

      case BFUnitTriangle:
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(3*PolynomialDegree-1);
        N_Edges = 3;
        break;
    }

    RefTrans = RefTransArray[FEId];

    IsIsoparametric = FALSE;
    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(j=0;j<N_Edges;j++)
      {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryEdge)
        {
          bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
        if(jointtype == InterfaceJoint)
        {
          bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
        if(jointtype == IsoInterfaceJoint ||
          jointtype == IsoBoundEdge)
          IsIsoparametric = TRUE;
      }
    }                                             // endif

    if(IsIsoparametric)
    {
      switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          break;
      }
    }

    switch(RefTrans)
    {
      case QuadAffin:
        rt = TFEDatabase2D::GetRefTrans2D(QuadAffin);
        ((TQuadAffin *)rt)->SetCell(cell);
        F_K = QuadAffin;
        break;
      case QuadBilinear:
        rt = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
        ((TQuadBilinear *)rt)->SetCell(cell);
        F_K = QuadBilinear;
        break;
      case QuadIsoparametric:
        rt = TFEDatabase2D::GetRefTrans2D(QuadIsoparametric);
        ((TQuadIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TQuadIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TQuadIsoparametric *)rt)->SetCell(cell);
        F_K = QuadIsoparametric;
        break;
      case TriaAffin:
        rt = TFEDatabase2D::GetRefTrans2D(TriaAffin);
        ((TTriaAffin *)rt)->SetCell(cell);
        F_K = TriaAffin;
        break;
      case TriaIsoparametric:
        rt = TFEDatabase2D::GetRefTrans2D(TriaIsoparametric);
        ((TTriaIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TTriaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TTriaIsoparametric *)rt)->SetCell(cell);
        F_K = TriaIsoparametric;
        break;
      default:
        cout << "unknown reftrans id: " << RefTrans << endl;
    }

    TFEDatabase2D::GetOrigFromRef(F_K, N_Points, xi, eta,  X, Y, AbsDetjk);

    for(j=0; j<N_Points; j++)
    {
      x_loc[N_AllLocalPoints] = X[j];
      y_loc[N_AllLocalPoints] = Y[j];
      N_AllLocalPoints++;
    }
  }                                               // for(i=0; i<N_Cells; i++)

  memcpy(x_loc_origOrder, x_loc,  N_AllLocalPoints*SizeOfDouble);
  memcpy(y_loc_origOrder, y_loc,  N_AllLocalPoints*SizeOfDouble);

  Sort(x_loc, y_loc, N_AllLocalPoints);

  lastx  = x_loc[0];
  N_RootPts = 0;
  N_XCommonPts = 0;
  disp = 0;

  for(i=0; i<N_AllLocalPoints; i++)
  {
    if( fabs(x_loc[i]-lastx)>1e-5 )
    {
      xi = x_loc+disp;
      eta = y_loc+disp;
      Sort(eta, xi, N_XCommonPts);

      lasty= eta[0];
      N_RootPts++;
      for(j=0; j<N_XCommonPts; j++)
      {
        if( fabs(eta[j]-lasty)>1e-4)
        {
          N_RootPts++;
          lasty = eta[j];
        }
      }                                           //  for(j=0; j<N_XCommonPts; j++)

      disp +=N_XCommonPts;
      N_XCommonPts = 1;
      lastx = x_loc[i];
    }
    else
      { N_XCommonPts++; }
    }                                             // for(i=0; i<N_AllLocalPoints; i++)
    // sort the final part
    xi = x_loc+disp;
  eta = y_loc+disp;
  Sort(eta, xi, N_XCommonPts);
  lasty= eta[0];
  N_RootPts++;
  for(j=0; j<N_XCommonPts; j++)
  {
    if( fabs(eta[j]-lasty)>1e-4)
    {
      N_RootPts++;
      lasty = eta[j];
    }
  }                                               //  for(j=0; j<N_XCommonPts; j++)

  cout<< "N_RootPts " << N_RootPts <<endl;
  // exit(0);
  IntlX = new double[N_RootPts];
  IntlY = new double[N_RootPts];

  lastx  = x_loc[0];
  N_RootPts = 0;
  N_XCommonPts = 0;
  disp = 0;
  for(i=0; i<N_AllLocalPoints; i++)
  {
    if( fabs(x_loc[i]-lastx)>1e-4 )
    {
      xi = x_loc+disp;
      eta = y_loc+disp;
      lasty= eta[0];
      IntlX[N_RootPts] = xi[0];
      IntlY[N_RootPts] = eta[0];
      N_RootPts++;

      for(j=0; j<N_XCommonPts; j++)
      {
        if( fabs(eta[j]-lasty)>1e-4)
        {
          IntlX[N_RootPts] = xi[j];
          IntlY[N_RootPts] = eta[j];
          N_RootPts++;
          lasty = eta[j];
        }
      }                                           //  for(j=0; j<N_XCommonPts; j++)
      disp +=N_XCommonPts;
      N_XCommonPts = 1;
      lastx = x_loc[i];
    }
    else
      { N_XCommonPts++; }
    }                                             //for(i=0; i<N_AllLocalPoints; i++)

    xi = x_loc+disp;
  eta = y_loc+disp;
  lasty= eta[0];
  IntlX[N_RootPts] = xi[0];
  IntlY[N_RootPts] = eta[0];
  N_RootPts++;
  for(j=1; j<N_XCommonPts; j++)
  {
    if( fabs(eta[j]-lasty)>1e-4)
    {
      IntlX[N_RootPts] = xi[j];
      IntlY[N_RootPts] = eta[j];
      N_RootPts++;
      lasty = eta[j];
    }
  }                                               //  for(j=0; j<N_XCommonPts; j++)

  //   for(i=0; i<N_RootPts; i++)
  //    cout<< i << " IntlX " << IntlX[i]<< " IntlY " << IntlY[i]  <<endl;

  delete [] x_loc;
  delete [] y_loc;

  RootPtIndex = new int[N_AllLocalPoints];

  // find the index for the local points in the root points
  for(i=0; i<N_AllLocalPoints; i++)
  {
    x0 = x_loc_origOrder[i];
    y0 = y_loc_origOrder[i];

    RootPtIndex[i]=GetIndex(IntlX, IntlY,  N_RootPts,  x0, y0);
  }

  delete [] x_loc_origOrder;
  delete [] y_loc_origOrder;

  FESpace2D->SetIntlPtIndexOfPts(RootPtIndex);

  N_IntlPts = N_RootPts;

  //    for(i=0; i<N_IntlPts; i++)
  //        cout <<i<< " IntlX[i] " <<IntlX[i]<< " IntlY[i] " <<IntlY[i] <<endl;
  //
  // cout<< "N_RootPts " << N_RootPts <<endl;
}



void GetQuadPtsADI_System(int &N_IntlPts, TFESpace2D *FESpace2D, double *&IntlX, double *&IntlY)
{
  int i,j,k,l, m;
  int N_Cells, PolynomialDegree, N_Points;

  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], AbsDetjk[MaxN_QuadPoints_2D];
  double *xi, *eta, *weights;

  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId;
  TFE2D *Element;
  BF2DRefElements RefElement;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;
  TRefTrans2D *F_K;

  Coll = FESpace2D->GetCollection();

  //   #ifdef _MPI
  //   N_Cells = Coll->GetN_OwnCells();
  //   #else
  N_Cells = Coll->GetN_Cells();
  //   #endif

  m = 0;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
    PolynomialDegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);

    switch(RefElement)
    {
      case BFUnitSquare:
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*PolynomialDegree);
        F_K = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
        ((TQuadBilinear *)F_K)->SetCell(cell);
        break;

      case BFUnitTriangle:
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(3*PolynomialDegree-1);
        F_K = TFEDatabase2D::GetRefTrans2D(TriaAffin);
        ((TTriaAffin *)F_K)->SetCell(cell);
        break;
    }                                             // endswitch

    //cout << "QuadFormula: " << QuadFormula << endl;
    qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
    qf2->GetFormulaData(N_Points, weights, xi, eta);

    if(i==0)
    {
      N_IntlPts = N_Points*N_Cells;
      //     cout<< "N_IntlPts " << N_IntlPts << endl;
      IntlX = new double[N_IntlPts];
      IntlY = new double[N_IntlPts];
    }

    switch(RefElement)
    {
      case BFUnitSquare:
        ((TQuadBilinear *)F_K)->GetOrigFromRef(N_Points, xi, eta, X, Y, AbsDetjk);
        break;

      case BFUnitTriangle:
        ((TTriaAffin *)F_K)->GetOrigFromRef(N_Points, xi, eta, X, Y, AbsDetjk);
        break;
    }                                             // endswitch

    for(j=0; j<N_Points; j++)
    {
      IntlX[m] = X[j];
      IntlY[m] = Y[j];
      m++;
    }
  }                                               // for(i=0; i<N_Cells; i++)

  //  cout<< "N_IntlPts " << m <<endl;
  //   for(i=0; i<N_IntlPts; i++)
  //    cout<< i << " IntlX " << IntlX[i]<< " IntlY " << IntlY[i]  <<endl;
  //  exit(0);

}


 

void GetInternalNodalPts(TFESpace1D * FeSpace_Intl, int &N_Intl_Levels, double *&IntlPosL)
{
  int i, j, k, l, m, r, N_Cells, N_RootPts, *RootPtIndex;
  int N_DOFs, N_LocalDOFs, N_Points;
  int N_AllLocalPoints;

  double L, L0, *xi, *eta, *L_loc, *L_loc_origOrder;
  double X[MaxN_PointsForNodal1D], Y[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;

  Coll = FeSpace_Intl->GetCollection();
  N_Cells = Coll->GetN_Cells();

  N_AllLocalPoints = 0;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FeSpace_Intl->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    N_AllLocalPoints +=N_Points;
  }

  L_loc = new double [N_AllLocalPoints];
  L_loc_origOrder = new double [N_AllLocalPoints];
  N_AllLocalPoints = 0;
  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    FEId = FeSpace_Intl->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();

    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);

    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, X, Y, AbsDetjk);

    for(j=0; j<N_Points; j++)
    {
      L_loc[N_AllLocalPoints] = X[j];
      N_AllLocalPoints++;
    }
  } // for(i=0; i<N_Cells; i++)

  memcpy(L_loc_origOrder, L_loc,  N_AllLocalPoints*SizeOfDouble);

  for(i=0; i<N_AllLocalPoints-1; i++)
   for(j=i; j<N_AllLocalPoints; j++)
    if(L_loc[i]> L_loc[j])
     {
      L= L_loc[i];
      L_loc[i]= L_loc[j];
      L_loc[j]= L;
     }

  L  = L_loc[0];
  N_RootPts = 1;

  for(i=1; i<N_AllLocalPoints; i++)
   {
    if( fabs(L_loc[i]-L)>1e-5 )
     {
      N_RootPts++;
      L = L_loc[i];
     }
   }

  IntlPosL= new double[N_AllLocalPoints];
  IntlPosL[0] = L_loc[0];
  N_RootPts = 1;
  L  = L_loc[0];

  for(i=1; i<N_AllLocalPoints; i++)
   {
    if( fabs(L_loc[i]-L)>1e-5 )
     {
      IntlPosL[N_RootPts] = L_loc[i];
      N_RootPts++;
      L = L_loc[i];
     }
   }

  delete [] L_loc;

  RootPtIndex = new int[N_AllLocalPoints];

  // find the index for the local points in the root points
  for(i=0; i<N_AllLocalPoints; i++)
   {
    L = L_loc_origOrder[i];
    l=0;
    r=N_RootPts;

    m = N_RootPts/2;
    L0 = IntlPosL[m];

    while(fabs(L-L0) > 1.e-8 )
     {
      if(L < L0)  //poin lies in left side
       {
        r = m;
       }
      else
       {
        l=m;
       }

      m= (l+r)/2;
      L0 = IntlPosL[m];
     } //  while ( 

    RootPtIndex[i] = m;
   }

  FeSpace_Intl->SetIntlPtIndexOfPts(RootPtIndex);
  FeSpace_Intl->SetN_RootNodalPts(N_RootPts);

  N_Intl_Levels = N_RootPts;

//  cout << N_AllLocalPoints << " N_RootPts  "  << N_RootPts << endl;
//   for(i=0; i<N_RootPts; i++)
//   cout << i << " L: "  << IntlPosL[i] << endl;
// exit(0);
} // void GetInternalNodalPts(T




void GetPtsValuesForIntl(int N_Levels, TFEFunction2D **ScalarFunctions, double *Sol, int N_U, double *ValuesT, int N_IntlPts)
{
  int i, ii, j, k, l, m, N_Cells, N_Points, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers, *DOF;
  int PolynomialDegree, ApproxOrder, *IntlPtIndexOfPts, *PtIndexLoc, disp;
  int *IncidentArray, *Incident;

  double *xi, *eta, *NodalValues, s, *Values_Level;
  double BasisValues[MaxN_BaseFunctions2D], maxval=0.;

  TFESpace2D *PBE_Spaces;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId_PBE;
  TFE2D *Element_PBE;
  TBaseFunct2D *bf_PBE;
  TNodalFunctional2D *nf_PBE;

  // assume all  ScalarFunctions use same fespace2D
  PBE_Spaces = ScalarFunctions[0]->GetFESpace2D();
  BeginIndex = PBE_Spaces->GetBeginIndex();
  GlobalNumbers = PBE_Spaces->GetGlobalNumbers();

  IntlPtIndexOfPts = PBE_Spaces->GetIntlPtIndexOfPts();
  IncidentArray = new int[N_Levels*N_IntlPts];

  memset(ValuesT, 0, N_Levels*N_IntlPts*SizeOfDouble);
  memset(IncidentArray, 0, N_Levels*N_IntlPts*SizeOfInt);

  // assume that both fespace2D and PBE_Spaces use same coll
  Coll = PBE_Spaces->GetCollection();
  N_Cells = Coll->GetN_Cells();

  disp = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId_PBE = PBE_Spaces->GetFE2D(i, cell);
    Element_PBE = TFEDatabase2D::GetFE2D(FEId_PBE);
    nf_PBE = Element_PBE->GetNodalFunctional2D();
    nf_PBE->GetPointsForAll(N_Points, xi, eta);
    bf_PBE = Element_PBE->GetBaseFunct2D();
    N_LocalDOFs = Element_PBE->GetN_DOF();

    DOF = GlobalNumbers + BeginIndex[i];
    PtIndexLoc = IntlPtIndexOfPts + disp;

    for(j=0;j<N_Points;j++)
    {
      bf_PBE->GetDerivatives(D00, xi[j], eta[j], BasisValues);

      k = PtIndexLoc[j];

      Incident = IncidentArray + k*N_Levels;
      Values_Level  =  ValuesT + k*N_Levels;

      for(ii=0;ii<N_Levels;ii++)
      {
        Incident[ii]++;
        NodalValues =  Sol + ii*N_U;
        for(l=0;l<N_LocalDOFs;l++)
        {
          m = DOF[l];
          s = NodalValues[m];
          Values_Level[ii] += BasisValues[l]*s;
        }
      }                                           // for(ii=0;ii<N_Levels;ii++)
    }                                             // for(j=0;j<N_Points;j++)

    disp +=N_Points;
  }                                               // for(i=0;i<N_Cells;i++)

  for(i=0;i<N_IntlPts;i++)
  {
    Incident = IncidentArray + i*N_Levels;
    Values_Level  =  ValuesT + i*N_Levels;

    for(j=0;j<N_Levels;j++)
    {
      Values_Level[j] /= (double)Incident[j];
      if(maxval < Values_Level[j]) maxval = Values_Level[j];
    }
  }

  delete [] IncidentArray;

  // cout <<  " GetPtsValuesForIntl_Nodal " << maxval <<endl;
  //exit(0);
}



#endif


int main(int argc, char* argv[])
{
  // ======================================================================
  // variable declaration
  // ======================================================================
 
  TDomain *Domain = new TDomain();
  TDomain *Domain_Intl = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *coll, *Coll_Intl, *mortarcoll = NULL;
  TJointCollection *JointColl;
  TBaseCell *cell;
  TFESpace2D ***Scalar_Spaces, **Velocity_Spaces, *pressure_space, **Grid_space;
  TOutput2D *Output, *Output_PBS;
  TAuxParam2D *aux;
  BoundCondFunct2D *BoundaryConditions[4];
  BoundCondFunct2D *BDCond[4];
  DoubleFunct2D *InitiaValues[4];
  BoundValueFunct2D *BoundValues[4];
  BoundValueFunct2D *BDValue[4];
  CoeffFct2D *Coefficients[4];  
  MatVecProc *MatVect;
  DefectProc *Defect;
  TDiscreteForm2D *DiscreteForm;
  TDiscreteForm2D *DiscreteFormMatrixMRhs[4], *DiscreteFormMatrixMRhs_SUPG[4];
  TDiscreteForm2D *DiscreteFormMatrixARhs[4], *DiscreteFormMatrixARhs_SUPG[4];
  TSquareStructure2D ***sqstructureA;
  TSquareMatrix2D *sqmatrixA, *SQMATRICES[3], *SQMATRICESSOLVER[3];
  TSquareMatrix2D *sqmatrixM , *sqmatrixK, *sqmatrixS;
  TSquareMatrix2D ***MatricesA, ***MatricesM, ***MatricesK;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;  
  TFEVectFunct2D **velocity;
  TFEFunction2D ***ScalarFunctions, **velo1, **velo2, *fefct[3];
  TMGLevel2D *MGLevel, *MGLevel_low;
  TMultiGrid2D *MG;
  TItMethod *itmethod, *prec;    
  TFESpace1D *FeSpace_Intl;
  TSquareStructure1D *SqStruct1D_Intl;
  TSquareMatrix1D *M_Intl, *A_Intl, *S_Intl, *K_Intl;
  TADISystem1D **ADI_System;
  TFEFunction1D *FeFunction_Intl;
  BoundCond cond_Lmin, cond_Lmax;
  TBoundEdge *Outlet_Joint;
  TBoundComp *BoundComp;
  TJoint *Joint;
  TFEFunction2D **ScalarPbeFunctions, *PbeFunctions;
  
  double impuls_residual,linredfac, values[5];
  double gamma, tau, oldtau, maxval=0., L0, L1;
  double Parameters[2], hmin, hmax, limit;
  double ***RhsArray, ***SolArray, **B, **B_Pbe, **RhsArray_Pbe, **OldRhsArray_Pbe, **UArray;
  double *IntlPosL, *IntlX, *IntlY;;
  double *Scalar_Entries, l2, H1, *PBE_IntlPtValuesT, *RhsPbe_Intl_Loc, *OldPBE_IntlPtValuesT, *C_IntlPtValues;
  double *PBE_LSolT;
  double *PBE_IntlPtValues, *PBE_IntlPtRhsValues, *T_IntlPtValues, *velo_IntlPts, *grad_velo_IntlPts;
  double len, *C_Sat_IntlPtValues, *SolPbe, *SolPbe_OneLevel, *PBE_IntlPtRhsValuesT, *OldPBE_IntlPtRhsValuesT;
  double *SolPbe_Lmin, *SolPbe_Lmax, q3_max, T;
  double *Sol_IntlLoc, *OldSol_IntlLoc, *B_IntlLoc, *defect_IntlLoc, *Sol_Loc, *SolPbe_Output, *MatValues;
  double *WArray, *gridpos, *gridpos_old;

  
  int i, j, k, l, m, N, ret, img=1, N_Cells, N_Intl_Levels, N_U, N_LDof, N_U_PBE;  
  int ORDER, VELOCITYORDER, mg_type, mg_level, LEVELS, N_Paramters, zerostart, FirstSolve;
  int **N_Uarray, **N_Active, N_Veloarray, PBE_INDEX, Max_It, Max_N_Unknowns=0;
  int N_ScalarEqns, N_PBEqns, N_IndepntScalarEqns, VeloFunction;
  int Out_Level, Disctypes[4], start_pbe, end_pbe, N_Q3Data=1, dGDisc=0;
  int Eq_pos, begin, end, N_Cells_Intl, N_IntlPts_All, N_IntlPts=0, N_IntlLevels=0, MaxN_PtsForNodal;

  
  char *PsBaseName, *GrapeBaseName, *GnuBaseName, *ReadGrapeBaseName, RankBaseName[1000];
  char *VtkBaseName, *MatlabBaseName, *GmvBaseName;
  char *PRM, *GEO;
  
  char WString[] = "w";
  char UString[] = "u";
  char PString[] = "p";
  char Readin[] = "readin.dat";
  char MassMatrix[] = "Mass matrix";
  char Mass[] = "Mass";
  char Name[] = "name";
  char NameStrings[5][10]= {"T", "C1", "PSD", "C3", "C4"};
  char VString[] = "V";
  char GString[] = "IntlGrid";
  char IntLevel[] = "PBEL";
  char PhySol[] = "";
  char SubID[] = "PBEL";

  std::ostringstream os;
  os << " ";

  bool UpdateConvection=TRUE, UpdateRhs=FALSE, Initialize_ScalarSolver=FALSE;
  bool ConvectionFirstTime=TRUE;

  bool UpdatePBEConvection=FALSE, UpdatePBERhs=FALSE, Initialize_PBEScalarSolver=TRUE;
  bool PBEConvectionFirstTime=TRUE;
  bool *DirichletBDPt; 
  
  
  //======================================================================
  // read parameter file
  //======================================================================
//   total_time = GetTime();
  if(argc>=2)
    { ret=Domain->ReadParam(argv[1]); }
    else
      { ret=Domain->ReadParam(Readin); }

      OpenFiles();
  OutFile.setf(std::ios::scientific);  
  
  ExampleFile();

  GetExampleFileData(BoundaryConditions, BoundValues, InitiaValues, Coefficients,
                     N_PBEqns, N_IndepntScalarEqns, Disctypes);
  N_ScalarEqns=N_PBEqns+N_IndepntScalarEqns;
  
  
  
  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;

  // assign names for output files
  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  GrapeBaseName = TDatabase::ParamDB->GRAPEBASENAME;
  GnuBaseName = TDatabase::ParamDB->GNUBASENAME;
  ReadGrapeBaseName = TDatabase::ParamDB->READGRAPEBASENAME;
  GmvBaseName = TDatabase::ParamDB->GMVBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  MatlabBaseName = TDatabase::ParamDB->MATLABBASENAME;

  // pointers to the routines which compute matrix-vector
  // products and the defect
  MatVect = MatVect_Scalar;
  Defect = Defect_Scalar;

  //======================================================================
  // initialize discrete forms
  //======================================================================
  for(i=0; i<N_ScalarEqns; i++)
  {
    
    InitializeDiscreteFormsScalar(DiscreteFormMatrixMRhs[i], DiscreteFormMatrixARhs[i], DiscreteFormMatrixMRhs_SUPG[i],
                                  DiscreteFormMatrixARhs_SUPG[i], Coefficients[i]);
  }
  
  
  //======================================================================
  // read boundary parameterization and initialize coarse grid
  //======================================================================
  Domain->Init(PRM, GEO);
  //   Domain->ConvertQuadToTri(TDatabase::ParamDB->CONVERT_QUAD_TO_TRI);

#ifdef __SIMPATURS__
//   for(i=0;i<TDatabase::ParamDB->P9;i++)
//     Domain->RefineallxDirection();
#endif
 
  
  // refine grid up to the coarsest level
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();

 
  if(TDatabase::ParamDB->WRITE_PS)
  {
    // write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
  }
    
  
  // initializ time
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = TDatabase::TimeDB->TIMESTEPLENGTH;
  SetTimeDiscParameters();
  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR;

//   t3 = GetTime();
//   total_time = t3 - total_time;
  SetPolynomialDegree();

  // check the example file, to activate
  VeloFunction=TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD; 
  
  
   // set type of multilevel
  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SCALAR;
  if(TDatabase::ParamDB->SOLVER_TYPE==AMG||
      TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
   mg_type = 0;  
  
  
  if(mg_type)
    mg_level = 1;
  else
    mg_level = 0;
    
  LEVELS = TDatabase::ParamDB->LEVELS;  
  
  // initialize multigrid
  if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
    Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
    i=1;
    MG = new TMultiGrid2D(i, N_Paramters, Parameters);
    
    switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
      {
          case 11:
           zerostart = 1;
          break;
          case 16:
           zerostart = 0;
          break;
      }    
    }  
  
  FirstSolve = TDatabase::ParamDB->SC_FIRST_SOLUTION_LEVEL_SCALAR;   
  
  
  
  Scalar_Spaces = new TFESpace2D**[LEVELS+1];
  sqstructureA = new TSquareStructure2D**[LEVELS+1];

  N_Uarray = new int*[LEVELS+1];
  N_Active = new int*[LEVELS+1];
  MatricesA = new TSquareMatrix2D**[LEVELS+1];
  MatricesM = new TSquareMatrix2D**[LEVELS+1];
  MatricesK = new TSquareMatrix2D**[LEVELS+1];
  SolArray = new double**[LEVELS+1];
  RhsArray = new double**[LEVELS+1];
//   IntlX = new double*[LEVELS+1];
//   IntlY = new double*[LEVELS+1];
  
  B = new double*[LEVELS+1];
  B_Pbe = new double*[LEVELS+1];
  OldRhsArray_Pbe = new double*[LEVELS+1];
  RhsArray_Pbe = new double*[LEVELS+1];
  
  ScalarFunctions = new TFEFunction2D**[LEVELS+1];
//   ADI_System = new TADISystem1D**[LEVELS+1];
  
  for(i=0;i<LEVELS+1;i++)
   {
    Scalar_Spaces[i] = new TFESpace2D*[N_ScalarEqns];
    sqstructureA[i]  = new TSquareStructure2D*[N_ScalarEqns];

    N_Uarray[i]  = new int[N_ScalarEqns];
    N_Active[i]  = new int[N_ScalarEqns];
    MatricesA[i]  = new TSquareMatrix2D*[N_ScalarEqns];
    MatricesM[i]  = new TSquareMatrix2D*[N_ScalarEqns];
    MatricesK[i]  = new TSquareMatrix2D*[N_ScalarEqns];
    SolArray[i]  = new double*[N_ScalarEqns];
    RhsArray[i]  = new double*[N_ScalarEqns];

 
    ScalarFunctions[i]  = new TFEFunction2D*[N_ScalarEqns];    
   } //  for(i=0;i<N_IndepntScalarEqns;i++)
  
  
  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  OutPut("N_Cells (space) : " << N_Cells <<endl);
 
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  VELOCITYORDER = TDatabase::ParamDB->VELOCITY_SPACE;

#ifdef __SIMPATURS__
//   BoundCondition_LminLMax(cond_Lmin, cond_Lmax);
#endif
  
#ifdef __PBS__
    PBE_INDEX = N_IndepntScalarEqns;
#endif
  //=========================================================================
  // construct all finite element spaces
  // loop over all levels (not a multigrid level but for convergence study)
  //======================================================================
  mg_level = LEVELS + mg_level;
  
  for(l=0;l<mg_level;l++)
  {
    if (l<LEVELS)
    {
      OutPut("*******************************************************" << endl);
      OutPut("******           GEOMETRY  LEVEL ");
      OutPut(l << "              *******" << endl);
    }
    else
    {
      OutPut("*******************************************************" << endl);
      OutPut("******           GEOMETRY  LEVEL ");
      OutPut(l-1 << "              *******" << endl);
    }
    OutPut("******           MULTIGRID LEVEL ");
    OutPut(l << "              *******" << endl);
    OutPut("*******************************************************" << endl);
    OutPut("memory before: " << setw(10) << GetMemory() << endl);

    // refine domain regularly
    if(l && (l<LEVELS)) Domain->RegRefineAll();

  if(TDatabase::ParamDB->WRITE_PS)
  {
    // write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain" << l <<".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
  }    
    
    
       coll=Domain->GetCollection(It_Finest, 0);
//     OutPut( "number of cells: " << coll->GetN_Cells() << endl);
//     coll->GetHminHmax(&hmin,&hmax);
//     OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
//     cout << endl << endl; 
  
      for(i=0;i<N_IndepntScalarEqns;i++)
       {
 
       // fespaces for scalar equations
      //       Scalar_Spaces[i] =  new TFESpace2D(coll, Name, NameStrings[i],  BoundaryConditions[i], ORDER, NULL); 	 
        OutPut("N_Cells (space) : " << N_Cells <<endl);

        exit(0);   
    
    
       }   // for(i=0;i<N_IndepntScalarEqns;i++)
   
   
#ifdef __PBS__
  //=========================================================================
  // MultDimEqn setting begin
  // construct FESpace and matrices for  equation
  //=========================================================================  
  
    //=========================================================================
    // fespaces for population balance equation in physical space
    // FeSpace and memory for all matrices
    // assume that the convection and the reaction coefficient terms are independent
    // of internal coordinates, so lhs matrices are same for all lelvels of internal coordinate
    //=========================================================================    
    
    if ((mg_type==1)&&(i<mg_level-1))
    {
     // nonconforming space       
     Scalar_Spaces[l][PBE_INDEX] = new TFESpace2D(coll, Name, NameStrings[PBE_INDEX], BoundaryConditions[PBE_INDEX], -1, NULL);  
    }  
    else
    {
      // standard multigrid or finest level
      // get fe space of high order disc on finest geo grid      
      ORDER  = TDatabase::ParamDB->ANSATZ_ORDER;
      Scalar_Spaces[l][PBE_INDEX] = new TFESpace2D(coll, Name, NameStrings[PBE_INDEX], BoundaryConditions[PBE_INDEX], ORDER, NULL);
    }     
  
     N_Uarray[l][PBE_INDEX] = Scalar_Spaces[l][PBE_INDEX]->GetN_DegreesOfFreedom();
     N_Active[l][PBE_INDEX] = Scalar_Spaces[l][PBE_INDEX]->GetActiveBound();  
     
     if(Max_N_Unknowns<N_Uarray[l][PBE_INDEX])  Max_N_Unknowns=N_Uarray[l][PBE_INDEX];
  
     OutPut("DOF of PBE physical space : "<< setw(10) << N_Uarray[l][PBE_INDEX] << endl);
  
     
    //=========================================================================
    // memory allocate all vectors and construction for PBS space
    //=========================================================================
    // build matrices
    // first build matrix structure
    sqstructureA[l][PBE_INDEX] = new TSquareStructure2D(Scalar_Spaces[l][PBE_INDEX]);
    sqstructureA[l][PBE_INDEX]->Sort();  
  
    // two matrices used, M is the mass matrix, also a system matrix
    // the iterative solver uses M
    sqmatrixM = new TSquareMatrix2D(sqstructureA[l][PBE_INDEX]);
    MatricesM[l][PBE_INDEX] = sqmatrixM;

    // A contains the non time dependent part of the discretization
    sqmatrixA = new TSquareMatrix2D(sqstructureA[l][PBE_INDEX]);
    MatricesA[l][PBE_INDEX] = sqmatrixA;
  
    //=========================================================================
    // Construct FeSpace and all data for the internal domain
    // FESpace is same in internal coordinate for all  QuadPts/Dof of
    // Physical FESapce and for all levels
    //=========================================================================    
    if(l==0)
     {
     //=========================================================================
     // Construct the internal domain for this mg level
     //=========================================================================    
      N = int(TDatabase::ParamDB->REACTOR_P11);
      L0 = TDatabase::ParamDB->REACTOR_P12;
      L1 = TDatabase::ParamDB->REACTOR_P13;

      Generate1DMesh(Domain_Intl, L0, L1, N);
      Coll_Intl = Domain_Intl->GetCollection(It_Finest, 0);
      N_Cells_Intl= Coll_Intl->GetN_Cells();           
       
       
      FeSpace_Intl = new TFESpace1D(Coll_Intl, VString, VString, TDatabase::ParamDB->TEST_ORDER);
      N_LDof = FeSpace_Intl->GetN_DegreesOfFreedom();

      if(TDatabase::ParamDB->TEST_ORDER<-9) // DG method
       {
        FeSpace_Intl->SetAsDGSpace();
        dGDisc = 1;
        TDatabase::ParamDB->P10=0;  // no supg method
       }

      N_Intl_Levels =  FeSpace_Intl->GetN_DegreesOfFreedom();
      GetInternalNodalPts(FeSpace_Intl, N_Intl_Levels, IntlPosL);
      OutPut("Dof PBE Internal space : "<< setw(14) << N_Intl_Levels << endl);
      
      SqStruct1D_Intl = new TSquareStructure1D(FeSpace_Intl);
      SqStruct1D_Intl->Sort();
      M_Intl = new TSquareMatrix1D(SqStruct1D_Intl);
      A_Intl = new TSquareMatrix1D(SqStruct1D_Intl);

      if(TDatabase::ParamDB->P10)
       {
        S_Intl = new TSquareMatrix1D(SqStruct1D_Intl);
        K_Intl = new TSquareMatrix1D(SqStruct1D_Intl);
       }
       
       
       TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
       TDatabase::IteratorDB[It_LE]->SetParam(Domain);
       TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
       TDatabase::IteratorDB[It_Between]->SetParam(Domain);
       TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);              

     } //  // if(l==0)
       
  
  
   if(l==TDatabase::ParamDB->REACTOR_P20)
     { 
      OutPut("*******************************************************" << endl);
      OutPut("*****  Internal solution at MULTIGRID LEVEL ");
      OutPut(l << "     *****" << endl);
      OutPut("*******************************************************" << endl);
      OutPut("memory before: " << setw(10) << GetMemory() << endl);           
       
      if(TDatabase::ParamDB->REACTOR_P5==0)
      { 
       GetNodalPtsADI_System(N_IntlPts, FeSpace_Intl, Scalar_Spaces[l][PBE_INDEX], MaxN_PtsForNodal, IntlX, IntlY); 
      }
     else
      {
       GetQuadPtsADI_System(N_IntlPts, Scalar_Spaces[l][PBE_INDEX], IntlX, IntlY); 
      }     
     
        //OutPut("Dof PBE Internal space N_IntlPts : "<<  N_IntlPts << endl);
      ADI_System = new TADISystem1D*[N_IntlPts];

//        for(i=0; i<N_IntlPts; i++)
//         {
// 	 if(fabs(IntlX[i]-200)<1.e-12  && fabs(IntlY[i]-0.5)<1.e-12)
// 	  {
//            OutPut(i << " IntlX[i] " <<IntlX[i]<< " IntlY[i] " <<IntlY[i] <<endl);
// 	   f_Integ_NodalPt = i;
// 	   break;
// 	  }
// 	 
// 	}
     
     
     Sol_IntlLoc = new double[N_Intl_Levels];
     OldSol_IntlLoc = new double[N_Intl_Levels];
     B_IntlLoc = new double[N_Intl_Levels];
     defect_IntlLoc = new double[N_Intl_Levels];

     memset(Sol_IntlLoc, 0, N_Intl_Levels*SizeOfDouble);
     memset(OldSol_IntlLoc, 0, N_Intl_Levels*SizeOfDouble);
     memset(B_IntlLoc, 0, N_Intl_Levels*SizeOfDouble);
     memset(defect_IntlLoc, 0, N_Intl_Levels*SizeOfDouble);

     FeFunction_Intl = new TFEFunction1D(FeSpace_Intl, NameStrings[PBE_INDEX], NameStrings[PBE_INDEX], Sol_IntlLoc,  N_Intl_Levels);

     for(i=0; i<N_IntlPts; i++)
      ADI_System[i] = new TADISystem1D(FeFunction_Intl, M_Intl, A_Intl, S_Intl, K_Intl, IntlX[i], IntlY[i],
                              InitialCondition_Psd_Intl,  BoundValue_LMin, BoundValue_LMax,
                              Sol_IntlLoc, OldSol_IntlLoc, B_IntlLoc, defect_IntlLoc);

      C_IntlPtValues = new double[N_IntlPts];
      T_IntlPtValues = new double[N_IntlPts];
      C_Sat_IntlPtValues = new double[N_IntlPts];
      velo_IntlPts = new double[2*N_IntlPts];
      grad_velo_IntlPts = new double[4*N_IntlPts];
      i = N_IntlPts*N_Intl_Levels;

      PBE_IntlPtValuesT = new double[i];
      OldPBE_IntlPtValuesT = new double[i];
      PBE_IntlPtValues = new double[i];
      PBE_IntlPtRhsValues = new double[i];
      DirichletBDPt = new bool[N_IntlPts];

      i= N_LDof*N_IntlPts;
      PBE_LSolT =  new double[i];
      PBE_IntlPtRhsValuesT = new double[i];
      OldPBE_IntlPtRhsValuesT = new double[i];
     
      Out_Level = N_Intl_Levels/4;

     if(cond_Lmin==DIRICHLET && !dGDisc)
      { start_pbe = 1; }
     else
      { start_pbe = 0; }

     if(cond_Lmax==DIRICHLET && !dGDisc) { end_pbe = N_Intl_Levels-1; }
     else { end_pbe = N_Intl_Levels; }
//        exit(0);    
     } //  if(l==TDatabase::ParamDB->REACTOR_P20)
    
   //========================================================================================
   // memory allocate for all vectors and construction of PBE fefunction in
   // physical space (for all internal levels)
   //========================================================================================
   N_U = N_Uarray[l][PBE_INDEX];
   N_U_PBE = N_Intl_Levels*N_U;

   SolArray[l][PBE_INDEX] = new double[N_U_PBE];
   SolPbe = SolArray[l][PBE_INDEX];
   RhsArray_Pbe[l] = new double[N_U_PBE];
   B_Pbe[l] = new double[N_U_PBE];
  
   OldRhsArray_Pbe[l] = new double[N_U_PBE];
   memset(RhsArray_Pbe[l], 0, N_U_PBE*SizeOfDouble);
   memset(OldRhsArray_Pbe[l], 0, N_U_PBE*SizeOfDouble); 
   memset(SolPbe, 0, N_U_PBE*SizeOfDouble);    
   
   ScalarPbeFunctions = new TFEFunction2D* [N_Intl_Levels];
   for(i=0; i<N_Intl_Levels; i++)
    ScalarPbeFunctions[i] = new  TFEFunction2D(Scalar_Spaces[l][PBE_INDEX], NameStrings[PBE_INDEX], NameStrings[PBE_INDEX], SolPbe+i*N_U, N_U);
   
   
  // Identify the dirichlet internal points -begin
    SolPbe_Lmin = new double[N_U];
    SolPbe_Lmax = new double[N_U];

    if(TDatabase::ParamDB->REACTOR_P5==0)
     {
      // set 1 for all DOF at L_Min
      for(i=0;i<N_U;i++)
        SolPbe[i] = 1.;

      // set 1 for Dirichlet DOF and else 0 at L_Max
      SolPbe_OneLevel = SolPbe+(N_Intl_Levels-1)*N_U;
      memset(SolPbe_OneLevel, 0, N_U*SizeOfDouble);
      
      for(i=N_Active[l][PBE_INDEX];i<N_U;i++)
        SolPbe_OneLevel[i] = 1.;

      GetPtsValuesForIntl(N_Intl_Levels, ScalarPbeFunctions, SolPbe, N_U, PBE_IntlPtValuesT, N_IntlPts);

//       for(i=0; i<N_IntlPts; i++)
//       {
//         if( fabs(PBE_IntlPtValuesT[i*N_Intl_Levels] - PBE_IntlPtValuesT[i*N_Intl_Levels + N_Intl_Levels-1])<1.e-8 )
//         {
//           //cout <<i<< " IntlX[i] " <<IntlX[i]<< " IntlY[i] " <<IntlY[i] <<endl;
//           // no computation need for these points in internal direction,
//           //i.e., for nodal points on inlet, outlet, walls
//           DirichletBDPt[i] = TRUE;
//         }
// 
//       }
     }
    else //if(TDatabase::ParamDB->REACTOR_P5==0)
     {
      // assumed that no quad points of X-direction lie on the boundary
      for(i=0; i<N_IntlPts; i++)
      DirichletBDPt[i] = FALSE;
     }   
   
 
   
   if(l==mg_level-1)
    {
     SolPbe_Output = new double[N_U];
     memset(SolPbe_Output, 0, N_U*SizeOfDouble);
     PbeFunctions = new  TFEFunction2D(Scalar_Spaces[l][PBE_INDEX], NameStrings[PBE_INDEX], NameStrings[PBE_INDEX], SolPbe_Output, N_U);

     // prepare output, only the concentration will be saved
     Output_PBS = new TOutput2D(0,  1, 1, 1, Domain);
     Output_PBS->AddFEFunction(PbeFunctions);   
   
    } // if(l==mg_level-1)
#endif  // __PBS__    
    
  }// for(l=0;l<LEVELS;l++)
  
  

  
  CloseFiles();
  return 0;
}






