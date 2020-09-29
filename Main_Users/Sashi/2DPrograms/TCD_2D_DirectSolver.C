// =======================================================================
//
// Purpose:     main program with parallel solver (no multigrid solver)
//              more than once scalar equations and with a given velocity
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 03.09.2009
//        :     MPI implementation started on 07.09.2009
//        :     operator splitting first and second order (08.03.2010)
//        :     MPI disabled on 12.08.2013, added OMP
// =======================================================================
// #ifdef _OMP
// #include <omp.h>
// #endif

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
//#include <ConvDiff2D_Routines.h>
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

#include <BrAgg.h>
#include <exported.h>

#include <TimeDiscRout.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdio.h>
#include <stdlib.h>
// #include <malloc.h>

double bound = 0;
 
#include <MainUtilities.h>
//#include <TimeUtilities.h>
#include <FEM_TVD_FCT.h>

#ifdef _MPI
#include <MeshPartition.h>
#include <ParFECommunicator2D.h>
#include <Scalar_ParSolver.h>
#include <ParVector2D.h>
#endif

#ifdef _OMP
#include <ParDirectSolver.h>
#endif

#define AMG 0
#define GMG 1
#define DIRECT 2

#include <sys/stat.h>
#include <sys/types.h>

#ifdef __MAC64__
// #include <Accelerate/Accelerate.h>
#endif

// =======================================================================
// include current example
// =======================================================================

// #include "../Examples/TCD_2D/ansatz1.h"
// #include "../Examples/TCD_2D/SinCos1.h"
//#include "../Examples/TCD_2D/SinCos2.h"
// #include "../Examples/TCD_2D/SinCos4.h"
// #include "../Examples/TCD_2D/Robin_indian.01.h"
//#include "../Examples/TCD_2D/Robin_indian.00.h"
//#include "../Examples/TCD_2D/Robin_indian.02.h"
//#include "../Examples/TCD_2D/Robin_indian.03.h"
//#include "../Examples/TCD_2D/Robin_indian.04.h"
//#include "../Examples/TCD_2D/Robin_scaled.00.h"
// #include "../Examples/TCD_2D/RotatingBodies.h"
//#include "../Examples/TCD_2D/Bulk_Academic_Test.h"
//#include "../Examples/TCD_2D/JohnMaubachTobiska1997inst.h"
//#include "../Examples/TCD_2D/JohnMaubachTobiska1997inst2.h"
//#include "../Examples/ConvDiff2D/Hemker1996.h"
// #include "../Examples/TCD_2D/Sin3.h"
//#include "../Examples/TCD_2D/Sin3_0.h"

// #include "../Examples/TCD_2D/Time1_Par.h"
// #include "../Examples/TCD_2D/Time2.h"
// #include "../Examples/TCD_2D/Gaussian-BenchMark.h"

// #include "../Examples/TCD_2D/SimPaTurS.h"
// #include "../Examples/TCD_2D/2D_1D_ADI_SimPaTurS.h"
// #include "../Examples_All/TCD_2D/2D_1D_ADI_SimPaTurS_Axial3D.h"
// #include "../Examples/TCD_2D/2D_1D_ConstT_UnitSqr.h"
 #include "../TCD_2D/2D_1D_Smmoth.h"
// #include "../Examples/TCD_2D/2D_1D_ConstT_UnitSqr_IntlOnly.h"
// #include "../Examples/TCD_2D/terahertz.h"
// #include "../Examples/TCD_2D/TeraHertz.h"
// #include  "../Examples/TCD_2D/Levelset.h"
// #include  "../Examples/TCD_2D/TimeDomian.h"
// ======================================================================
// utilities for main program
// ======================================================================

void BoundCondition_NSE(int i, double t, BoundCond &cond)
{

  switch(i)
  {
    case 0:
    case 2:
    case 3:
    case 4:
    case 5:
      cond = DIRICHLET;
      break;
    case 1:
      cond = NEUMANN;
      TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=0;
      break;
    default:
      cout << "wrong boundary part number" << endl;
      exit(0);
      break;
  }
}

// update the LHS in the Robin BC to the stiffness matrix     
#ifdef __ROBINBC__    

void RobinBCInt(TSquareMatrix2D *A, BoundCondFunct2D *BDCond)
{
  int i, j, k, l, m, N_Edges, IJoint, N_NJoints, comp;
  int *BeginIndex, *GlobalNumbers, *DOF, *RowPtr, *KCol;
  int N_Active, N_Cells, N_BaseFunct, *N_BaseFuncts;
  int JointNumbers[MAXN_JOINTS];
  int N_LinePoints, AnsatzDOF, TestDOF, index1, index2;
   
  double *LineWeights, *zeta, r_axial;   
  double t0, t1, X_B[MaxN_QuadPoints_1D], Y_B[MaxN_QuadPoints_1D];
  double **uref, x0, y0, x1, y1, hE, Mult;
  double uorig[MaxN_BaseFunctions2D], *MatValues, area;
  
  BaseFunct2D *BaseFuncts; 
  TFESpace2D *fespace;
  TCollection *Coll;
  TFEDesc2D *FeDesc;  
  TBaseCell *cell;  
  TJoint *joint;
  TBoundEdge *boundedge;
  TBoundComp *BoundComp;
  TIsoBoundEdge *isoboundedge;
  BoundCond Cond0, Cond1;
  FE2D FEId;
  TFE2D *ele;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  BF2DRefElements RefElement;  
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  
// heat stress number (or) BIOT NUMBER
  double C0 = TDatabase::ParamDB->P6;
  
  
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();  
  
  
  fespace = A->GetFESpace();
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  N_Active = fespace->GetActiveBound();
  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();  

  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();  
  MatValues = A->GetEntries();

  for(i=0;i<N_Cells;i++)
   {
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    IJoint = 0;
    for(j=0;j<N_Edges;j++)
     {
      joint = cell->GetJoint(j);
      
      if(joint->GetType() == BoundaryEdge||
         joint->GetType() == InterfaceJoint ||
         joint->GetType() == IsoBoundEdge)
       {Example
         if(joint->GetType() == BoundaryEdge||
            joint->GetType() == InterfaceJoint)
          {
            boundedge = (TBoundEdge *)joint;
            BoundComp = boundedge->GetBoundComp();
            boundedge->GetParameters(t0, t1);
          }
          else
          {
            isoboundedge = (TIsoBoundEdge *)joint;
            BoundComp = isoboundedge->GetBoundComp();
            isoboundedge->GetParameters(t0, t1);
          }      
      
       }// if(joint->GetType() == BoundaryEdge||
       
      // get id of the boundary component
      comp = BoundComp->GetID();
      // get type of the boundary condition at the beginning
      // and at the end of the current edge
      BDCond(comp, t0, Cond0);       
       
        if(Cond0 == ROBIN)
        { 
         // cout << "comp " << comp  <<endl;
          JointNumbers[IJoint] = j;
          IJoint++;
        }    
      } // endfor j
      
      
     N_NJoints = IJoint;
     if(N_NJoints > 0)
      {
       FEId = fespace->GetFE2D(i, cell);
       N_BaseFunct = N_BaseFuncts[FEId];
       ele = TFEDatabase2D::GetFE2D(FEId);
       RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);   
       
       DOF = GlobalNumbers + BeginIndex[i];
 
       l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
       LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
       qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
       qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
       TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);

        switch(RefElement)
        {
         case BFUnitTriangle:

            RefTrans = TriaAffin;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaAffin *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:

            RefTrans = QuadAffin;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadAffin *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowed" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
          } // endswitch
          
          
       for(j=0;j<N_NJoints;j++)
        {     
         IJoint = JointNumbers[j]; 

         switch(RefElement)
         {
         case BFUnitTriangle:
            ((TTriaAffin *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
          break;

          case BFUnitSquare:
            ((TQuadAffin *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
          break;
         } // endswitch


         uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId], LineQuadFormula, IJoint);


         // multiply value with weights from quadrature formula
         // and determinant from integral transformation to the
         // unit edge (-1,1)       
         cell->GetVertex(IJoint)->GetCoords(x0, y0);
         cell->GetVertex((IJoint+1) % N_Edges)->GetCoords(x1, y1);
         // compute (half of the) length of the boundary edge
         hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;

        for(k=0;k<N_LinePoints;k++)
         {
           // D00
           for(l=0;l<N_BaseFunct;l++)
            uorig[l] = uref[k][l];

          if(TDatabase::ParamDB->Axial3D)
           { Mult = C0*fabs(X_B[k])*hE*LineWeights[k]; }
          else
           { Mult = C0*r_axial*hE*LineWeights[k]; }
 

          for(l=0;l<N_BaseFunct;l++)
           {
            TestDOF = DOF[l];
            index2 = RowPtr[TestDOF+1];
             for(m=0;m<N_BaseFunct;m++)
              {
               AnsatzDOF = DOF[m];
               index1 = RowPtr[TestDOF];
               if(index1+1 == index2) continue;
               while(KCol[index1] != AnsatzDOF) index1++;

               MatValues[index1] += (Mult*uorig[m]*uorig[l]);
               //cout << "A: (" << TestDOF << " " << AnsatzDOF << ") =  " << Mult*uorig[m]*uorig[l] << endl;
              } // endfor m
          } // endfor l 
         } //  for(k=0;k<N_LinePoints;k
        } // for(j=0;j<N_NJoints
      }//  if(IJoint > 0)  
  
  } // endfor i  
  
//   cout << " RobinBCInt " << endl;
//   exit(0);
} // RobinBCInt

#endif




#ifdef __SIMPATURS__
void Trans_XLtoLX(double *PBE_IntlPtValuesT, double *PBE_IntlPtValues, int N_PBE_SpatialPts, int N_PBE_IntlLevels)
{
  int i, j;
  double *sol;

  for(i=0; i<N_PBE_SpatialPts; i++)
    for(j=0; j<N_PBE_IntlLevels; j++)
      PBE_IntlPtValues[j*N_PBE_SpatialPts + i] = PBE_IntlPtValuesT[i*N_PBE_IntlLevels + j];
}


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


void GetQuadPtsADI_System(int &N_PBE_SpatialPts, TFESpace2D *FESpace2D, double *&IntlX, double *&IntlY)
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
  N_Cells = Coll->GetN_Cells();

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
      N_PBE_SpatialPts = N_Points*N_Cells;
      //     cout<< "N_PBE_SpatialPts " << N_PBE_SpatialPts << endl;
      IntlX = new double[N_PBE_SpatialPts];
      IntlY = new double[N_PBE_SpatialPts];
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

  //  cout<< "N_PBE_SpatialPts " << m <<endl;
  //   for(i=0; i<N_PBE_SpatialPts; i++)
  //    cout<< i << " IntlX " << IntlX[i]<< " IntlY " << IntlY[i]  <<endl;
  //  exit(0);

}


void GetPtsValuesForIntl_Quad(int N_Levels, TFEFunction2D **ScalarFunctions, double *ValuesT, int N_PBE_SpatialPts)
{
  int i, j, k, l, m, N_Cells, N_Points, N_LocalDOFs, N_Sets=1;
  int *BeginIndex, *GlobalNumbers, *DOF;
  int PolynomialDegree;

  double *xi, *eta, *weights, maxval=0.;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], AbsDetjk[MaxN_QuadPoints_2D];
  double **origvalues, *sol, *org, val;

  bool Needs2ndDer[1];

  TFESpace2D *PBE_Spaces;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId_PBE;
  TFE2D *Element_PBE;
  BF2DRefElements RefElement;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;
  TRefTrans2D *F_K;
  BaseFunct2D BaseFunct;

  // assume all  ScalarFunctions use same fespace2D
  PBE_Spaces = ScalarFunctions[0]->GetFESpace2D();
  BeginIndex = PBE_Spaces->GetBeginIndex();
  GlobalNumbers = PBE_Spaces->GetGlobalNumbers();

  Needs2ndDer[0] = FALSE;

  // assume that both fespace2D and PBE_Spaces use same coll
  Coll = PBE_Spaces->GetCollection();
  N_Cells = Coll->GetN_Cells();
  m = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId_PBE = PBE_Spaces->GetFE2D(i, cell);
    Element_PBE = TFEDatabase2D::GetFE2D(FEId_PBE);
    N_LocalDOFs = Element_PBE->GetN_DOF();
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId_PBE);
    PolynomialDegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_PBE);
    BaseFunct = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D(FEId_PBE);

    switch(RefElement)
    {
      case BFUnitSquare:
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*PolynomialDegree);
        F_K = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
        ((TQuadBilinear *)F_K)->SetCell(cell);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_Points, weights, xi, eta);
        ((TQuadBilinear *)F_K)->GetOrigValues(BaseFunct, N_Points, xi, eta, N_LocalDOFs, QuadFormula);
        break;

      case BFUnitTriangle:
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(3*PolynomialDegree-1);
        F_K = TFEDatabase2D::GetRefTrans2D(TriaAffin);
        ((TTriaAffin *)F_K)->SetCell(cell);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_Points, weights, xi, eta);
        ((TTriaAffin *)F_K)->GetOrigValues(BaseFunct, N_Points, xi, eta, N_LocalDOFs, QuadFormula);
        break;
    }                                             // endswitch

    //cout << "QuadFormula: " << QuadFormula << endl;

    origvalues=TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);
    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0; j<N_Levels; j++)
    {
      sol = ScalarFunctions[j]->GetValues();

      for(k=0; k<N_Points; k++)
      {
        org = origvalues[k];
        val = 0.;
        for(l=0;l<N_LocalDOFs;l++)
        {
          val += org[l]*sol[ DOF[l] ];            // for M
        }

        ValuesT[(m+k)*N_Levels + j] = val;

        if(maxval<val) maxval=val;

      }                                           //  for(k=0; k<N_Points; k++)
    }                                             // for(j=0; j<N_Levels; j++)

    m += N_Points;
  }                                               //  for(i=0;i<N_Cells;i++)

  // cout <<  " GetPtsValuesForIntl_Quad " << maxval <<endl;

  //    cout <<  m << " N_PBE_SpatialPts " << N_PBE_SpatialPts <<endl;
  // exit(0);
}


void GetQuadPtsValuesForIntl(TFEFunction2D *Heat, TFEFunction2D *Concentration, TFESpace2D *PBE_Spaces,
double *Values, double *Sat_Values, int N_PBE_SpatialPts)
{
  int i, j, l, m, N_Cells, N_Points, N_LocalDOFs, T_N_LocalDOFs, N_Pts;
  int *BeginIndex, *GlobalNumbers, *DOF;
  int *T_BeginIndex, *T_GlobalNumbers, *T_DOF;
  int PolynomialDegree, ApproxOrder;

  double *xi, *eta, *weights, *C_NodalValues, s, *T_NodalValues;
  double BasisValues[MaxN_BaseFunctions2D], T_val;
  double T_ref = TDatabase::ParamDB->REACTOR_P23;
  double C_ref = TDatabase::ParamDB->REACTOR_P25;

  TFESpace2D *fespace2D, *T_fespace2D;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId, T_FEId, FEId_PBE;
  TFE2D *Element, *T_Element, *Element_PBE;
  TBaseFunct2D *bf, *T_bf;
  TNodalFunctional2D *nf_PBE;
  BF2DRefElements RefElement;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;

  T_NodalValues =  Heat->GetValues();
  T_fespace2D = Heat->GetFESpace2D();
  T_BeginIndex = T_fespace2D->GetBeginIndex();
  T_GlobalNumbers = T_fespace2D->GetGlobalNumbers();

  C_NodalValues =  Concentration->GetValues();
  fespace2D = Concentration->GetFESpace2D();
  BeginIndex = fespace2D->GetBeginIndex();
  GlobalNumbers = fespace2D->GetGlobalNumbers();

  memset(Values, 0, N_PBE_SpatialPts*SizeOfDouble);
  memset(Sat_Values, 0, N_PBE_SpatialPts*SizeOfDouble);

  // assume that all fespace2Ds and PBE_Spaces use same coll
  Coll = fespace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  N_Pts=0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId_PBE = PBE_Spaces->GetFE2D(i, cell);
    Element_PBE = TFEDatabase2D::GetFE2D(FEId_PBE);
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId_PBE);
    PolynomialDegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_PBE);

    switch(RefElement)
    {
      case BFUnitSquare:
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*PolynomialDegree);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_Points, weights, xi, eta);
        break;

      case BFUnitTriangle:
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(3*PolynomialDegree-1);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_Points, weights, xi, eta);
        break;
    }                                             // endswitch

    FEId = fespace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    bf = Element->GetBaseFunct2D();
    N_LocalDOFs = Element->GetN_DOF();
    DOF = GlobalNumbers + BeginIndex[i];

    T_FEId = T_fespace2D->GetFE2D(i, cell);
    T_Element = TFEDatabase2D::GetFE2D(T_FEId);
    T_bf = T_Element->GetBaseFunct2D();
    T_N_LocalDOFs = T_Element->GetN_DOF();
    T_DOF = T_GlobalNumbers + T_BeginIndex[i];

    for(j=0;j<N_Points;j++)
    {
      // C value at this point
      bf->GetDerivatives(D00, xi[j], eta[j], BasisValues);
      for(l=0;l<N_LocalDOFs;l++)
      {
        m = DOF[l];
        s = C_NodalValues[m];
        Values[N_Pts] += BasisValues[l]*s;
      }

      // first find T value at this point, then compute C_sat
      T_bf->GetDerivatives(D00, xi[j], eta[j], BasisValues);
      for(l=0;l<T_N_LocalDOFs;l++)
      {
        m = T_DOF[l];
        s = T_NodalValues[m];
        Sat_Values[N_Pts] += BasisValues[l]*s;
      }

      N_Pts++;
    }                                             // for(j=0;j<N_Points;j++)
  }                                               // for(i=0;i<N_Cells;i++)

  for(i=0;i<N_PBE_SpatialPts;i++)
  {
    // now compute C_Sat from T
    Sat_Values[i] = (1.3045*(T_ref*Sat_Values[i] - 273.15) + 35.3642)/C_ref;
    //     cout << i <<  "Sat_Values[i] " << Sat_Values[i]<< endl;
  }

  // cout << " GetQuadPtsValuesForIntl " <<endl;
  // exit(0);
}


void GetSolFromQuadPtVales(int N_Levels, int MaxN_PtsForNodal, TFEFunction2D **ScalarFunctions,
double *Sol, int N_U, double *ValuesT, int N_PBE_SpatialPts)
{
  int i, j, k, l, m, n, N_Cells, N_Points, N_LocalDOFs, N_Sets=1;
  int *BeginIndex, *GlobalNumbers, *DOF;
  int PolynomialDegree, *N_IncidentArray;

  double *xi, *eta, *weights, *RHS, maxval=0.;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], AbsDetjk[MaxN_QuadPoints_2D];
  double **origvalues, *sol, *org, w, val, G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];;

  bool Needs2ndDer[1];

  TFESpace2D *PBE_Spaces;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId_PBE;
  TFE2D *Element_PBE;
  BF2DRefElements RefElement;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;
  TRefTrans2D *F_K;
  BaseFunct2D BaseFunct;

  // assume all  ScalarFunctions use same fespace2D
  PBE_Spaces = ScalarFunctions[0]->GetFESpace2D();
  BeginIndex = PBE_Spaces->GetBeginIndex();
  GlobalNumbers = PBE_Spaces->GetGlobalNumbers();

  Needs2ndDer[0] = FALSE;
  RHS = new double[N_Levels*MaxN_QuadPoints_2D];
  N_IncidentArray = new int[N_U];
  memset(N_IncidentArray, 0, N_U*SizeOfInt);
  memset(Sol, 0, N_U*N_Levels*SizeOfDouble);

  // assume that both fespace2D and PBE_Spaces use same coll
  Coll = PBE_Spaces->GetCollection();
  N_Cells = Coll->GetN_Cells();
  
  m=0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId_PBE = PBE_Spaces->GetFE2D(i, cell);
    Element_PBE = TFEDatabase2D::GetFE2D(FEId_PBE);

    N_LocalDOFs = Element_PBE->GetN_DOF();
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId_PBE);
    PolynomialDegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_PBE);
    BaseFunct = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D(FEId_PBE);

    switch(RefElement)
    {
      case BFUnitSquare:
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*PolynomialDegree);
        F_K = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
        ((TQuadBilinear *)F_K)->SetCell(cell);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_Points, weights, xi, eta);
        ((TQuadBilinear *)F_K)->GetOrigValues(BaseFunct, N_Points, xi, eta, N_LocalDOFs, QuadFormula);
        ((TQuadBilinear *)F_K)->GetOrigFromRef(N_Points, xi, eta, X, Y, AbsDetjk);
        break;

      case BFUnitTriangle:
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(3*PolynomialDegree-1);
        F_K = TFEDatabase2D::GetRefTrans2D(TriaAffin);
        ((TTriaAffin *)F_K)->SetCell(cell);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_Points, weights, xi, eta);
        ((TTriaAffin *)F_K)->GetOrigValues(BaseFunct, N_Points, xi, eta, N_LocalDOFs, QuadFormula);
        ((TTriaAffin *)F_K)->GetOrigFromRef(N_Points, xi, eta, X, Y, AbsDetjk);
        break;
    }                                             // endswitch

    //cout << "QuadFormula: " << QuadFormula << endl;

    origvalues=TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);
    DOF = GlobalNumbers + BeginIndex[i];

    memset(RHS, 0, N_Levels*N_LocalDOFs*SizeOfDouble);
    memset(G, 0, N_LocalDOFs*N_LocalDOFs*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      w = AbsDetjk[j]*weights[j];
      org = origvalues[j];
      sol = ValuesT+(m*N_Levels);
      m++;

      for(k=0;k<N_LocalDOFs;k++)
      {
        val = w*org[k];

        // multiple r.h.s
        for(l=0;l<N_Levels;l++)
        {
          RHS[k*N_Levels + l] += val*sol[l];
          if(maxval<sol[l]) maxval=sol[l];
        }

        //single matrix
        for(n=0;n<N_LocalDOFs;n++)
        {
          G[k*N_LocalDOFs + n] += val*org[n];
        }                                         // for(n=0;n<N_LocalDOFs;n++)
      }                                           //for(k=0;k<N_Levels;k++)
    }                                             // for(j=0;j<N_Points;j++)

    //    for(k=0;k<N_LocalDOFs;k++)
    //     for(n=0;n<N_LocalDOFs;n++)
    //     cout<< "(" << k<< " , " << n << ") = " << G[k*N_LocalDOFs + n] << " AT " <<  G[n*N_LocalDOFs + k] <<endl;

    SolveMultipleSystemsNew(G, RHS, N_LocalDOFs, N_LocalDOFs, N_Levels, N_Levels);

    for(j=0;j<N_Levels;j++)
    {
      sol = Sol+(N_U*j) ;

      for(k=0;k<N_LocalDOFs;k++)
      {
        l = DOF[k];
        if(j==0)
        {
          N_IncidentArray[l]++;
          //           OwnDof[l] = TRUE;
        }
        sol[l] += RHS[k*N_Levels + j];
      }
    }                                             // for(j=0;j<N_Levels;j++)
  }                                               // for(i=0;i<N_Cells;i++)

  for(j=0;j<N_Levels;j++)
  {
    sol = Sol+(N_U*j);

    for(k=0;k<N_U;k++)
    {
      //         if(j==0 && OwnDof[k])
      if(j==0)
      {
        if(N_IncidentArray[k]==0)
        {
          cout<< "error in GetSolFromQuadPtVales " <<endl;
          exit(0);
        }
      }
      sol[k] /= double(N_IncidentArray[k]);

    }                                             // for(k=0;k<N_LocalDOFs;k++)
  }                                               // for(j=0;j<N_Levels;j++)

  delete [] RHS;
  delete [] N_IncidentArray;

  //  cout << " GetSolFromQuadPtVales " << maxval << endl;
}   //void GetSolFromQuadPtVales(


void GetSolFromQuadPtVales(TFEFunction2D *ScalarFunction, double *Values, TFESpace2D *PBE_Spaces)
{
  int i, j, k, l, m, n, N_Cells, N_Points, N_LocalDOFs, N_U;
  int *BeginIndex, *GlobalNumbers, *DOF, PolynomialDegree, *N_IncidentArray;

  double *xi, *eta, *weights, RHS[MaxN_BaseFunctions2D], G[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
  double val, sol, PtValues[MaxN_PointsForNodal2D], *Sol;
  double FunctionalValues[MaxN_PointsForNodal2D];
  double **origvalues, *org;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], AbsDetjk[MaxN_QuadPoints_2D], w;

  //   bool  *OwnDof;

  TFESpace2D *FE_Space;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId, FEId_PBE;
  TFE2D *Element, *Element_PBE;
  TNodalFunctional2D *nf;
  BF2DRefElements RefElement;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;
  BaseFunct2D BaseFunct;
  TRefTrans2D *F_K;

  // assume all  ScalarFunctions use same fespace2D
  Sol =  ScalarFunction->GetValues();
  FE_Space = ScalarFunction->GetFESpace2D();
  BeginIndex = FE_Space->GetBeginIndex();
  GlobalNumbers = FE_Space->GetGlobalNumbers();
  N_U = FE_Space->GetN_DegreesOfFreedom();
  N_IncidentArray = new int[N_U];
  memset(N_IncidentArray, 0, N_U*SizeOfInt);
  memset(Sol, 0, N_U*SizeOfDouble);

  // assume that both fespace2D and FE_Space use same coll
  Coll = FE_Space->GetCollection();
  N_Cells = Coll->GetN_Cells();

  m = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId_PBE = PBE_Spaces->GetFE2D(i, cell);
    Element_PBE = TFEDatabase2D::GetFE2D(FEId_PBE);
    RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId_PBE);
    PolynomialDegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId_PBE);

    switch(RefElement)
    {
      case BFUnitSquare:
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*PolynomialDegree);
        F_K = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
        ((TQuadBilinear *)F_K)->SetCell(cell);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_Points, weights, xi, eta);
        ((TQuadBilinear *)F_K)->GetOrigFromRef(N_Points, xi, eta, X, Y, AbsDetjk);
        break;

      case BFUnitTriangle:
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(3*PolynomialDegree-1);
        F_K = TFEDatabase2D::GetRefTrans2D(TriaAffin);
        ((TTriaAffin *)F_K)->SetCell(cell);
        qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_Points, weights, xi, eta);
        ((TTriaAffin *)F_K)->GetOrigFromRef(N_Points, xi, eta, X, Y, AbsDetjk);
        break;
    }                                             // endswitch

    FEId = FE_Space->GetFE2D(i, cell);
    BaseFunct = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D(FEId);

    origvalues=TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);
    DOF = GlobalNumbers + BeginIndex[i];

    memset(RHS, 0, N_LocalDOFs*SizeOfDouble);
    memset(G, 0, N_LocalDOFs*N_LocalDOFs*SizeOfDouble);

    for(j=0;j<N_Points;j++)
    {
      w = AbsDetjk[j]*weights[j];
      org = origvalues[j];
      sol = Values[m];
      m++;

      for(k=0;k<N_LocalDOFs;k++)
      {
        val = w*org[k];
        RHS[k] +=sol*val;

        //single matrix
        for(n=0;n<N_LocalDOFs;n++)
        {
          G[k*N_LocalDOFs + n] += val*org[n];
        }                                         // for(n=0;n<N_LocalDOFs;n++)
      }                                           // for(k=0;k<N_LocalDOFs;k++)
    }                                             // for(j=0;j<N_Points;j++)

    SolveMultipleSystemsNew(G, RHS, N_LocalDOFs, N_LocalDOFs, 1, 1);

    for(k=0;k<N_LocalDOFs;k++)
    {
      l = DOF[k];
      N_IncidentArray[l]++;
      //       OwnDof[l] = TRUE;
      Sol[l] += RHS[k];
    }
  }                                               // for(i=0;i<N_Cells;i++)

  for(k=0;k<N_U;k++)
  {
    //     if(N_IncidentArray[k]==0 && OwnDof[k])
    if(N_IncidentArray[k]==0)
    {
      cout<< "error in GetSolFromQuadPtVales " <<endl;
      exit(0);
    }

    Sol[k] /= double(N_IncidentArray[k]);
  }                                               // for(k=0;k<N_LocalDOFs;k++)


} //void GetSolFromQuadPtVales(


void GetNodalPtsADI_System(int &N_PBE_SpatialPts, TFESpace1D *FeSpace_Intl,
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

  bool IsIsoparametric, update;

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
  N_Cells = Coll->GetN_Cells();

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

  N_PBE_SpatialPts = N_RootPts;

  //    for(i=0; i<N_PBE_SpatialPts; i++)
  //        cout <<i<< " IntlX[i] " <<IntlX[i]<< " IntlY[i] " <<IntlY[i] <<endl;
  //
  // cout<< "N_RootPts " << N_RootPts <<endl;
}


void GetPtsValuesForIntl(TFEFunction2D *Heat, TFEFunction2D *Concentration, TFESpace2D *PBE_Spaces,
double *Values, double *TValues, double *Sat_Values, int N_PBE_SpatialPts)
{
  int i, j, k, l, m, N_Cells, N_Points, N_LocalDOFs, T_N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers, *DOF;
  int *T_BeginIndex, *T_GlobalNumbers, *T_DOF;
  int PolynomialDegree, ApproxOrder, *IntlPtIndexOfPts, *PtIndexLoc, disp;
  int *IncidentArray;

  double *xi, *eta, *C_NodalValues, s, *T_NodalValues;
  double BasisValues[MaxN_BaseFunctions2D], T_val;
  double T_ref = TDatabase::ParamDB->REACTOR_P23;
  double C_ref = TDatabase::ParamDB->REACTOR_P25;

  TFESpace2D *fespace2D, *T_fespace2D;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId, T_FEId, FEId_PBE;
  TFE2D *Element, *T_Element, *Element_PBE;
  TBaseFunct2D *bf, *T_bf;
  TNodalFunctional2D *nf_PBE;

  T_NodalValues =  Heat->GetValues();
  T_fespace2D = Heat->GetFESpace2D();
  T_BeginIndex = T_fespace2D->GetBeginIndex();
  T_GlobalNumbers = T_fespace2D->GetGlobalNumbers();

  C_NodalValues =  Concentration->GetValues();
  fespace2D = Concentration->GetFESpace2D();
  BeginIndex = fespace2D->GetBeginIndex();
  GlobalNumbers = fespace2D->GetGlobalNumbers();

  IntlPtIndexOfPts = PBE_Spaces->GetIntlPtIndexOfPts();
  IncidentArray = new int[N_PBE_SpatialPts];

  memset(Values, 0, N_PBE_SpatialPts*SizeOfDouble);
  memset(TValues, 0, N_PBE_SpatialPts*SizeOfDouble);
  memset(Sat_Values, 0, N_PBE_SpatialPts*SizeOfDouble);
  memset(IncidentArray, 0, N_PBE_SpatialPts*SizeOfInt);

  // assume that all fespace2Ds and PBE_Spaces use same coll
  Coll = fespace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  disp = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId_PBE = PBE_Spaces->GetFE2D(i, cell);
    Element_PBE = TFEDatabase2D::GetFE2D(FEId_PBE);
    nf_PBE = Element_PBE->GetNodalFunctional2D();
    nf_PBE->GetPointsForAll(N_Points, xi, eta);
    PtIndexLoc = IntlPtIndexOfPts + disp;

    FEId = fespace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    bf = Element->GetBaseFunct2D();
    N_LocalDOFs = Element->GetN_DOF();
    DOF = GlobalNumbers + BeginIndex[i];

    T_FEId = T_fespace2D->GetFE2D(i, cell);
    T_Element = TFEDatabase2D::GetFE2D(T_FEId);
    T_bf = T_Element->GetBaseFunct2D();
    T_N_LocalDOFs = T_Element->GetN_DOF();
    T_DOF = T_GlobalNumbers + T_BeginIndex[i];

    for(j=0;j<N_Points;j++)
    {
      k = PtIndexLoc[j];
      IncidentArray[k]++;

      // C value at this point
      bf->GetDerivatives(D00, xi[j], eta[j], BasisValues);
      for(l=0;l<N_LocalDOFs;l++)
      {
        m = DOF[l];
        s = C_NodalValues[m];
        Values[k] += BasisValues[l]*s;
      }

      // first find T value at this point, then compute C_sat
      T_bf->GetDerivatives(D00, xi[j], eta[j], BasisValues);
      for(l=0;l<T_N_LocalDOFs;l++)
      {
        m = T_DOF[l];
        s = T_NodalValues[m];
        TValues[k] += BasisValues[l]*s;
      }

    }                                             // for(j=0;j<N_Points;j++)

    disp +=N_Points;
  }                                               // for(i=0;i<N_Cells;i++)

  for(i=0;i<N_PBE_SpatialPts;i++)
  {
    if(IncidentArray[i]==0)
    {
      cout << i<<  " IncidentArray[i] " << IncidentArray[i] <<endl;
      exit(0);
    }
    Values[i] /=double(IncidentArray[i]);
    TValues[i] /=double(IncidentArray[i]);        // contains T value

    // now compute C_Sat from T
    Sat_Values[i] = (1.3045*(T_ref*TValues[i] - 273.15) + 35.3642)/C_ref;

  }
  // cout << " GetPtsValuesForIntl " <<endl;
  // MPI_Finalize;
  // exit(0);
  delete [] IncidentArray;
}


void GetVeloGradForIntlPts(TFESpace2D *velospace2D, TFEFunction2D *velo1, TFEFunction2D *velo2,
TFESpace2D *PBE_Spaces, double *Velo_PBE_SpatialPts, double *GradVelo_PBE_SpatialPts,
int N_PBE_SpatialPts)
{
  int i, j, k, l, m, N_Cells, N_Points, N_LocalDOFs, T_N_LocalDOFs, N;
  int *BeginIndex, *GlobalNumbers, *DOF;
  int *T_BeginIndex, *T_GlobalNumbers, *T_DOF;
  int PolynomialDegree, ApproxOrder, *IntlPtIndexOfPts, *PtIndexLoc, disp;
  int *IncidentArray;

  double *xi, *eta, *C_NodalValues, s, *T_NodalValues;
  double BasisValues[MaxN_BaseFunctions2D], BasisValuesX[MaxN_BaseFunctions2D], T_val;
  double BasisValuesY[MaxN_BaseFunctions2D], *u1_NodalValues, *u2_NodalValues, u1, u2;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  
  TFESpace2D *fespace2D, *T_fespace2D;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId, T_FEId, FEId_PBE;
  TFE2D *Element, *T_Element, *Element_PBE;
  TBaseFunct2D *bf, *T_bf;
  TNodalFunctional2D *nf_PBE;
  BaseFunct2D *BaseFuncts;
  RefTrans2D RefTrans;
  
  IntlPtIndexOfPts = PBE_Spaces->GetIntlPtIndexOfPts();
  IncidentArray = new int[N_PBE_SpatialPts];

  BeginIndex = velospace2D->GetBeginIndex();
  GlobalNumbers = velospace2D->GetGlobalNumbers();
  u1_NodalValues =  velo1->GetValues();
  u2_NodalValues =  velo2->GetValues();
  
  memset(IncidentArray, 0, N_PBE_SpatialPts*SizeOfInt);
  memset(Velo_PBE_SpatialPts, 0, 2*N_PBE_SpatialPts*SizeOfDouble);
  
// #ifdef __AXIAL3D__     
//   memset(GradVelo_PBE_SpatialPts, 0, 9*N_PBE_SpatialPts*SizeOfDouble);  // radial velo is zero
// #else 
  memset(GradVelo_PBE_SpatialPts, 0, 4*N_PBE_SpatialPts*SizeOfDouble);
// #endif  

  
   BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
     
  // assume that all fespace2Ds and PBE_Spaces use same coll
  Coll = velospace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  disp = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId_PBE = PBE_Spaces->GetFE2D(i, cell);
    Element_PBE = TFEDatabase2D::GetFE2D(FEId_PBE);
    nf_PBE = Element_PBE->GetNodalFunctional2D();
    nf_PBE->GetPointsForAll(N_Points, xi, eta);
    PtIndexLoc = IntlPtIndexOfPts + disp;

    FEId = velospace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    RefTrans = Element->GetRefTransID();    
  
    // set cell for reference transformation
    TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);    
    
    bf = Element->GetBaseFunct2D();
    N_LocalDOFs = Element->GetN_DOF();
    DOF = GlobalNumbers + BeginIndex[i];
    
  
    for(j=0;j<N_Points;j++)
    {
      k = PtIndexLoc[j];
      IncidentArray[k]++;

      // u value at this point
      bf->GetDerivatives(D00, xi[j], eta[j], BasisValues);
      bf->GetDerivatives(D10, xi[j], eta[j], BasisValuesX);
      bf->GetDerivatives(D01, xi[j], eta[j], BasisValuesY);
      
      TFEDatabase2D::GetOrigValues(RefTrans, xi[j], eta[j], bf, Coll, (TGridCell *)cell,
                                   BasisValues, BasisValuesX, BasisValuesY, uorig, uxorig, uyorig);        
      
      
      
      for(l=0;l<N_LocalDOFs;l++)
      {
        m = DOF[l];
        u1 = u1_NodalValues[m];
        u2 = u2_NodalValues[m];

        // u1 and u2
        Velo_PBE_SpatialPts[2*k] += uorig[l]*u1;
        Velo_PBE_SpatialPts[2*k + 1] += uorig[l]*u2;

        // u1grad and u2grad
// #ifdef __AXIAL3D__   
//         GradVelo_PBE_SpatialPts[9*k] += BasisValuesX[l]*u1;   radia velo is zero, so still 2d arrays can be used in pipe flow
// //         GradVelo_PBE_SpatialPts[9*k + 1] += 0;
//         GradVelo_PBE_SpatialPts[9*k + 2] += BasisValuesY[l]*u1;
// 
// //         GradVelo_PBE_SpatialPts[9*k + 3] += 0;
//         GradVelo_PBE_SpatialPts[9*k + 4] += BasisValues[l]*u1/r_axial;
// //         GradVelo_PBE_SpatialPts[9*k + 5] += 0;
// 	
//         GradVelo_PBE_SpatialPts[9*k + 6] += BasisValuesX[l]*u2;
// //         GradVelo_PBE_SpatialPts[9*k + 7] += 0;
//         GradVelo_PBE_SpatialPts[9*k + 8] += BasisValuesY[l]*u2;
// #else
        GradVelo_PBE_SpatialPts[4*k] += uxorig[l]*u1;
        GradVelo_PBE_SpatialPts[4*k + 1] += uyorig[l]*u1;
        GradVelo_PBE_SpatialPts[4*k + 2] += uxorig[l]*u2;
        GradVelo_PBE_SpatialPts[4*k + 3] += uyorig[l]*u2;
// #endif
      }                                           //for(l=0;l<N_LocalDOFs;l++)
    }                                             // for(j=0;j<N_Points;j++)
    disp +=N_Points;
  }                                               // for(i=0;i<N_Cells;i++)

   OutPut("Velo_PBE_SpatialPts " << Ddot(2*N_PBE_SpatialPts, Velo_PBE_SpatialPts, Velo_PBE_SpatialPts ) << endl);

  for(i=0;i<N_PBE_SpatialPts;i++)
  {
    s = (double)IncidentArray[i];
    if(s==0.0)
    {
      cout << i<<  " IncidentArray[i] " << IncidentArray[i] <<endl;
      exit(0);
    }

    Velo_PBE_SpatialPts[2*i] /=s;
    Velo_PBE_SpatialPts[2*i + 1] /=s;

// #ifdef __AXIAL3D__   
//     N=9;
// #else
    N=4;
// #endif

    GradVelo_PBE_SpatialPts[4*i ] /=s; 
    GradVelo_PBE_SpatialPts[4*i + 1] /=s;
    GradVelo_PBE_SpatialPts[4*i + 2] /=s;
    GradVelo_PBE_SpatialPts[4*i + 3] /=s;

//   exit(0);
  }
}

 /** transfer sol from spatial to internal */
void GetPtsValuesForIntl(int N_Levels, TFEFunction2D **ScalarFunctions, double *Sol, int N_U, double *ValuesT, int N_PBE_SpatialPts)
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
  IncidentArray = new int[N_Levels*N_PBE_SpatialPts];

  memset(ValuesT, 0, N_Levels*N_PBE_SpatialPts*SizeOfDouble);
  memset(IncidentArray, 0, N_Levels*N_PBE_SpatialPts*SizeOfInt);

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

  for(i=0;i<N_PBE_SpatialPts;i++)
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

void AssembleFD_M_StartPtDirichlet(TSquareMatrix1D *M)
{
  int i, j, N, N_Values; 
  int begin, end, *Row, *KCol;
  double *ValuesM;
  
  
  N = M->GetN_Columns();
  Row = M->GetRowPtr();
  KCol = M->GetKCol();
  ValuesM = M->GetEntries();  
   
  M->Reset();
    
    for(i=0;i<N;i++)
    {
      begin=Row[i];
      end=Row[i+1];
      for(j=begin;j<end;j++)
      {
       if(i==KCol[j])
        ValuesM[j] = 1.;
      }   // endfor j
    }    
  
}


void AssembleM_StartPtDirichlet(TFESpace1D *FeSpace, TSquareMatrix1D *M)
{
  int i, j, k, l, N_Cells, N_BaseFunct, N_U, dGDisc;
  int N_Points, N_Sets=1, *GlobalNumbers, *BeginIndex, *DOF;
  int TestDOF, begin, end, *RowPtr, *KCol;

  double *Weights, *zeta, X[20], AbsDetjk[20];
  double LocMatrixM[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
  double **origvaluesD0, **origvaluesD1, Mult;
  double *orgD0, *orgD1, test0, test1, ansatz0, ansatz1, *ValuesM;
  double x=0., len =0.;

  bool Needs2ndDer[1];

  TBaseCell *Cell;
  FE1D FEId;
  TFE1D *Element;
  TBaseFunct1D *bf;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TRefTrans1D *F_K;
  BaseFunct1D BaseFunct_ID, BaseFunct[1];
  TCollection *Coll;
  BoundCond BDType;
  BoundCond cond_Lmin, cond_Lmax;

  BoundCondition_LminLMax(cond_Lmin, cond_Lmax);

  dGDisc=FeSpace->IsDGSpace();
  Coll = FeSpace->GetCollection();
  N_U = FeSpace->GetN_DegreesOfFreedom();
  GlobalNumbers = FeSpace->GetGlobalNumbers();
  BeginIndex = FeSpace->GetBeginIndex();
  RowPtr = M->GetRowPtr();
  KCol = M->GetKCol();
  ValuesM = M->GetEntries();

  // all QuadPts in a cell use same FEspace in internal direction
  N_Cells = Coll->GetN_Cells();
  Needs2ndDer[0] = FALSE;

  for(i=0; i<N_Cells; i++)
  {
    Cell = Coll->GetCell(i);
    FEId = FeSpace->GetFE1D(i, Cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();

    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_Points, Weights, zeta);

    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(Cell);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, X, AbsDetjk);

    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta,  LineQuadFormula,  Needs2ndDer);

    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);

    memset(LocMatrixM, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);

    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Points;j++)
    {
      Mult = Weights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];

      //cout<< " zeta[j] " << zeta[j]  <<endl;
      //len +=Mult;
      for(k=0;k<N_BaseFunct;k++)
      {
        test0  = orgD0[k];
        //cout<< " uref " << test0  <<endl;
        for(l=0;l<N_BaseFunct;l++)
        {
          ansatz0  = orgD0[l];
          LocMatrixM[k*N_BaseFunct + l] += (Mult*ansatz0*test0);
        }
      }
    }

    //   add to global matrices
    for(j=0;j<N_BaseFunct;j++)
    {
      TestDOF = DOF[j];

      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(k=begin;k<end;k++)
      {
        for(l=0;l<N_BaseFunct;l++)
        {
          if(KCol[k] == DOF[l])
          {
            ValuesM[k] +=LocMatrixM[j*N_BaseFunct + l];
            break;
          }
        }                                         // for(m=0;m<N_BaseFunct_low
      }                                           // for(n=begin;n<end;n++)
    }                                             // for(l=0;l<N_BaseFunct_low
  }                                               // for(i=0; i<N_Cells; i++)

  //update boundary data
  // starting point:  in PBS
  // end point: in PBS

  if(cond_Lmin==DIRICHLET && !dGDisc)
  {
    begin = RowPtr[0];
    end = RowPtr[1];

    for(k=begin;k<end;k++)
    {
      if(KCol[k] == 0 )
        { ValuesM[k] = 1.; }
        else
          { ValuesM[k] = 0.; }
    }
  }                                             //if(cond_Lmin==DIRICHLET)

  if(cond_Lmax==DIRICHLET && !dGDisc)
  {
    begin = RowPtr[N_U-1];
    end = RowPtr[N_U];

    for(k=begin;k<end;k++)
    {
      if(KCol[k] == N_U-1 )
        { ValuesM[k] = 1.; }
        else
          { ValuesM[k] = 0.; }
    }
  }                                             //if(cond_Lmin==DIRICHLET)

  // cout<< " len " << len << endl;
  // exit(0);
  // //print matrix
  //   for(j=0;j<N_U;j++)
  //    {
  //     begin = RowPtr[j];
  //     end = RowPtr[j+1];
  //     for(k=begin;k<end;k++)
  //      {
  //       cout << "M(" << j << ", "<< KCol[k] << ") = " << ValuesM[k] <<endl;
  //      }
  //     cout<<endl;
  //    }
  //    exit(0);
}

/** transfer sol from internal to space */
void GetSolFromNodalPtVales(int N_Levels, int MaxN_PtsForNodal, TFEFunction2D **ScalarFunctions, double *Sol, int N_U, double *ValuesT, int N_PBE_SpatialPts)
{
  int i, ii, j, k, l, m, N_Cells, N_Points, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers, *DOF, disp;
  int PolynomialDegree, ApproxOrder, *IntlPtIndexOfPts, *PtIndexLoc;

  double *xi, *eta, *NodalValues, s, *Values_Level, *Sol_Level;
  double *val, *PtValues, maxval=0.;
  double FunctionalValues[MaxN_PointsForNodal2D];

  TFESpace2D *PBE_Spaces;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId_PBE;
  TFE2D *Element_PBE;
  TNodalFunctional2D *nf_PBE;

  // assume all  ScalarFunctions use same fespace2D
  PBE_Spaces = ScalarFunctions[0]->GetFESpace2D();
  BeginIndex = PBE_Spaces->GetBeginIndex();
  GlobalNumbers = PBE_Spaces->GetGlobalNumbers();
  IntlPtIndexOfPts = PBE_Spaces->GetIntlPtIndexOfPts();
  val = new double[MaxN_PtsForNodal*N_Levels];

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
    DOF = GlobalNumbers + BeginIndex[i];
    N_LocalDOFs = Element_PBE->GetN_DOF();
    PtIndexLoc = IntlPtIndexOfPts + disp;

    // collect point values for all level
    for(j=0;j<N_Points;j++)
    {
      k = PtIndexLoc[j];
      Values_Level  =  ValuesT + k*N_Levels;

      for(ii=0;ii<N_Levels;ii++)
      {
        val[ii*N_Points + j] = Values_Level[ii];
        if(maxval<Values_Level[ii]) maxval=Values_Level[ii];
      }
    }                                             // for(j=0;j<N_Points;j++)

    for(ii=0;ii<N_Levels;ii++)
    {
      PtValues = val + ii*N_Points;
      Sol_Level = Sol+ ii*N_U;

      nf_PBE->GetAllFunctionals(Coll, (TGridCell *)cell, PtValues, FunctionalValues);

      for(j=0;j<N_LocalDOFs;j++)
        Sol_Level[DOF[j]] = FunctionalValues[j];
    }

    disp +=N_Points;
  }                                               // for(i=0;i<N_Cells;i++)

  delete [] val;
  //    cout << " GetSolFromQuadPtVales " << maxval << endl;
}    //void GetSolFromNodalPtVales(


void GetSolFromNodalPtVales(TFEFunction2D *ScalarFunction, double *Values, TFESpace2D *PBE_Spaces)
{
  int i, ii, j, k, l, m, N_Cells, N_Points, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers, *DOF, disp, *IntlPtIndexOfPts, *PtIndexLoc;

  double *xi, *eta;
  double val[MaxN_PointsForNodal2D], PtValues[MaxN_PointsForNodal2D], *Sol;
  double FunctionalValues[MaxN_PointsForNodal2D];

  TFESpace2D *FE_Space;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId;
  TFE2D *Element;
  TNodalFunctional2D *nf;

  // assume all  ScalarFunctions use same fespace2D
  Sol =  ScalarFunction->GetValues();
  FE_Space = ScalarFunction->GetFESpace2D();
  BeginIndex = FE_Space->GetBeginIndex();
  GlobalNumbers = FE_Space->GetGlobalNumbers();

  IntlPtIndexOfPts = PBE_Spaces->GetIntlPtIndexOfPts();

  // assume that both fespace2D and FE_Space use same coll
  Coll = FE_Space->GetCollection();
  N_Cells = Coll->GetN_Cells();

  disp = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId = FE_Space->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    nf = Element->GetNodalFunctional2D();
    nf->GetPointsForAll(N_Points, xi, eta);
    DOF = GlobalNumbers + BeginIndex[i];
    N_LocalDOFs = Element->GetN_DOF();

    PtIndexLoc = IntlPtIndexOfPts+disp;

    for(j=0;j<N_LocalDOFs;j++)
    {
      k = PtIndexLoc[j];
      PtValues[j] = Values[k];
    }

    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PtValues, FunctionalValues);

    for(j=0;j<N_LocalDOFs;j++)
      Sol[DOF[j]] = FunctionalValues[j];

    disp +=N_Points;
  }                                               // for(i=0;i<N_Cells;i++)

}                                                 //void GetSolFromNodalPtVales(

void  L_Nodal2Sol(TFESpace1D *FESpace1D, double *Sol_QuadIntl, int N_XLocPoints, 
                  int N_Levels, double *Sol_AllL)
{

  int i,j,k,l;
  int N_Cells;
  int N_DOFs, N_LocalDOFs, *IncidentArray;
  int *BeginIndex, *GlobalNumbers;
  int N_Points, ii;
  int *DOF, *Index, *IndexOfNodalPts, disp;

  double *xi, *eta, *Values_Level, *sol;
  double X[MaxN_PointsForNodal1D], Y[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double *PointValues, *PtVal;
  double FunctionalValues[MaxN_PointsForNodal1D];
  double FctVal[4];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TNodalFunctional1D *nf;

  Coll = FESpace1D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  IndexOfNodalPts = FESpace1D->GetIntlPtIndexOfPts();

  BeginIndex = FESpace1D->GetBeginIndex();
  GlobalNumbers = FESpace1D->GetGlobalNumbers();
  N_DOFs = FESpace1D->GetN_DegreesOfFreedom();

  PointValues = new double [MaxN_PointsForNodal1D*N_XLocPoints];
  IncidentArray = new int [N_DOFs];
  memset(Sol_AllL, 0, SizeOfDouble*N_DOFs*N_XLocPoints);
  memset(IncidentArray, 0, SizeOfInt*N_DOFs);

  disp = 0; 

  for(i=0;i<N_Cells;i++)
   {
    cell = Coll->GetCell(i);
    FEId = FESpace1D->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    DOF = GlobalNumbers+BeginIndex[i];

    Index = IndexOfNodalPts+disp;

    for(ii=0;ii<N_XLocPoints;ii++)
     {
      Values_Level  =  Sol_QuadIntl + ii*N_Levels;

      for(j=0;j<N_Points;j++)
      {
       k = Index[j]; // corresponding L level
       PointValues[ii*N_Points + j] = Values_Level[k];
      }
     }

    for(ii=0;ii<N_XLocPoints;ii++)
     {
      PtVal = PointValues + ii*N_Points;
      nf->GetAllFunctionals(PtVal, FunctionalValues);

      sol = Sol_AllL + ii*N_DOFs;

      for(j=0;j<N_LocalDOFs;j++)
       {
        sol[DOF[j]] += FunctionalValues[j];

        if(ii==0)
         IncidentArray[DOF[j]] ++;
       }
     }

    disp +=N_Points;
  } // for(i=0;i<N_Cells

  for(ii=0;ii<N_XLocPoints;ii++)
    {
     for(i=0;i<N_DOFs;i++)
      {
       if(ii==0)
        if(IncidentArray[i] == 0)
         {
          cout << "Error in L_Nodal2Sol : "<< IncidentArray[i] << endl;
          exit(0);
         }
       Sol_AllL[ii*N_DOFs + i] /= (double)IncidentArray[i];
      } // for(i=0;i<N_DOFs;i
    } // for(ii=0;ii<N_XLocPoints;i

  delete [] PointValues;
  delete [] IncidentArray;
}


void  L_Sol2Nodal(TFESpace1D *FESpace1D, double *Sol_AllL, int N_Levels,  double *Sol_NodalPts, int N_XPoints)
{
  int i,j,k,l, ii;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points, disp;
  int *DOF, *IndexArray, *NodalPtIndex;

  double *xi, *eta, *sol, *sol_Nodal, val;
  double Z[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double BasisValues[MaxN_PointsForNodal1D][MaxN_BaseFunctions1D];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;

  Coll = FESpace1D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace1D->GetBeginIndex();
  GlobalNumbers = FESpace1D->GetGlobalNumbers();
  N_DOFs = FESpace1D->GetN_DegreesOfFreedom();
  NodalPtIndex = FESpace1D->GetIntlPtIndexOfPts();

  IndexArray = new int[N_Levels];
  memset(IndexArray, 0, SizeOfInt*N_Levels);
  memset(Sol_NodalPts , 0, SizeOfDouble*N_XPoints*N_Levels);

  disp = 0;
  for(i=0; i<N_Cells; i++)
   {
    cell = Coll->GetCell(i);
    FEId = FESpace1D->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);
    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, Z, AbsDetjk);
    DOF = GlobalNumbers + BeginIndex[i];

     for(j=0;j<N_Points;j++)
       bf->GetDerivatives(D0, xi[j], BasisValues[j]);

     for(ii=0;ii<N_XPoints;ii++)
      {
       sol = Sol_AllL + ii*N_DOFs;
       sol_Nodal =  Sol_NodalPts + ii*N_Levels; 

//        for(l=0;l<N_LocalDOFs;l++)
//         cout << " val " << sol[DOF[l]] << endl;

       for(j=0;j<N_Points;j++)
        {
         k = NodalPtIndex[disp + j]; // find the right L level
         val = 0.;

          for(l=0;l<N_LocalDOFs;l++)
           val += sol[DOF[l]]*BasisValues[j][l];

         sol_Nodal[k] += val;
         if(ii==0)
          IndexArray[k]++;
        } // for(j=0;j<N_Points
      } // for(ii=0

     disp +=N_Points;
   } // for(i=0; i<N_Cells; i


   for(ii=0;ii<N_XPoints;ii++)
    {
     sol_Nodal =  Sol_NodalPts + ii*N_Levels; 
     for(i=0;i<N_Levels;i++)
      {
       if(ii==0)
        if(IndexArray[i] == 0)
         {
          cout << "Error in L_Sol2Nodal : "<< IndexArray[i] << endl;
          exit(0);
         }

       sol_Nodal[i] /= (double)IndexArray[i];
      }
    }
}




void GetOSError(int N_PBE_IntlLevels, TFEFunction2D **ScalarFunctions,  TFESpace1D *FESpace1D,
double *SolPbe, int N_U, double *errors, DoubleFunctND *ExactFunct)
{
  int i,j,k,l, m, n, N_Cells, N_V, N_LDof;
  int N_DOFs, N_LocalDOFs, N_BaseFunct_Intl;
  int *BeginIndex, *GlobalNumbers, *DOF, *BeginIndex_Intl, *GlobalNumbers_Intl, *DOF_Intl;
  int N_LocalUsedElements, N_Points, *N_BaseFunct, N_Cells_Intl;
  int L, N_LinePoints, N_Sets=1;

  double *weights, *xi, *eta, *Sol_AllL;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double *Sol_QuadIntl, *sol;
  double **OrigFEValues, *Orig, value, *InternalError, *InternalErrorgrad;
  double *LineWeights, *zeta, Z[MaxN_QuadPoints_1D], LineAbsDetjk[MaxN_QuadPoints_1D];;
  double **origvaluesD0, **origvaluesD1, *orgD0, *orgD1, Mult;
  double x, y, z, Exact[5], *sol_QI, valuegrad, Xi[3];

  TBaseCell *cell, *Cell_Intl;
  TCollection *Coll, *Coll_Intl;
  TFESpace2D *FESpace2D;
  TFE2D *Element;
  FE2D LocalUsedElements[1], CurrentElement;
  BaseFunct2D BaseFunct, *BaseFuncts;
  FE1D FEId_Intl;
  TFE1D *Element_Intl;
  TBaseFunct1D *bf_Intl;
  BaseFunct1D BaseFunct_ID_Intl, BaseFunct_Intl[1];
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TRefTrans1D *F_K;

  bool *SecondDer;
  bool Needs2ndDer[1];

  FESpace2D = ScalarFunctions[0]->GetFESpace2D();
  Coll = FESpace2D->GetCollection();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  N_Cells = Coll->GetN_Cells();

  Coll_Intl = FESpace1D->GetCollection();
  BeginIndex_Intl = FESpace1D->GetBeginIndex();
  GlobalNumbers_Intl = FESpace1D->GetGlobalNumbers();
  N_Cells_Intl = Coll_Intl->GetN_Cells();
  N_LDof = FESpace1D->GetN_DegreesOfFreedom();
  N_LocalUsedElements = 1;
  SecondDer = new bool[1];
  SecondDer[0] = FALSE;
  Needs2ndDer[0] = FALSE;

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  //first check how many quad pts  in the cell
  cell = Coll->GetCell(0);
  LocalUsedElements[0] = FESpace2D->GetFE2D(0, cell);
  Element = TFEDatabase2D::GetFE2D(LocalUsedElements[0]);
  TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, Coll, cell, SecondDer,
    N_Points, xi, eta, weights, X, Y, AbsDetjk);

  Sol_QuadIntl = new double[N_Points*N_PBE_IntlLevels];
  Sol_AllL = new double[N_Points*N_LDof];
  InternalError = new double[N_Points];
  InternalErrorgrad = new double[N_Points];

  errors[0] = 0.;
  errors[1] = 0.;

  for(i=0; i<N_Cells; i++)
  {
    cell = Coll->GetCell(i);
    LocalUsedElements[0] = FESpace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(LocalUsedElements[0]);

  // ####################################################################
  // calculate values on original element
  // ####################################################################
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, Coll, cell, SecondDer,
      N_Points, xi, eta, weights, X, Y, AbsDetjk);
    // calculate all needed derivatives of this FE function
    CurrentElement = FESpace2D->GetFE2D(i, cell);
    BaseFunct = BaseFuncts[CurrentElement];
    N_LocalDOFs = N_BaseFunct[CurrentElement];
    DOF = GlobalNumbers+BeginIndex[i];

    OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);

    // find values at all quad points
    memset(Sol_QuadIntl, 0, N_PBE_IntlLevels*N_Points*SizeOfDouble);
    for(j=0; j<N_PBE_IntlLevels; j++)
    {
      sol = SolPbe+j*N_U;

      for(k=0; k<N_Points; k++)
      {
        Orig = OrigFEValues[k];
        value = 0.;

        for(l=0; l<N_LocalDOFs; l++)
          value += sol[DOF[l]]*Orig[l];

        Sol_QuadIntl[k*N_PBE_IntlLevels  +  j] = value;
      }
    }      //for(j=0; j<N_PBE_IntlLevels; j++)

   L_Nodal2Sol(FESpace1D, Sol_QuadIntl, N_Points, N_PBE_IntlLevels, Sol_AllL);

    memset(InternalError, 0, N_Points*SizeOfDouble);
    memset(InternalErrorgrad, 0, N_Points*SizeOfDouble);

    for(j=0; j<N_Points; j++)
    {
      x = X[j];
      y = Y[j];
//       sol_QI = Sol_QuadIntl+j*N_PBE_IntlLevels;
      sol_QI = Sol_AllL + j*N_LDof;

      for(k=0; k<N_Cells_Intl; k++)
      {
        Cell_Intl = Coll_Intl->GetCell(k);
        FEId_Intl = FESpace1D->GetFE1D(k, Cell_Intl);
        Element_Intl = TFEDatabase2D::GetFE1D(FEId_Intl);
        bf_Intl = Element_Intl->GetBaseFunct1D();
        N_BaseFunct_Intl = Element_Intl->GetN_DOF();
        BaseFunct_ID_Intl = Element_Intl->GetBaseFunct1D_ID();
        DOF_Intl = GlobalNumbers_Intl+BeginIndex_Intl[k];

        L = bf_Intl->GetPolynomialDegree();
        LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*L);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
        ((TLineAffin *)F_K)->SetCell(Cell_Intl);
        ((TLineAffin *)F_K)->GetOrigFromRef(N_LinePoints, zeta, Z, LineAbsDetjk);

        BaseFunct_Intl[0] = BaseFunct_ID_Intl;
        ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct_Intl, N_LinePoints, zeta,  LineQuadFormula,  Needs2ndDer);

        origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID_Intl, D0);
        origvaluesD1=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID_Intl, D1);
        DOF_Intl = GlobalNumbers_Intl + BeginIndex_Intl[k];

        for(l=0; l<N_LinePoints; l++)
        {
          Mult = 0.5*LineWeights[l]*LineAbsDetjk[l];
          orgD0 = origvaluesD0[l];
          orgD1 = origvaluesD1[l];
          z = Z[l];

          //find sol at this point
          value = 0.; valuegrad= 0.;
          for(m=0; m<N_BaseFunct_Intl; m++)
          {
            value += sol_QI[DOF_Intl[m]]*orgD0[m];
            valuegrad  += sol_QI[DOF_Intl[m]]*orgD1[m];
          }
          Xi[0] = x;
          Xi[1] = y;
	  Xi[2] = z;
          ExactFunct(3,  Xi, Exact);
          InternalError[j] += Mult*(Exact[0]-value)*(Exact[0]-value);
          InternalErrorgrad[j] += Mult*(Exact[1]-valuegrad)*(Exact[1]-valuegrad);

//           if(j==0 && l==0)
//            cout <<  j << " " << k << " " << Exact[0] << " " << value << endl;

        }  //for(l=0; l
      } //for(k=0;

      Mult = weights[j]*AbsDetjk[j];
      errors[0] += Mult*InternalError[j];
      errors[1] += Mult*InternalErrorgrad[j];
    }     // for(j=0;
  }      //for(i=0; i<N_Cells; i++)

  delete [] Sol_QuadIntl;
  delete [] InternalError;

}  //void GetOSError


void SetLMinLMaxBoundValue(int N_PBE_IntlLevels, TFEFunction2D **PbeFunctionsXdir, int N_U, double *RhsArray_Pbe)
{
  double *sol;
  BoundCond cond_Lmin, cond_Lmax;

  BoundCondition_LminLMax(cond_Lmin, cond_Lmax);

  //set the boundary condition at the L_Min
  if(cond_Lmin==NEUMANN)
  {
    // do nothing, it should be zero Neumann
  }
  else if(cond_Lmin==DIRICHLET)
  {
    PbeFunctionsXdir[0]->Interpolate(BoundValue_LMin);
    sol = PbeFunctionsXdir[0]->GetValues();
    memcpy(RhsArray_Pbe,  sol,  N_U*SizeOfDouble);

  }
  else
  {
    cout<< " Only DIRICHLET and zero NEUMANN are allowed " <<endl;
    exit(0);
  }

  //set the boundary condition at the L_Min
  if(cond_Lmax==NEUMANN)
  {
    // do nothing, it should be zero Neumann
  }
  else if(cond_Lmax==DIRICHLET)
  {
    PbeFunctionsXdir[N_PBE_IntlLevels-1]->Interpolate(BoundValue_LMax);
    sol = PbeFunctionsXdir[N_PBE_IntlLevels-1]->GetValues();
    memcpy(RhsArray_Pbe+((N_PBE_IntlLevels-1)*N_U),  sol,  N_U*SizeOfDouble);
  }
  else
  {
    cout<< " Only DIRICHLET and zero NEUMANN are allowed " <<endl;
    exit(0);
  }
}

void CompOutletAvgFunct(double *solpbe, int N_U, TFEFunction2D *PbeFunctions, 
			int *Out_CellIndex, int *Out_EdgeIndex, int N_OutEdges)
{
  int i, IJoint, k, l, N, N_LinePoints, TestDOF;
  int *BeginIndex, *GlobalNumbers, *DOF, N_BaseFunct, *N_BaseFuncts;
  
  double *Values1, *Values2, *Values3, *Values4, **uref, *uorig;
  double *LineWeights, *zeta, t0, t1, normn, S1, S2, S3, S4;
  double  X_B[100], Y_B[100], val, val1, val2, val3, val4;

  TFESpace2D *FeSpace;
  TCollection *Coll;
  TBaseCell *cell;
  FE2D FEId;
  BaseFunct2D *BaseFuncts;  
  TFE2D *ele;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  BF2DRefElements RefElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  
  FeSpace = PbeFunctions->GetFESpace2D();
  Coll = FeSpace->GetCollection();

  BeginIndex = FeSpace->GetBeginIndex();
  GlobalNumbers = FeSpace->GetGlobalNumbers();

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  

  Values1 = solpbe+6*N_U;
  Values2 = solpbe+13*N_U; 
  Values3 = solpbe+18*N_U; 
  Values4 = solpbe+4*N_U; 

  S1=0;
  S2=0;
  S3=0;
  S4=0;
  // outlet int
   for(i=0;i<N_OutEdges;i++)
    {
     N = Out_CellIndex[i];
     IJoint=Out_EdgeIndex[i];
     
     cell = Coll->GetCell(N);

     FEId = FeSpace->GetFE2D(N, cell);   
     DOF = GlobalNumbers + BeginIndex[N];
     N_BaseFunct = N_BaseFuncts[FEId];
     ele = TFEDatabase2D::GetFE2D(FEId);
     RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);

     l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
     LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
     qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
     qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
     TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);
     
      switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadAffin;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadAffin *)F_K)->SetCell(cell);
        break;

        case BFUnitTriangle:
          RefTrans = TriaAffin;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaAffin *)F_K)->SetCell(cell);
        break;
      } // endswitch
     

      uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId], LineQuadFormula, IJoint);
      
      for(k=0;k<N_LinePoints;k++)
       {
        F_K->GetTangent(IJoint, zeta[k], t0, t1);  
        normn = sqrt(t0*t0+t1*t1);
        uorig = uref[k];

	// find value at this quad point
	val1 = 0;  val2 = 0;  val3 = 0;
        for(l=0;l<N_BaseFunct;l++)
         {
          TestDOF = DOF[l];
	  val=uorig[l];
	  val1 += val*Values1[TestDOF];
	  val2 += val*Values2[TestDOF];
	  val3 += val*Values3[TestDOF];  
	  val4 += val*Values4[TestDOF];  
         }
         
         val=normn*LineWeights[k];
         S1 += val*val1;
         S2 += val*val2;         
         S3 += val*val3;    
         S4 += val*val3;    
       }  
    }// for(i=0;i<N_C
 
//    if(fabs(S1)<1.0e-100 || S1<0) S1 = 0;
//    if(fabs(S2)<1.0e-100 || S2<0) S2 = 0;
//    if(fabs(S3)<1.0e-100 || S3<0) S3 = 0;  
   OutPut(" S OutBdInt : j=4,6,13,18: " << S4 << "   " << S1 << "   " << S2 << "    " << S3 <<endl); 
}


void GetInternalNodalPts(TFESpace1D * FeSpace_Intl, int &N_PBE_IntlLevels, double *&IntlPosL)
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

  N_PBE_IntlLevels = N_RootPts;

//  cout << N_AllLocalPoints << " N_RootPts  "  << N_RootPts << endl;
//   for(i=0; i<N_RootPts; i++)
//   cout << i << " L: "  << IntlPosL[i] << endl;
// exit(0);
}



#endif




void Get_FeFunction2DMass(TFEFunction2D *fefunction, double *parameters)
 {
  int i,j,k,l, polydegree, N_QFPoints, ORDER;
  int N_Cells, N_Joints, N_Vertices;
  int *BeginIndex, *GlobalNumbers, *DOF, N_BF;

  double *U, Mult, r_axial, val, mass, volume, Concentration;
  double *weights, *xi, *eta;
  double values[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];
  double AbsDetjk[MaxN_QuadPoints_2D], X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];

  TJoint *joint;
  TFESpace2D *FeSpace;
  TBaseCell *cell;
  TCollection *coll;
  JointType jointtype;
  BoundTypes bdtype;
  RefTrans2D RefTrans;
  bool IsIsoparametric;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;
  FE2D FEid;
  TBaseFunct2D *bf;
  TRefTrans2D *F_K;

  FeSpace = fefunction->GetFESpace2D();
  BeginIndex = FeSpace->GetBeginIndex();
  GlobalNumbers = FeSpace->GetGlobalNumbers();
  U = fefunction->GetValues();

  coll = FeSpace->GetCollection();
  N_Cells = coll->GetN_Cells();

  mass = 0.;
  volume = 0.;
  Concentration = 0.;
  for(i=0;i<N_Cells;i++)
   {
    cell = coll->GetCell(i);
    FEid = FeSpace->GetFE2D(i, cell);

    RefTrans = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(FEid);
    N_Joints = cell->GetN_Joints();
    IsIsoparametric = FALSE;
    if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
     {
     for(j=0;j<N_Joints;j++)
      {
       joint = cell->GetJoint(j);
       jointtype = joint->GetType();
       if(jointtype == BoundaryEdge)
        {
         bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
         if(bdtype != Line)  IsIsoparametric = TRUE;
        }
       if(jointtype == InterfaceJoint)
        {
          bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
        if(jointtype == IsoInterfaceJoint || jointtype == IsoBoundEdge)
         IsIsoparametric = TRUE;

      } // for(j=0;j<
     } // if(TDatabase::ParamDB->USE_ISOPARAMETRIC)

   if(IsIsoparametric)
    {
      switch(N_Joints)
      {
        case 4:
          RefTrans = QuadIsoparametric;
        break;

        case 3:
          RefTrans = TriaIsoparametric;
        break;
      }
    } // endif IsIsoparametric

    F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
    ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEid);
    switch(RefTrans)
    {
      case TriaAffin:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(9);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TTriaAffin *)F_K)->SetCell(cell);
//         locvol = ((TTriaAffin *)rt)->GetVolume();
        ((TTriaAffin *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta, X, Y, AbsDetjk);
      break;

      case TriaIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(9);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
        ((TTriaIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
        ((TTriaIsoparametric *)F_K)->SetCell(cell);
//         locvol = ((TTriaIsoparametric *)F_K)->GetVolume();
        ((TTriaIsoparametric *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case QuadAffin:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TQuadAffin *)F_K)->SetCell(cell);
//         locvol = ((TQuadAffin *)rt)->GetVolume();
        ((TQuadAffin *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta, X, Y, AbsDetjk);
      break;

      case QuadBilinear:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TTriaIsoparametric *)F_K)->SetApproximationOrder(polydegree);
        ((TQuadBilinear *)F_K)->SetCell(cell);
        ((TQuadBilinear *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta, X, Y, AbsDetjk);
      break;

      case QuadIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TQuadIsoparametric *)F_K)->SetApproximationOrder(polydegree);
        ((TQuadIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
        ((TQuadIsoparametric *)F_K)->SetCell(cell);
//         locvol = ((TQuadIsoparametric *)rt)->GetVolume();
        ((TQuadIsoparametric *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta, X, Y, AbsDetjk);
      break;
    }

    // find basis functions on cell i
    bf = TFEDatabase2D::GetBaseFunct2DFromFE2D(FEid);
    N_BF = bf->GetDimension();
    DOF = GlobalNumbers + BeginIndex[i];

    for(k=0;k<N_QFPoints;k++)
     {
      bf->GetDerivatives(D00, xi[k], eta[k], values[k]);
//       r_axial = fabs(X[k]);
      
     if(TDatabase::ParamDB->Axial3DAxis==1)
      {
       r_axial = Y[k];  // r in axial3D (X: symmetric) problems   
      }
     else
      {
       r_axial = X[k];  // r in axial3D (Y: symmetric) problems     
      }     
      
      if(r_axial<=0)
       {
        cout << k << " xi, eta: " << xi[k] << " "  << eta[k] << endl;
        cout << k << " xi, eta: " << X[k] << " "  << eta[k] << endl;
        cout << " r_axial negative in Get_FeFunction2DMass change Quad rule " <<  r_axial <<endl;
        exit(0);
       }

      Mult = r_axial*weights[k]*AbsDetjk[k];
      val = 0.;
      for(l=0;l<N_BF;l++)
       {
        j = DOF[l];
        val += U[j]*values[k][l];
       }

     mass += val*Mult;
     volume += Mult;
    } //  for(k=0;k<N_QFPoints;
   } //  for(i=0;i<N_Cells;i++)

   parameters[0] = 2.*Pi*mass;
   parameters[1] = 2.*Pi*volume;
   Concentration = mass/volume;
   parameters[2] = Concentration;


//    OutPut( "Time, C Mass " <<TDatabase::TimeDB->CURRENTTIME<< " " <<parameters[0]<< " "<<endl);
//    OutPut( "Time, C Surface area " <<TDatabase::TimeDB->CURRENTTIME<< " " <<parameters[1]<< " "<<endl);
// exit(0);
 }


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
  TFESpace2D **Scalar_Spaces, **Scalar_Spaces_low, *Velocity_Spaces, *pressure_space, *Grid_space;
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
  TSquareStructure2D **sqstructureA;
  TSquareMatrix2D *sqmatrixA, *SQMATRICES[3], *SQMATRICESSOLVER[3];
  TSquareMatrix2D *sqmatrixM , *sqmatrixK, *sqmatrixS;
  TSquareMatrix2D **MatricesA, **MatricesM, **MatricesK;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TFEVectFunct2D *velocity, *GridVelo, *GridPos, *Pbe_xdir, *Pbe_xdir_low;
  TFEFunction2D **ScalarFunctions, *velo1, *velo2, *fefct[3], *gridvelo1, *gridvelo2,  **ScalarFunctions_low;
  TFEFunction2D **ScalarPbeFunctions, *PbeFunctions, **ScalarPbeFunctions_low, *PSD_mass;
  TFESpace2D *fesp[2], *ferhs[1];
  TFESpace1D *FeSpace_Intl;
  TSquareStructure1D *SqStruct1D_Intl;
  TSquareMatrix1D *M_Intl, *A_Intl, *S_Intl, *K_Intl;
  TADISystem1D **ADI_System;
  TFEFunction1D *FeFunction_Intl;
  BoundCond cond_Lmin, cond_Lmax;
  TBoundEdge *Outlet_Joint;
  TBoundComp *BoundComp;
  TJoint *Joint;
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
  
#ifdef _OMP
  TParDirectSolver **Scalar_SMPSolver;
#endif

  double x,y,max,min,sum, L, L0, L1, *IntlPosL;
  double tau1, tau2, tau_array[4];
  double errors[9], olderror=0., olderror1=0., p1, p2, L2error_Max=-1.e8, L2error_Max_t;
  double t1, t2, res, res2, oldres, solver_time, solver_time_curr, residual, oldresidual;
  double total_time, t3, t4, **rhs_edge, **oldrhs_fem_fct0, **oldrhs_fem_fct1;
  double aggr_time=0, pbs_t1, pbs_t2, pbs_t3, pbs_t4, pbs_step1=0, pbs_step2=0;
  double impuls_residual,linredfac, values[5];
  double gamma, tau, oldtau, maxval=0.;
  double Parameters[2], hmin, hmax, limit;
  double **RhsArray, **SolArray, **SolArray_low, **B, *B_Pbe, *RhsArray_Pbe, *OldRhsArray_Pbe, *UArray;
  double *SolPbe, *IntlX, *IntlY, *SolPbe_low, *PbeAux;
  double *RHSs[3], *defect, end_time, *BDUpdatedSolPbe_Intl_Loc, C, C_Sat;
  double *Scalar_Entries, l2, H1, *PBE_IntlPtValuesT, *RhsPbe_Intl_Loc, *OldPBE_IntlPtValuesT, *C_IntlPtValues;
  double *PBE_LSolT;
  double *PBE_IntlPtValues, *PBE_IntlPtRhsValues, *T_IntlPtValues, *Velo_PBE_SpatialPts, *GradVelo_PBE_SpatialPts;
  double len, *C_Sat_IntlPtValues, *SolPbe_OneLevel, *PBE_IntlPtRhsValuesT, *OldPBE_IntlPtRhsValuesT;
  double *SolPbe_Lmin, *SolPbe_Lmax, q3_max, T;
  double *Sol_IntlLoc, *OldSol_IntlLoc, *B_IntlLoc, *defect_IntlLoc, *Sol_Loc, *SolPbe_Output, *MatValues, *Oldrhs_IntlLoc;
  double *WArray, *gridpos, *gridpos_old, *layers_mass, params[8], *pbeq3mass, *pbeq3massAtNodalPts, param[3];
  double l_infty, u_infty, f_infty; 

  double **MatValues_Solver;
  int **KCol_Solver, **Row_Solver;
  void **Symbolic_Solver, **Numeric_Solver;  
    
 
  int i, j, k, l, m, N, ret, img=1, N_Cells, N_PBE_IntlLevels, N_U, N_U_low, N_U_Max, N_LDof, ii;
  int *N_Uarray, *N_Active, *N_Uarray_low, *N_Active_low, N_Veloarray, PBE_INDEX, PBE_INDEX_LOW;
  int N_ScalarEqns, N_PBEqns, N_IndepntScalarEqns;
  int sold_parameter_type, VeloFunction, MaxN_PtsForNodal;
  int N_SquareMatrices, N_Rhs, N_FESpaces;
  int FirstSolve, Max_It, ORDER, PBE_ORDER, VELOCITYORDER, N_SubSteps;
  int N_Paramters=1, time_discs, N_ScalarSpaces, pressure_space_code;
  int very_first_time=0, Max_N_Unknowns=-1;
  int *ScalarKCol, *ScalarRowPtr, N_Entries, N_Eqn, *KCol, *RowPtr;
  int Eq_pos, begin, end, N_Cells_Intl, N_IntlPts_All, N_PBE_SpatialPts=0, N_IntlLevels=0;
  int Out_Level, Disctypes[4], start_pbe, end_pbe, N_Q3Data=1, dGDisc=0, N_OutData=1;
  int N_Edges, comp, N_OutEdges, *Out_CellIndex, *Out_EdgeIndex, f_Integ_NodalPt[5],  f_Integ_NodalPt_rank=-1;
  int *GlobalNumbers_Grid, *BeginIndex_Grid, N_GridDOFs, N_f_OutletIndex, *f_OutletIndex;
  int rb_flag=0, rb_flag_pbe=0, f_Inlet_SpaceIndex=-1;
  
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
  const char Q3dir[] = "Q3Data";
   const char vtkdir[] = "VTK"; 
  std::ostringstream os;
  os << " ";

  mkdir(Q3dir, 0777);  
  mkdir(vtkdir, 0777);  
  
  bool UpdateConvection=FALSE, UpdateRhs=FALSE, Initialize_ScalarSolver=FALSE;
  bool ConvectionFirstTime=TRUE;

  bool UpdatePBEConvection=FALSE, UpdatePBERhs=FALSE, Initialize_PBEScalarSolver=TRUE;
  bool PBEConvectionFirstTime=TRUE;
  bool *DirichletBDPt,  FACTORIZE=TRUE;
  //======================================================================
  // read parameter file
  //======================================================================
  total_time = GetTime();
  if(argc>=2)
   { ret=Domain->ReadParam(argv[1]); }
  else
   { ret=Domain->ReadParam(Readin); }

      OpenFiles();
  OutFile.setf(std::ios::scientific);

  if(TDatabase::ParamDB->DISCTYPE != LOCAL_PROJECTION)
  {
    // switch off all local projection terms
    TDatabase::ParamDB->LP_FULL_GRADIENT = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT = 1;

    TDatabase::ParamDB->LP_STREAMLINE = 0;
    TDatabase::ParamDB->LP_STREAMLINE_COEFF = 0;
    TDatabase::ParamDB->LP_STREAMLINE_EXPONENT = 1;
  }

  if(TDatabase::ParamDB->LP_FULL_GRADIENT)
  {
    if(TDatabase::ParamDB->LP_STREAMLINE)
    {
      TDatabase::ParamDB->LP_STREAMLINE = 0;
      TDatabase::ParamDB->LP_STREAMLINE_COEFF = 0;
      TDatabase::ParamDB->LP_STREAMLINE_EXPONENT = 1;
      OutPut("local projection stabilisation in streamline direction ");
      OutPut("is switched off due to stabilisation of full gradient." << endl);
    }
  }

  if(TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE == -123)
    TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;

  if(TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE == -123)
    TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;

  // write parameters into outfile
  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();

  GetExampleFileData(BoundaryConditions, BoundValues, InitiaValues, Coefficients,
                     N_PBEqns, N_IndepntScalarEqns, Disctypes);
  N_ScalarEqns=N_PBEqns+N_IndepntScalarEqns;
 
  
  if(TDatabase::ParamDB->DISCTYPE == GALERKIN)
    if (TDatabase::ParamDB->SOLD_TYPE)
  {
    TDatabase::ParamDB->SOLD_TYPE = 0;
    OutPut("SOLDTYPE set to 0 !!!" << endl);
  }

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
  for(i=0;i<TDatabase::ParamDB->P9;i++)
    Domain->RefineallxDirection();
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
  
// exit(0);

  // ORDER should be same for all scalar and Physical PBE spaces 
  // due to coupling, but PBE X- dir can have multiple disc
  PBE_ORDER = TDatabase::ParamDB->PBE_P8;
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  VELOCITYORDER = TDatabase::ParamDB->VELOCITY_SPACE;

  // initializ time
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = TDatabase::TimeDB->TIMESTEPLENGTH;
  SetTimeDiscParameters();
  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR;

  t3 = GetTime();
  total_time = t3 - total_time;
  SetPolynomialDegree();

  // check the example file, to activate
  VeloFunction=TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD;

  // +1 space for multiple disc fespace for PBE in X-direction
  Scalar_Spaces = new TFESpace2D*[N_ScalarEqns + 1];
  Scalar_Spaces_low = new TFESpace2D*[N_ScalarEqns + 1];
  sqstructureA = new TSquareStructure2D*[N_ScalarEqns];

  N_Uarray = new int[N_ScalarEqns + 1];
  N_Active = new int[N_ScalarEqns + 1];
  N_Uarray_low = new int[N_ScalarEqns + 1];
  N_Active_low = new int[N_ScalarEqns + 1];
    
  MatricesA = new TSquareMatrix2D*[N_ScalarEqns];
  MatricesM = new TSquareMatrix2D*[N_ScalarEqns];
  MatricesK = new TSquareMatrix2D*[N_ScalarEqns];
  SolArray = new double*[N_ScalarEqns + 1];
  SolArray_low = new double*[N_ScalarEqns + 1];
  RhsArray = new double*[N_ScalarEqns];

  B = new double*[N_ScalarEqns];

  ScalarFunctions = new TFEFunction2D*[N_ScalarEqns + 1];
  ScalarFunctions_low = new TFEFunction2D*[N_ScalarEqns + 1];

  
   MatValues_Solver = new double *[N_ScalarEqns];
   KCol_Solver = new int *[N_ScalarEqns];
   Row_Solver = new int *[N_ScalarEqns];
   Symbolic_Solver = new void *[N_ScalarEqns];
   Numeric_Solver = new void *[N_ScalarEqns];    
  
  
#ifdef _OMP
  if(TDatabase::ParamDB->SOLVER_TYPE==100)
    Scalar_SMPSolver = new TParDirectSolver *[N_ScalarEqns];
#endif

  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  OutPut("N_Cells (space) : " << N_Cells <<endl);

#ifdef __SIMPATURS__
  BoundCondition_LminLMax(cond_Lmin, cond_Lmax);
#endif
  
// #pragma omp parallel for
// for(i=0;i<100;i++)
// {
//    int tid = omp_get_thread_num();
//     cout << i << " test " << tid << endl;   
// }
// 
// exit(0);

  //=========================================================================
  // construct all finite element spaces
  //=========================================================================
  for(i=0;i<N_IndepntScalarEqns;i++)
  {
    // fespaces for scalar equations
    Scalar_Spaces[i] =  new TFESpace2D(coll, Name, NameStrings[i], BoundaryConditions[i], ORDER, NULL);  
    
    
    N_Uarray[i] = Scalar_Spaces[i]->GetN_DegreesOfFreedom();
    N_Active[i] = Scalar_Spaces[i]->GetActiveBound();

    if(Max_N_Unknowns<N_Uarray[i])  Max_N_Unknowns=N_Uarray[i];

    OutPut("DOF Scalar : " << i << setw(10) << N_Uarray[i] << endl);      
    
    //=========================================================================
    // memory allocate all vectors and construction of all fefunction
    //=========================================================================
    SolArray[i] =  new double[N_Uarray[i]];
    B[i] = new double [N_Uarray[i]];

    memset(SolArray[i], 0, N_Uarray[i]*SizeOfDouble);
    memset(B[i], 0, N_Uarray[i]*SizeOfDouble);

    ScalarFunctions[i] = new TFEFunction2D(Scalar_Spaces[i], NameStrings[i], NameStrings[i], SolArray[i], N_Uarray[i]);

    ScalarFunctions[i]->Interpolate(InitiaValues[i]);

    RhsArray[i] = new double[N_Uarray[i]];
    memset(RhsArray[i], 0, N_Uarray[i]*SizeOfDouble);

    
    
     if(PBE_ORDER != ORDER)
     {
      Scalar_Spaces_low[i] =  new TFESpace2D(coll, Name, NameStrings[i], BoundaryConditions[i], PBE_ORDER, NULL); 
      N_Uarray_low[i] = Scalar_Spaces_low[i]->GetN_DegreesOfFreedom();
      N_Active_low[i] = Scalar_Spaces_low[i]->GetActiveBound();       
       
      SolArray_low[i] =  new double[N_Uarray_low[i]];
      memset(SolArray_low[i], 0, N_Uarray_low[i]*SizeOfDouble);
      
      ScalarFunctions_low[i] = new TFEFunction2D(Scalar_Spaces_low[i], NameStrings[i], NameStrings[i], SolArray_low[i], N_Uarray_low[i]);
     }
    else
     {
      Scalar_Spaces_low[i] = Scalar_Spaces[i];
      N_Uarray_low[i] = N_Uarray[i];
      N_Active_low[i] = N_Active[i];       
      SolArray_low[i] = SolArray[i];
      ScalarFunctions_low[i] =  ScalarFunctions[i];
     }  
    
    
    //=========================================================================
    // allocate memory for all matrices
    //=========================================================================
    // build matrices
    // first build matrix structure
    sqstructureA[i] = new TSquareStructure2D(Scalar_Spaces[i]);
    sqstructureA[i]->Sort();   // sort column numbers: numbers are increasing

      // two matrices used
      // M is the mass matrix, also a system matrix
      // the iterative solver uses M
      sqmatrixM = new TSquareMatrix2D(sqstructureA[i]);
      MatricesM[i] = sqmatrixM;

      // A contains the non time dependent part of the discretization
      sqmatrixA = new TSquareMatrix2D(sqstructureA[i]);
      MatricesA[i] = sqmatrixA;

      if(Disctypes[i] == SDFEM)
      {
        // stabilisation matrix K
        sqmatrixK = new TSquareMatrix2D(sqstructureA[i]);
        MatricesK[i] = sqmatrixK;

        // stabilisation matrix S, inconsistent in time, so no need S
//         sqmatrixK = new TSquareMatrix2D(sqstructureA[i]);
//         MatricesK[i] = sqmatrixK;
      }
  }       // for(i=0;i<N_IndepntScalarEqns;i++)

  //=========================================================================
  // PBS setting begin
  // construct FESpace and matrices for population balance equation
  //=========================================================================
#ifdef __PBS__
  i = N_IndepntScalarEqns;
  PBE_INDEX = N_IndepntScalarEqns;
  
  if(PBE_ORDER != ORDER)
   { PBE_INDEX_LOW = PBE_INDEX + 1; }
  else
   { PBE_INDEX_LOW = PBE_INDEX; }
  
  N_PBEqns = 1;
  N_ScalarEqns = N_IndepntScalarEqns + N_PBEqns;
  // no. of population balance equations (PBE)
  //=========================================================================
  // fespaces for population balance equation in physical space
  // FeSpace and memory for all matrices
  // assume that the convection and the reaction coefficient terms are independent
  // of internal coordinates, so lhs matrices are same for all lelvels of internal coordinate
  //=========================================================================
  Scalar_Spaces[PBE_INDEX] = new TFESpace2D(coll, Name, NameStrings[PBE_INDEX],
    BoundaryConditions[PBE_INDEX], ORDER, NULL);
  
  N_Uarray[PBE_INDEX] = Scalar_Spaces[PBE_INDEX]->GetN_DegreesOfFreedom();
  N_Active[PBE_INDEX] = Scalar_Spaces[PBE_INDEX]->GetActiveBound();
  
  if(Max_N_Unknowns<N_Uarray[PBE_INDEX])  Max_N_Unknowns=N_Uarray[PBE_INDEX];
  
  OutPut("DOF of PBE physical space : "<< setw(10) << N_Uarray[PBE_INDEX] << endl);

  N_OutEdges=0;
  //collect information of outlet boundary  
  for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);
      N_Edges = cell->GetN_Edges();

      for(j=0;j<N_Edges;j++)
       { 
        Joint = cell->GetJoint(j);
        if(Joint->GetType() == BoundaryEdge)
         {
          Outlet_Joint = (TBoundEdge *)Joint;
          BoundComp = Outlet_Joint->GetBoundComp();
          comp=BoundComp->GetID();
	  if(comp==1)
	   N_OutEdges ++;
         }     
       }// for(j=0;j    
     } // for(i=0;i<N_Ce

    OutPut(" N_Outlet Edges " << N_OutEdges << endl);
   
    Out_CellIndex = new int[N_OutEdges];
    Out_EdgeIndex = new int[N_OutEdges];
    N_OutEdges=0;
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);
      N_Edges = cell->GetN_Edges();

      for(j=0;j<N_Edges;j++)
       { 
        Joint = cell->GetJoint(j);
        if(Joint->GetType() == BoundaryEdge)
         {
          Outlet_Joint = (TBoundEdge *)Joint;
          BoundComp = Outlet_Joint->GetBoundComp();
          comp=BoundComp->GetID();
	  if(comp==1)
	  {
	   Out_CellIndex[N_OutEdges] = i;
	   Out_EdgeIndex[N_OutEdges] = j;
// 	       OutPut(" N_Outlet Edges " << N_OutEdges << " " << i << " " << j  <<  endl);
	   N_OutEdges ++;
	  }
         }     
       }// for(j=0;j    
     } // for(i=0;i<N_Ce   
   
  //=========================================================================
  // memory allocate all vectors and construction for PBS space
  //=========================================================================
  // build matrices
  // first build matrix structure
  sqstructureA[PBE_INDEX] = new TSquareStructure2D(Scalar_Spaces[PBE_INDEX]);
  sqstructureA[PBE_INDEX]->Sort();  // sort column numbers: numbers are increasing


    // two matrices used
    // M is the mass matrix, also a system matrix
    // the iterative solver uses M
    sqmatrixM = new TSquareMatrix2D(sqstructureA[PBE_INDEX]);
    MatricesM[PBE_INDEX] = sqmatrixM;

    // A contains the non time dependent part of the discretization
    sqmatrixA = new TSquareMatrix2D(sqstructureA[PBE_INDEX]);
    MatricesA[PBE_INDEX] = sqmatrixA;

    if(Disctypes[PBE_INDEX] == SDFEM)
     {
      // stabilisation matrix K
      sqmatrixK = new TSquareMatrix2D(sqstructureA[PBE_INDEX]);
      MatricesK[PBE_INDEX] = sqmatrixK;
     }    
    
    //=========================================================================
    // Construct FeSpace and all data for the internal domain
    // FESpace is same in internal coordinate for all  QuadPts
    //=========================================================================
   /** Pbe X - fe space for L-direction is disc space,
       but x-direction will be solved with ORDER fespace*/
  if(PBE_ORDER != ORDER)
   {
    // low order space
    Scalar_Spaces[PBE_INDEX_LOW] = new TFESpace2D(coll, Name, NameStrings[PBE_INDEX],
                                  BoundaryConditions[PBE_INDEX], PBE_ORDER, NULL); 
  
    N_Uarray[PBE_INDEX_LOW] = Scalar_Spaces[PBE_INDEX_LOW]->GetN_DegreesOfFreedom();
    N_Active[PBE_INDEX_LOW] = Scalar_Spaces[PBE_INDEX_LOW]->GetActiveBound();      
   }
      
    
    N = int(TDatabase::ParamDB->REACTOR_P11);
    L0 = TDatabase::ParamDB->REACTOR_P12;
    L1 = TDatabase::ParamDB->REACTOR_P13;

    Generate1DMesh(Domain_Intl, L0, L1, N);
    Coll_Intl = Domain_Intl->GetCollection(It_Finest, 0);
    N_Cells_Intl= Coll_Intl->GetN_Cells();

    cout<< "N_Cells_Internal " << N_Cells_Intl <<endl;


//     FE1D_List_Intl = new FE1D[N_Cells_Intl];
//     for(j=0;j<N_Cells_Intl;j++)
//       FE1D_List_Intl[j] = (FE1D)TDatabase::ParamDB->TEST_ORDER;

//     FeSpace_Intl = new TFESpace1D(Coll_Intl, VString, VString, FE1D_List_Intl);
//     if(TDatabase::ParamDB->TEST_ORDER<-10)

    if(TDatabase::ParamDB->INTL_DISCTYPE==FD)
     { 
      // Finite difference in internal 
      TDatabase::ParamDB->TEST_ORDER = 1;
     }
    
    
    FeSpace_Intl = new TFESpace1D(Coll_Intl, VString, VString, TDatabase::ParamDB->TEST_ORDER);

    N_LDof = FeSpace_Intl->GetN_DegreesOfFreedom();

    if(TDatabase::ParamDB->TEST_ORDER<-9)
     {
      FeSpace_Intl->SetAsDGSpace();
      dGDisc = 1;
      TDatabase::ParamDB->INTL_DISCTYPE=DG;  // no supg method
     }

    N_PBE_IntlLevels =  FeSpace_Intl->GetN_DegreesOfFreedom();
    GetInternalNodalPts(FeSpace_Intl, N_PBE_IntlLevels, IntlPosL);

    OutPut("Dof PBE Internal space : "<< setw(14) << N_PBE_IntlLevels << endl);

    SqStruct1D_Intl = new TSquareStructure1D(FeSpace_Intl);
    SqStruct1D_Intl->Sort();
    M_Intl = new TSquareMatrix1D(SqStruct1D_Intl);
    A_Intl = new TSquareMatrix1D(SqStruct1D_Intl);

    if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
    {
      S_Intl = new TSquareMatrix1D(SqStruct1D_Intl);
      K_Intl = new TSquareMatrix1D(SqStruct1D_Intl);
    }

    if(TDatabase::ParamDB->REACTOR_P5==0)
      { 
       GetNodalPtsADI_System(N_PBE_SpatialPts, FeSpace_Intl, Scalar_Spaces[PBE_INDEX_LOW], MaxN_PtsForNodal, IntlX, IntlY); 
      }
    else
       { 
        GetQuadPtsADI_System(N_PBE_SpatialPts, Scalar_Spaces[PBE_INDEX_LOW], IntlX, IntlY); 
       }

        OutPut("Dof PBE Internal space N_PBE_SpatialPts : "<<  N_PBE_SpatialPts << endl);
        ADI_System = new TADISystem1D*[N_PBE_SpatialPts];

       N_f_OutletIndex =0;
       for(i=0; i<N_PBE_SpatialPts; i++)
        {
// 	 if(fabs(IntlX[i]-210)<1.e-12  && fabs(IntlY[i]-0.25)<1.e-12)
// 	  {
//            OutPut(i << " IntlX[i] " <<IntlX[i]<< " IntlY[i] " <<IntlY[i] <<endl);
// 	   f_Integ_NodalPt = i;
// 	   break;
// 	  }
	 if(fabs(IntlX[i] -200.)<1.e-12  && fabs(IntlY[i]-0.0)<1.e-12)
	  {
           OutPut(i << " IntlX[i] " <<IntlX[i]<< " IntlY[i] " <<IntlY[i] <<endl);
	   f_Integ_NodalPt[0] = i;
// 	   break;
	  }

	 if(fabs(IntlX[i] -200.)<1.e-12  && fabs(IntlY[i]-0.125)<1.e-12)
	  {
           OutPut(i << " IntlX[i] " <<IntlX[i]<< " IntlY[i] " <<IntlY[i] <<endl);
	   f_Integ_NodalPt[1] = i;
// 	   break;
	  }

	 if(fabs(IntlX[i] -200.)<1.e-12  && fabs(IntlY[i]-0.25)<1.e-12)
	  {
           OutPut(i << " IntlX[i] " <<IntlX[i]<< " IntlY[i] " <<IntlY[i] <<endl);
	   f_Integ_NodalPt[2] = i;
// 	   break;
	  }

	 if(fabs(IntlX[i] -200.)<1.e-12  && fabs(IntlY[i]-0.375)<1.e-12)
	  {
           OutPut(i << " IntlX[i] " <<IntlX[i]<< " IntlY[i] " <<IntlY[i] <<endl);
	   f_Integ_NodalPt[3] = i;
	  }
  
	 if(fabs(IntlX[i] -200.)<1.e-12  && fabs(IntlY[i]-0.5)<1.e-12)
	  {
           OutPut(i << " IntlX[i] " <<IntlX[i]<< " IntlY[i] " <<IntlY[i] <<endl);
	   f_Integ_NodalPt[4] = i;
	  }

	 if(fabs(IntlX[i] -0.)<1.e-12  && fabs(IntlY[i]-0.5)<1.e-12)
	  {
           OutPut(i << " IntlX[i] " <<IntlX[i]<< " IntlY[i] " <<IntlY[i] <<endl);
           f_Inlet_SpaceIndex = i;
	  }

         // for calculating avgQ_3 at outlet
	 if(fabs(IntlX[i] - 200.)<1.e-12)
	  {
           OutPut(i << " IntlX[i] " <<IntlX[i]<< " IntlY[i] " <<IntlY[i] <<endl);
	   N_f_OutletIndex++;
	  }


	}
//    cout << "N_f_OutletIndex " <<N_f_OutletIndex << endl;

      if(N_f_OutletIndex)
       { f_OutletIndex = new int[N_f_OutletIndex]; }
      N_f_OutletIndex = 0;
       for(i=0; i<N_PBE_SpatialPts; i++)
        {
  	 if(fabs(IntlX[i] - 200.)<1.e-12)
	  {
//            OutPut(i << " IntlX[i] " <<IntlX[i]<< " IntlY[i] " <<IntlY[i] <<endl);
	   f_OutletIndex[N_f_OutletIndex] = i;
	   N_f_OutletIndex++;
	  }
	}

//        exit(0);

    Sol_IntlLoc = new double[N_PBE_IntlLevels];
    OldSol_IntlLoc = new double[N_PBE_IntlLevels];
    B_IntlLoc = new double[N_PBE_IntlLevels];
    defect_IntlLoc = new double[N_PBE_IntlLevels];


    memset(Sol_IntlLoc, 0, N_PBE_IntlLevels*SizeOfDouble);
    memset(OldSol_IntlLoc, 0, N_PBE_IntlLevels*SizeOfDouble);
    memset(B_IntlLoc, 0, N_PBE_IntlLevels*SizeOfDouble);
    memset(defect_IntlLoc, 0, N_PBE_IntlLevels*SizeOfDouble);

    FeFunction_Intl = new TFEFunction1D(FeSpace_Intl, NameStrings[PBE_INDEX], NameStrings[PBE_INDEX], 
                                        Sol_IntlLoc,  N_PBE_IntlLevels);

    for(i=0; i<N_PBE_SpatialPts; i++)
     {
      ADI_System[i] = new TADISystem1D(FeFunction_Intl, M_Intl, A_Intl, S_Intl, K_Intl, IntlX[i], IntlY[i],
                                       InitialCondition_Psd_Intl,  BoundValue_LMin, BoundValue_LMax,
                                       Sol_IntlLoc, OldSol_IntlLoc, B_IntlLoc, defect_IntlLoc, IntlPosL,
#ifdef __PBSConstT__   
                                       NULL
#else
                                       Get_GrowthAndB_Nuc
#endif
                                       );
     }
     
    C_IntlPtValues = new double[N_PBE_SpatialPts];
    T_IntlPtValues = new double[N_PBE_SpatialPts];
    C_Sat_IntlPtValues = new double[N_PBE_SpatialPts];
    
// #ifdef __AXIAL3D__ 
//     Velo_PBE_SpatialPts = new double[3*N_PBE_SpatialPts];    
//     GradVelo_PBE_SpatialPts = new double[9*N_PBE_SpatialPts];
// #else
    Velo_PBE_SpatialPts = new double[2*N_PBE_SpatialPts];    
    GradVelo_PBE_SpatialPts = new double[4*N_PBE_SpatialPts];
// #endif    
    
    
    i = N_PBE_SpatialPts*N_PBE_IntlLevels;

    PBE_IntlPtValuesT = new double[i];
    OldPBE_IntlPtValuesT = new double[i];
    PBE_IntlPtValues = new double[i];
    PBE_IntlPtRhsValues = new double[i];
    DirichletBDPt = new bool[N_PBE_SpatialPts];

    i= N_LDof*N_PBE_SpatialPts;
    PBE_LSolT =  new double[i];
    PBE_IntlPtRhsValuesT = new double[i];
    OldPBE_IntlPtRhsValuesT = new double[i];
    Oldrhs_IntlLoc = new double[N_LDof];

    memset(OldPBE_IntlPtRhsValuesT, 0, i*SizeOfDouble);
    
//     L_Cube = new double[N_PBE_IntlLevels];
// 
//     for(i=0; i<N_PBE_IntlLevels; i++)
//       L_Cube[i] = pow(IntlPosL[i], 3.0);

//        for(i=0; i<N_PBE_IntlLevels; i++)
//         cout <<i<< " IntlPosL[i] " <<IntlPosL[i]  << endl;
// //         cout <<i<< " IntlPosL[i] " <<IntlPosL[i] << " "<< L_Cube[i] << endl;

//        exit(0);

  Out_Level = N_PBE_IntlLevels/4;

  if(cond_Lmin==DIRICHLET && !dGDisc)
   { start_pbe = 1; }
  else
   { start_pbe = 0; }

  if(cond_Lmax==DIRICHLET && !dGDisc) { end_pbe = N_PBE_IntlLevels-1; }
  else { end_pbe = N_PBE_IntlLevels; }

//   MPI_Finalize();
//   exit(0);
  //========================================================================================
  // memory allocate for all vectors and construction of PBE fefunction in
  // physical space (for all internal levels)
  //========================================================================================
  N_U = N_Uarray[PBE_INDEX];
  N_U_low = N_Uarray[PBE_INDEX_LOW] ;
  SolPbe = new double[N_PBE_IntlLevels*N_U];
  B_Pbe = new double[N_PBE_IntlLevels*N_U];
  
  if(N_U>N_U_low) N_U_Max = N_U;
  else N_U_Max = N_U_low;
  
  PbeAux =  new double[N_U_Max];
  
  if(PBE_ORDER != ORDER)
   {
    SolPbe_low = new double[N_PBE_IntlLevels*N_U_low]; 
    memset(SolPbe_low, 0, N_PBE_IntlLevels*N_U_low*SizeOfDouble);    
   }
  else
   {
    SolPbe_low = SolPbe; 
    memset(SolPbe, 0, N_PBE_IntlLevels*N_U*SizeOfDouble);    
   }
         
  RhsArray_Pbe = new double[N_PBE_IntlLevels*N_U];
  OldRhsArray_Pbe = new double[N_PBE_IntlLevels*N_U];
  memset(RhsArray_Pbe, 0, N_PBE_IntlLevels*N_U*SizeOfDouble);
  memset(OldRhsArray_Pbe, 0, N_PBE_IntlLevels*N_U*SizeOfDouble);
 

  SolPbe_Output = new double[N_U];
  memset(SolPbe_Output, 0, N_U*SizeOfDouble);
  PbeFunctions = new  TFEFunction2D(Scalar_Spaces[PBE_INDEX], NameStrings[PBE_INDEX], NameStrings[PBE_INDEX], SolPbe_Output, N_U);

  // prepare output, only the concentration will be saved
  Output_PBS = new TOutput2D(0,  1, 1, 1, Domain);
  Output_PBS->AddFEFunction(PbeFunctions);

  ScalarPbeFunctions = new TFEFunction2D* [N_PBE_IntlLevels];
   
  Pbe_xdir = new TFEVectFunct2D(Scalar_Spaces[PBE_INDEX], NameStrings[PBE_INDEX], NameStrings[PBE_INDEX], SolPbe, N_U, N_PBE_IntlLevels);
  
  for(j=0; j<N_PBE_IntlLevels; j++)
   ScalarPbeFunctions[j] = Pbe_xdir->GetComponent(j);
     
  if(PBE_ORDER != ORDER)
   {
    ScalarPbeFunctions_low = new TFEFunction2D* [N_PBE_IntlLevels];
    Pbe_xdir_low = new TFEVectFunct2D(Scalar_Spaces[PBE_INDEX_LOW], NameStrings[PBE_INDEX], NameStrings[PBE_INDEX], SolPbe_low, N_U_low, N_PBE_IntlLevels);
    
   
    for(j=0; j<N_PBE_IntlLevels; j++)
     ScalarPbeFunctions_low[j] = Pbe_xdir_low->GetComponent(j);   
   }
  else
   {
     ScalarPbeFunctions_low = ScalarPbeFunctions;
   }

   // calculate q3 at every dof, and then calculate the mass
   pbeq3massAtNodalPts = new double[N_PBE_SpatialPts];
   pbeq3mass = new double[N_U_low];
   PSD_mass =  new TFEFunction2D(Scalar_Spaces[PBE_INDEX_LOW], NameStrings[PBE_INDEX], NameStrings[PBE_INDEX],  pbeq3mass, N_U_low);
   
    
//     ScalarPbeFunctions[j] = new  TFEFunction2D(Scalar_Spaces[PBE_INDEX], NameStrings[PBE_INDEX], NameStrings[PBE_INDEX], SolPbe+j*N_U, N_U);

  // Identify the dirichlet internal points -begin
    SolPbe_Lmin = new double[N_U];
    SolPbe_Lmax = new double[N_U];

    if(TDatabase::ParamDB->REACTOR_P5==0 && PBE_ORDER == ORDER)
     {
      // set 1 for all DOF at L_Min
      for(i=0;i<N_U;i++)
        SolPbe[i] = 1.;

      // set 1 for Dirichlet DOF and else 0 at L_Max
      SolPbe_OneLevel = SolPbe+(N_PBE_IntlLevels-1)*N_U;
      memset(SolPbe_OneLevel, 0, N_U*SizeOfDouble);
      for(i=N_Active[PBE_INDEX];i<N_U;i++)
        SolPbe_OneLevel[i] = 1.;

      /** transfer sol from spatial to internal */
      GetPtsValuesForIntl(N_PBE_IntlLevels, ScalarPbeFunctions_low, SolPbe_low, N_U_low, PBE_IntlPtValuesT, N_PBE_SpatialPts);

      for(i=0; i<N_PBE_SpatialPts; i++)
      {
        if( fabs(PBE_IntlPtValuesT[i*N_PBE_IntlLevels] - PBE_IntlPtValuesT[i*N_PBE_IntlLevels + N_PBE_IntlLevels-1])<1.e-8 )
        {
          //cout <<i<< " IntlX[i] " <<IntlX[i]<< " IntlY[i] " <<IntlY[i] <<endl;
          // no computation need for these points in internal direction,
          //i.e., for nodal points on inlet, outlet, walls
          DirichletBDPt[i] = TRUE;
        }

      }
     }
    else
     {
      // assumed that no quad points of X-direction lie on the boundary
      for(i=0; i<N_PBE_SpatialPts; i++)
      DirichletBDPt[i] = FALSE;
     }

  memset(SolPbe, 0, N_PBE_IntlLevels*N_U*SizeOfDouble);      
//   memset(SolPbe, 0, N_PBE_IntlLevels*N_U*SizeOfDouble);
  //========================================================================================
  // Identify the dirichlet internal points - end
  // PBS setting end
  //========================================================================================
  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);
#endif
  
  //  each scalar equation including PBE in space will use same velo function
  if(VeloFunction)
  {
    GetVelocityAndPressureSpace(coll, BoundCondition_NSE,
      mortarcoll, Velocity_Spaces,
      pressure_space, &pressure_space_code,
      TDatabase::ParamDB->VELOCITY_SPACE,
      TDatabase::ParamDB->PRESSURE_SPACE);

     N_Veloarray = Velocity_Spaces->GetN_DegreesOfFreedom();
     OutPut("DOF Velocity : "<< setw(10) << N_Veloarray << endl);


//       strcat(ReadGrapeBaseName, "0.Sol");
      
      UArray = new double[2*N_Veloarray];

      memset(UArray, 0, 2*N_Veloarray*SizeOfDouble);

      // vector fe function
      velocity = new TFEVectFunct2D(Velocity_Spaces, UString, UString, UArray, N_Veloarray, 2);
      // individual components of velocity
      velo1 = velocity->GetComponent(0);
      velo2 = velocity->GetComponent(1);


//       velocity->ReadSol(ReadGrapeBaseName);
                velo1->Interpolate(InitialU1);
//                 velo2->Interpolate(InitialU2);

  } //if(VeloFunction)
  

  
// #ifdef __MOVINGMESH__    
//    TDatabase::TimeDB->CURRENTTIME=0.;
//    Grid_space = new TFESpace2D(coll, VString, VString, GridBoundCondition, 1, NULL); 
//    N_GridDOFs = Grid_space->GetN_DegreesOfFreedom(); 
//    
//    GlobalNumbers_Grid = Grid_space->GetGlobalNumbers();
//    BeginIndex_Grid    = Grid_space->GetBeginIndex();   
//    
//    WArray = new double[2*N_GridDOFs];
//    memset(WArray, 0, 2*N_GridDOFs*SizeOfDouble);
// 
//    // vector fe function
//    GridVelo = new TFEVectFunct2D(Grid_space, WString, WString, WArray, N_GridDOFs, 2);
//    // individual components of velocity
//    gridvelo1 = GridVelo->GetComponent(0);
//    gridvelo2 = GridVelo->GetComponent(1); 
//    
//    gridpos = new double[2*N_GridDOFs];
//    gridpos_old = new double[2*N_GridDOFs];   
//    memset(gridpos, 0, 2*N_GridDOFs*SizeOfDouble);
//    GridPos = new TFEVectFunct2D(Grid_space, WString, WString, gridpos, N_GridDOFs, 2);       
//    GridPos->GridToData();
//    memcpy(gridpos_old, gridpos, 2*N_GridDOFs*SizeOfDouble);
// #endif  //__MOVINGMESH__ 
//    
  if(N_IndepntScalarEqns)
  {
    // prepare output, only the concentration will be saved
    Output = new TOutput2D(2, N_IndepntScalarEqns, 1, 1, Domain);

    if(VeloFunction)
      Output->AddFEVectFunct(velocity);
    
// #ifdef __MOVINGMESH__  
//       Output->AddFEVectFunct(GridPos);
// #endif        

    for(i=0;i<N_IndepntScalarEqns;i++)
      Output->AddFEFunction(ScalarFunctions[i]);

    os.seekp(std::ios::beg);
    Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());
  }

  //======================================================================
  // assembling of mass matrix and initial rhs (f_0) for IndepntScalarEqns
  //======================================================================
    // Assemble matrices
    for(i=0;i<N_IndepntScalarEqns;i++)
    {
      // set parameters
      N_Rhs = 1;
      N_FESpaces = 1;
      fesp[0] = Scalar_Spaces[i];

      if (VeloFunction)
      {
        N_FESpaces = 2;
        fesp[1] = Velocity_Spaces;
        fefct[0] = velo1;
        fefct[1] = velo2;

        // defined in TCD2D.h
        aux =  new TAuxParam2D(TimeCDParamsVeloFieldN_FESpaces,
          TimeCDParamsVeloFieldN_Fct,
          TimeCDParamsVeloFieldN_ParamFct,
          TimeCDParamsVeloFieldN_FEValues,
          fesp+1, fefct,
          TimeCDParamsVeloFieldFct,
          TimeCDParamsVeloFieldFEFctIndex,
          TimeCDParamsVeloFieldFEMultiIndex,
          TimeCDParamsVeloFieldN_Params,
          TimeCDParamsVeloFieldBeginParam);
      }
      else
        { aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 1, NULL); }

        // reset matrices
      N_SquareMatrices = 1;
      SQMATRICES[0] = MatricesM[i];
      SQMATRICES[0]->Reset();    
      DiscreteForm = DiscreteFormMatrixMRhs[i];
      
//       if(Disctypes[i] == SDFEM)
//         {
//          N_SquareMatrices = 2;
//          SQMATRICES[1] = MatricesS[i];
//          SQMATRICES[1]->Reset();
//          DiscreteForm = DiscreteFormMatrixMRhs_SUPG[i];
//         }
 
      BDCond[0] = BoundaryConditions[i];
      BDValue[0] = BoundValues[i];

      memset(RhsArray[i], 0, N_Uarray[i]*SizeOfDouble);
      RHSs[0] = RhsArray[i];
      ferhs[0] = Scalar_Spaces[i];

      Assemble2D(N_FESpaces, fesp,
        N_SquareMatrices, SQMATRICES,
        0, NULL,
        N_Rhs, RHSs, ferhs,
        DiscreteForm,
        BDCond,
        BDValue,
        aux);
      delete aux;

      // copy Dirichlet values from rhs into sol
      memcpy(SolArray[i]+N_Active[i], RHSs[0]+N_Active[i], (N_Uarray[i]-N_Active[i])*SizeOfDouble);
    }  // for(i=0;i<N_IndepntScalarEqns;i++)


#ifdef __PBS__
    // interpolate initial solution at L direction nodal-interpolation point levels
    TDatabase::TimeDB->CURRENTTIME=0.;
    memset(PBE_IntlPtValuesT, 0, N_PBE_SpatialPts*N_PBE_IntlLevels*SizeOfDouble);
    for(i=0;i<N_PBE_SpatialPts;i++)
    {
      BDUpdatedSolPbe_Intl_Loc = PBE_IntlPtValuesT + i*N_PBE_IntlLevels;
      //        ADI_System[i]->Interpolate(BDUpdatedSolPbe_Intl_Loc, Exact_Psd_Intl);
      ADI_System[i]->Interpolate(BDUpdatedSolPbe_Intl_Loc, InitialCondition_Psd_Intl);
    }
 

    /** transfer sol from internal to spatial */
    //  ScalarPbeFunctions_low= ScalarPbeFunctions if PBE_ORDER == ORDER
    if(TDatabase::ParamDB->REACTOR_P5==0)
     {
      GetSolFromNodalPtVales(N_PBE_IntlLevels, MaxN_PtsForNodal, ScalarPbeFunctions_low, SolPbe_low, N_U_low, PBE_IntlPtValuesT, N_PBE_SpatialPts);
     }
    else
     {
      GetSolFromQuadPtVales(N_PBE_IntlLevels, MaxN_PtsForNodal, ScalarPbeFunctions_low, SolPbe_low, N_U_low, PBE_IntlPtValuesT, N_PBE_SpatialPts);
     }  
     
     
   /** Prolongate sol from low fespace to high fespace */
   if(PBE_ORDER != ORDER)
    Prolongate(Scalar_Spaces[PBE_INDEX_LOW], Scalar_Spaces[PBE_INDEX], N_PBE_IntlLevels,
               SolPbe_low, SolPbe, PbeAux);
    
//     //set the boundary value at L_min and L_max
       // not necessar here, as DIRICHLET BC at L_min is not solved spatially
       // it is enough to set the L_min 
//     if(!dGDisc)
//     {
//      
//       
//       
//      SetLMinLMaxBoundValue(N_PBE_IntlLevels, ScalarPbeFunctions, N_U, RhsArray_Pbe);
//     }
    
    RhsArray[PBE_INDEX] = new double[N_U];
    memset(RhsArray[PBE_INDEX], 0, N_U*SizeOfDouble);
  
    for(j=start_pbe; j<end_pbe; j++)
    {
      TDatabase::ParamDB->REACTOR_P29=IntlPosL[j];

      // set parameters
      N_Rhs = 1;
      N_FESpaces = 1;
      fesp[0] = Scalar_Spaces[PBE_INDEX];
      ferhs[0] = Scalar_Spaces[PBE_INDEX];

      // reset matrices
      N_SquareMatrices = 1;
      SQMATRICES[0] = MatricesM[PBE_INDEX];
      SQMATRICES[0]->Reset();
      DiscreteForm = DiscreteFormMatrixMRhs[PBE_INDEX];
      
/*       if(Disctypes[PBE_INDEX] == SDFEM)
        {
         N_SquareMatrices = 2;
         SQMATRICES[1] = MatricesS[PBE_INDEX];
         SQMATRICES[1]->Reset();
         DiscreteForm = DiscreteFormMatrixMRhs_SUPG[PBE_INDEX];
        }     
  */    
      
      BDCond[0] = BoundaryConditions[PBE_INDEX];
      BDValue[0] = BoundValues[PBE_INDEX];

      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 1, NULL);

      memset(RhsArray_Pbe+j*N_U, 0, N_U*SizeOfDouble);
      RHSs[0] = RhsArray_Pbe+j*N_U;
    
      Assemble2D(N_FESpaces, fesp,
        N_SquareMatrices, SQMATRICES,
        0, NULL,
        N_Rhs, RHSs, ferhs,
        DiscreteForm,
        BDCond,
        BDValue,
        aux);
      delete aux;
 
      // copy Dirichlet values to the solution
      SolPbe_OneLevel = SolPbe+j*N_U;
      memcpy(SolPbe_OneLevel+N_Active[PBE_INDEX], RhsArray_Pbe+(j*N_U+N_Active[PBE_INDEX]),
             (N_U-N_Active[PBE_INDEX])*SizeOfDouble);
    }  //  for(j=start_pbe; j<end_pbe; j++)
    
 
    //    memcpy(OldSolPbe, SolPbe, N_U*N_PBE_IntlLevels*SizeOfDouble);
    memcpy(OldRhsArray_Pbe, RhsArray_Pbe, N_U*N_PBE_IntlLevels*SizeOfDouble);
#endif                                        //ifdef __PBS__

    // allocate arrays for solver
    defect = new double[Max_N_Unknowns];

  // parameters for time stepping scheme
  gamma = 0;
  m = 0;
  N_SubSteps = GetN_SubSteps();
  oldtau = 1;
  end_time = TDatabase::TimeDB->ENDTIME;

  // not active : TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL = 0
  if (TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL)
    { time_discs = 2; }
    else
      { time_discs = 1; }


//      GetPtsValuesForIntl(N_PBE_IntlLevels, ScalarPbeFunctions, SolPbe, N_U, PBE_IntlPtValuesT, N_PBE_SpatialPts);
// 
//      L_Nodal2Sol(FeSpace_Intl, PBE_IntlPtValuesT, N_PBE_SpatialPts, N_PBE_IntlLevels, PBE_LSolT);
// 
//      L_Sol2Nodal(FeSpace_Intl, PBE_LSolT, N_PBE_IntlLevels, PBE_IntlPtValuesT, N_PBE_SpatialPts);
// 
//      GetSolFromNodalPtVales(N_PBE_IntlLevels, MaxN_PtsForNodal, ScalarPbeFunctions, SolPbe, 
//                             N_U, PBE_IntlPtValuesT, N_PBE_SpatialPts);
  
// // interpolate exact solution - begin
//     TDatabase::TimeDB->CURRENTTIME = 0.;  
 
//        memset(PBE_IntlPtValuesT, 0, N_PBE_SpatialPts*N_PBE_IntlLevels*SizeOfDouble);
//         for(i=0;i<N_PBE_SpatialPts;i++)
//           {
//            BDUpdatedSolPbe_Intl_Loc = PBE_IntlPtValuesT + i*N_PBE_IntlLevels;
//            ADI_System[i]->Interpolate(BDUpdatedSolPbe_Intl_Loc, Exact_Psd_Intl);
//  
//           }
//          if(TDatabase::ParamDB->REACTOR_P5==0)
//           {
//            GetSolFromNodalPtVales(N_PBE_IntlLevels, MaxN_PtsForNodal, ScalarPbeFunctions, SolPbe, N_U, PBE_IntlPtValuesT, N_PBE_SpatialPts);
// 
//           }
//          else
//           {
//            GetSolFromQuadPtVales(N_PBE_IntlLevels, MaxN_PtsForNodal, ScalarPbeFunctions, SolPbe, N_U, PBE_IntlPtValuesT, N_PBE_SpatialPts);
// 
//           }
//        // interpolate exact solution - end
 
#ifdef __SIMPATURS__
 
  GetOSError(N_PBE_IntlLevels, ScalarPbeFunctions,  FeSpace_Intl, SolPbe, N_U, errors,  Exact_Psd_Intl);

  l2 = errors[0];
  H1 = errors[1];
  N_IntlPts_All = N_PBE_SpatialPts;


  olderror = l2 ;
  olderror1 =  H1 ;
    
    OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
    OutPut(" L2: " << sqrt(l2));
    OutPut(" H1-semi: " << sqrt(H1)<< endl);

    OutPut("Number of Internal Points: " << N_IntlPts_All<<endl);
#endif

//     printf( "Rank : %d,   %e\n", rank,  errors[0]  );

  if(TDatabase::ParamDB->WRITE_VTK)
  {
    //      # ifdef _MPI
    //      if(N_IndepntScalarEqns)
    //       Output->Write_ParVTK(MPI_COMM_WORLD, img, PhySol);
    //
    //      #ifdef __PBS__
    //      if(N_PBEqns)
    //       {
    //        if(ActiveProcess)
    //         GetSolFromNodalPtVales(N_PBE_IntlLevels, MaxN_PtsForNodal, ScalarPbeFunctions, SolPbe, N_U, PBE_IntlPtValuesT, N_PBE_SpatialPts);
    //
    //        for(j=0; j<5; j++)
    //         {
    //          if(j)
    //          {  memcpy(SolPbe_Output, SolPbe+(j*Out_Level -1)*N_U, N_U*SizeOfDouble); }
    //          else
    //           { memcpy(SolPbe_Output, SolPbe+j*N_U, N_U*SizeOfDouble); }
    //
    //          strcpy (SubID, IntLevel);
    //          sprintf (RankBaseName, "%d", j);
    //          strcat(SubID, RankBaseName);
    //          Output_PBS->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
    //         } //for(j=0; j<5; j+
    //
    //        memcpy(SolPbe_Output, SolPbe+N_U, N_U*SizeOfDouble);
    //        strcpy (SubID, IntLevel);
    //        sprintf (RankBaseName, "%d", j);
    //        strcat(SubID, RankBaseName);
    //        Output_PBS->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
    //
    //        } // if(N_PBEqns)
    //       #endif
    //      img++;
    //     #else

      if(N_IndepntScalarEqns)
      {
        os.seekp(std::ios::beg);
        if(img<10) os <<"VTK/"<< VtkBaseName<<".0000"<<img<<".vtk" << ends;
        else if(img<100) os <<"VTK/"<< VtkBaseName<<".000"<<img<<".vtk" << ends;
        else if(img<1000) os <<"VTK/"<< VtkBaseName<<".00"<<img<<".vtk" << ends;
        else if(img<10000) os <<"VTK/"<< VtkBaseName<<".0"<<img<<".vtk" << ends;
        else  os <<"VTK/"<< VtkBaseName<<"."<<img<<".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
      }

//     cout<<"test main " << N_IndepntScalarEqnsendl;
// exit(0);
#ifdef __PBS__
//     j=0;
      for(j=0; j<N_PBE_IntlLevels; j++)
      {
        TDatabase::ParamDB->REACTOR_P29=IntlPosL[j];
        memcpy(SolPbe_Output, SolPbe+j*N_U, N_U*SizeOfDouble);


        os.seekp(std::ios::beg);
        if(img<10) os <<"VTK/"<<  GnuBaseName<<j<<".0000"<<img<<".vtk" << ends;
        else if(img<100) os <<"VTK/"<<  GnuBaseName<<j<<".000"<<img<<".vtk" << ends;
        else if(img<1000) os <<"VTK/"<<  GnuBaseName<<j<<".00"<<img<<".vtk" << ends;
        else if(img<10000) os <<"VTK/"<<  GnuBaseName<<j<<".0"<<img<<".vtk" << ends;
        else  os <<"VTK/"<<  GnuBaseName<<j<<"."<<img<<".vtk" << ends;
        Output_PBS->WriteVtk(os.str().c_str());

//                  //testing
//                  os.seekp(std::ios::beg);
//                  if(img<10) os <<"VTK/"<<  GnuBaseName<<".0000"<<j+1<<".vtk" << ends;
//                  else if(img<100) os <<"VTK/"<<  GnuBaseName<<".000"<<j+1<<".vtk" << ends;
//                  else if(img<1000) os <<"VTK/"<<  GnuBaseName<<".00"<<j+1<<".vtk" << ends;
//                  else if(img<10000) os <<"VTK/"<<  GnuBaseName<<".0"<<j+1<<".vtk" << ends;
//                  else  os <<"VTK/"<<  GnuBaseName<<"."<<j+1<<".vtk" << ends;
//                  Output_PBS->WriteVtk(os.str().c_str());
      }
#endif
      img++;
  }

#ifdef __PBS__  
        BDUpdatedSolPbe_Intl_Loc = PBE_IntlPtValuesT + f_Integ_NodalPt[0]*N_LDof;

       q3_max=ADI_System[f_Integ_NodalPt[0]]->GetQ3Max(BDUpdatedSolPbe_Intl_Loc);
       OutPut( endl);
       OutPut( "q3_max " <<  setprecision(8) << q3_max << endl);
#endif
    
  
  
#ifdef __PBS__
    l_infty = TDatabase::ParamDB->UREA_l_infty;
    u_infty = TDatabase::ParamDB->UREA_u_infty;
    f_infty = TDatabase::ParamDB->UREA_f_infty; 

    if (TDatabase::ParamDB->PBE_P0 == 2 || TDatabase::ParamDB->PBE_P0 == 4 )
     { 
      aggr_time = 0;
      
       
      layers_mass = new double[N_PBE_IntlLevels];
      
      for (i=0 ; i<N_PBE_IntlLevels ; i++)
       {
        layers_mass[i] = pow(IntlPosL[i],3.0);
//      cout << "layers_mass " << layers_mass[i]<< endl;
       }

/*#ifdef __AXIAL3D__  
      params[0]= 3;//TDatabase::ParamDB->WINDTUNNEL_SPATIAL;//space : dim      
#else    */  
      params[0]= 2;//TDatabase::ParamDB->WINDTUNNEL_SPATIAL;//space : dim
// #endif      
      params[1]= 1;//TDatabase::ParamDB->WINDTUNNEL_BROWNIAN; // brown kernel included or not
      params[2]= 1;//TDatabase::ParamDB->WINDTUNNEL_POL_ORDER;//pol or constant 0:const 1:pol
      params[3]= 1;//TDatabase::ParamDB->WINDTUNNEL_BROWNIAN_TEMP;//brown kernel depends of
      //temperature or not no:0 yes:
      params[4]= TDatabase::ParamDB->UREA_AGGR_BROWNIAN;//scal param for brown kernel
//       OutPut("brown: "<< params[4]<<endl);
      params[5]=1;//TDatabase::ParamDB->WINDTUNNEL_SHEAR_FACTOR_TYPE;//shear depends of velocity or not
      params[6]= TDatabase::ParamDB->UREA_AGGR_SHEAR_FACTOR;// scal param for shear kernel
      params[7]= 1.074e-3;  // dynamic viscosity of the fluid
          
     } //if (TDatabase::ParamDB->PBE_P0

     
 
     
     if(!UpdateConvection)
      {
       GetVeloGradForIntlPts(Velocity_Spaces, velo1, velo2, Scalar_Spaces[PBE_INDEX_LOW], Velo_PBE_SpatialPts, 
                             GradVelo_PBE_SpatialPts, N_PBE_SpatialPts);          

        Dscal(2*N_PBE_SpatialPts, u_infty, Velo_PBE_SpatialPts); 
        Dscal(4*N_PBE_SpatialPts, u_infty/l_infty, GradVelo_PBE_SpatialPts); //  u_infty/l_infty =1, so no need
       }
 
// assemble the mass matrix for the internal space once for all time steps
  if(TDatabase::ParamDB->INTL_DISCTYPE==FD)
   { AssembleFD_M_StartPtDirichlet(M_Intl); }
  else
   { AssembleM_StartPtDirichlet(FeSpace_Intl, M_Intl);}

   /** store the M Mat values, so that no need to restore after every sptl pt. calculatio*/
   i = M_Intl->GetN_Entries();
   MatValues = new double[i];
   memcpy(MatValues,  M_Intl->GetEntries(),  i*SizeOfDouble);

    // get the initial rhs for PBE
   /** initial growth uses the c,T after spatial update., so it is really the first interation mat in the l update */
//     for(i=0;i<N_PBE_SpatialPts;i++)
//     {
//       if(DirichletBDPt[i])
//         continue;
// 
// 
//       /** Restrict sol from high fespace to low fespace */
//        if(PBE_ORDER != ORDER)
//        {
// #ifndef __PBSConstT__   
//         RestrictFunction(Scalar_Spaces_low[0], Scalar_Spaces[0], SolArray_low[0], SolArray[0], PbeAux);
// 
//         RestrictFunction(Scalar_Spaces_low[1], Scalar_Spaces[1], SolArray_low[1], SolArray[1], PbeAux);
// #endif
//        }
// 
// 
//        if( TDatabase::ParamDB->SC_VERBOSE > 0 )
//           t2 = GetTime();
	
        // get the init concentration
        // T, C, FESpace_PBE, C, C_sat, N_PBE_SpatialPts
        // get the current timestep BD updated PB for all internal levels
//         if(TDatabase::ParamDB->REACTOR_P5==0)
//         {
//           if(N_IndepntScalarEqns)
//             GetPtsValuesForIntl(ScalarFunctions_low[0], ScalarFunctions_low[1], Scalar_Spaces[PBE_INDEX_LOW],
//                                 C_IntlPtValues, T_IntlPtValues, C_Sat_IntlPtValues, N_PBE_SpatialPts);
// 
//           /** transfer sol from spatial to internal */
//           GetPtsValuesForIntl(N_PBE_IntlLevels, ScalarPbeFunctions_low, SolPbe_low, N_U_low, PBE_IntlPtValuesT, N_PBE_SpatialPts);
// 
//           /** evaluate nodal values from nodal points */
//           L_Nodal2Sol(FeSpace_Intl, PBE_IntlPtValuesT, N_PBE_SpatialPts, N_PBE_IntlLevels, PBE_LSolT);
//         }
//         else
//         {
//           if(N_IndepntScalarEqns)
// 	   {
//             GetQuadPtsValuesForIntl(ScalarFunctions[0], ScalarFunctions[1], Scalar_Spaces[PBE_INDEX],  C_IntlPtValues,  C_Sat_IntlPtValues, N_PBE_SpatialPts);
// 	    cout << " T_IntlPtValues for heat coupling with PBE is not yet implemented for QuadPts" <<endl;
// 	    exit(0);
// 	   }
// 	   /** transfer sol from spatial to internal */
//           GetPtsValuesForIntl_Quad(N_PBE_IntlLevels, ScalarPbeFunctions_low,  PBE_IntlPtValuesT, N_PBE_SpatialPts);
//         }
// 
//       if(N_IndepntScalarEqns)
//           {
//             C = C_IntlPtValues[i];
//             C_Sat = C_Sat_IntlPtValues[i];
// 	    T = T_IntlPtValues[i];
//             //  cout<< " C " << C << " C_Sat " << C_Sat;
//           }
//       else
//           {
//             T=0;
//             C = 0.;
//             C_Sat = 0.;
//           }
// 
//      /** assemble the Growth dept. intl A and SUPG matrices, for rhs calculation in the first time step */
//      ADI_System[i]->AssembleInitRhs(C, T, BilinearCoeffs_Psd_Intl, cond_Lmin, cond_Lmax);
//     }  //  for(i=0;i<N_PBE_SpatialPts;i++)


//   RowPtr = M_Intl->GetRowPtr();
//   KCol = M_Intl->GetKCol();
//   MatValues = M_Intl->GetEntries();

  coll->GetHminHmax(&hmin,&hmax);
  //     printf( "Rank : %d, h_min :  %e,  h_max : %e\n", rank, hmin, hmax);

//   tau = pow(hmin, TDatabase::ParamDB->REACTOR_P30);
#endif

  
#ifdef __PBSConstT__
  OutPut("TIMESTEPLENGTH : " << tau << endl;)
  // exit(0);
  TDatabase::TimeDB->TIMESTEPLENGTH = tau;
#endif // #ifdef __PBSConstT__


  OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
  OutPut("TIMESTEPLENGTH : " << TDatabase::TimeDB->TIMESTEPLENGTH << endl;)
 

#ifndef __PBSConstT__
#ifdef __PBS__
   CompOutletAvgFunct(SolPbe, N_U, PbeFunctions, Out_CellIndex, Out_EdgeIndex, N_OutEdges);
#endif 
#endif 


#ifdef __PBSConstT__
 UpdatePBERhs = TRUE;
#endif   
     OutPut("Solution f_Init " << Ddot(N_PBE_SpatialPts*N_LDof, SolPbe, SolPbe)<<endl);
//      exit(0);  

     // WriteQ3AtInlet
       BDUpdatedSolPbe_Intl_Loc = PBE_IntlPtValuesT+(f_Inlet_SpaceIndex*N_LDof);
       q3_max=ADI_System[f_OutletIndex[0]]->GetQ3Max(BDUpdatedSolPbe_Intl_Loc);
     
        os.seekp(std::ios::beg);
        os << "Q3Data/InletQ3-0"<<N_Q3Data<<".data" << ends;
        std::ofstream dat(os.str().c_str());
        dat.setf(std::ios::fixed);
        dat << setprecision(12);

        if (!dat)
         {
          cerr << "cannot open file for output" << endl;
          return -1;
         }
        dat << "%% Outlet Avg. Q3 data created for droplet by ParMooN" << endl;
        dat << "%% Current Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
        dat << "%% x=200 line at outlet"<< endl;
        dat << "%% q3_max : " <<  setprecision(12) << q3_max << endl;
        dat << "%% L , Avg.f, Avg. q3: " << endl;	 


        if( fabs(q3_max)>1.0e-100)
         {
           // N_LDof should be equal to N_PBE_IntlLevels !!
          for(j=0; j<N_PBE_IntlLevels; j++)
           {
            L=IntlPosL[j];
            dat <<L <<  " " <<   setprecision(12) << BDUpdatedSolPbe_Intl_Loc[j] << " " <<
                        setprecision(12) <<  L*L*L*BDUpdatedSolPbe_Intl_Loc[j]/q3_max<<endl;
           }
         }
        else
         {
          for(j=0; j<N_PBE_IntlLevels; j++)
           {
            L=IntlPosL[j];
            dat <<L << " " <<  BDUpdatedSolPbe_Intl_Loc[j] << " " << 0. <<endl;
           }
         }

        dat.close();
        cout << endl;
        OutPut( "OutletAvgQ3 data wrote into file Q3Data/OutletAvgQ3*" << N_Q3Data << ".data" <<endl);
        N_Q3Data++;     
     
//    exit(0);
   
  //======================================================================
  // start of time cycle
  //======================================================================
  tau = TDatabase::TimeDB->CURRENTTIME;
  while(TDatabase::TimeDB->CURRENTTIME< end_time)
  {                                               // time cycle
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

    for(l=0;l<N_SubSteps;l++)                     // sub steps of fractional step theta
    {
      
     pbs_t1 = GetTime();    

      if (!very_first_time)
      {
        SetTimeDiscParameters();
      }

      if(m==1)
      {
        OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
        OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
        OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
        OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
      }

      tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      if (!very_first_time)
        TDatabase::TimeDB->CURRENTTIME += tau;

 
      
//      if(TDatabase::TimeDB->CURRENTTIME>1500 && TDatabase::TimeDB->CURRENTTIME<=1501)
//       {
//        TDatabase::ParamDB->REACTOR_P24 =283.15+18.*(TDatabase::TimeDB->CURRENTTIME-1500);
//        UpdateRhs = TRUE;
//       } 
//      else if (TDatabase::TimeDB->CURRENTTIME>3500 && TDatabase::TimeDB->CURRENTTIME<=3501)
//       {
//        TDatabase::ParamDB->REACTOR_P24 = 301.15-18.*(TDatabase::TimeDB->CURRENTTIME-3500);
//        UpdateRhs = TRUE;
//       }

//         OutPut("T_Wall " << TDatabase::ParamDB->REACTOR_P24 <<endl);
        OutPut(endl << "CURRENT TIME: ");
        OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

      for(i=0;i<N_IndepntScalarEqns;i++)
      {
        // working array for rhs is B, initialize B
        memset(B[i], 0, N_Uarray[i]*SizeOfDouble);

        // compute terms with data from previous time step
        // old rhs multiplied with current subtime step and theta3 on B
        Daxpy(N_Active[i], tau*TDatabase::TimeDB->THETA3, RhsArray[i], B[i]);

          if(UpdateConvection || UpdateRhs ||  ConvectionFirstTime )
          {
            rb_flag = 0;
            // assemble A and rhs
            N_Rhs = 1;
            N_FESpaces = 1;
            fesp[0] = Scalar_Spaces[i];

            if (VeloFunction)
            {
              N_FESpaces = 2;
              fesp[1] = Velocity_Spaces;
              fefct[0] = velo1;
              fefct[1] = velo2;

              // defined in TCD2D.h
              aux =  new TAuxParam2D(TimeCDParamsVeloFieldN_FESpaces,
                TimeCDParamsVeloFieldN_Fct,
                TimeCDParamsVeloFieldN_ParamFct,
                TimeCDParamsVeloFieldN_FEValues,
                fesp+1, fefct,
                TimeCDParamsVeloFieldFct,
                TimeCDParamsVeloFieldFEFctIndex,
                TimeCDParamsVeloFieldFEMultiIndex,
                TimeCDParamsVeloFieldN_Params,
                TimeCDParamsVeloFieldBeginParam);

            }
            else
            {
              aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 1, NULL);     
      
            }

            N_SquareMatrices = 1;
            SQMATRICES[0] = MatricesA[i];
            SQMATRICES[0]->Reset();
            DiscreteForm = DiscreteFormMatrixARhs[i];

            if(Disctypes[i] == SDFEM)
            {
              N_SquareMatrices = 2;
              SQMATRICES[1] = MatricesK[i];
              SQMATRICES[1]->Reset();
              DiscreteForm = DiscreteFormMatrixARhs_SUPG[i];
            }

            BDCond[0] = BoundaryConditions[i];
            BDValue[0] = BoundValues[i];

            memset(RhsArray[i], 0, N_Uarray[i]*SizeOfDouble);
            RHSs[0] = RhsArray[i];
            ferhs[0] = Scalar_Spaces[i];

            Assemble2D(N_FESpaces, fesp,
              N_SquareMatrices, SQMATRICES,
              0, NULL,
              N_Rhs, RHSs, ferhs,
              DiscreteForm,
              BDCond,
              BDValue,
              aux);
            delete aux;

#ifdef __ROBINBC__    
          // update the LHS in the Robin BC to the stiffness matrix     
          RobinBCInt(SQMATRICES[0], BDCond[0]);

#endif
    
            if(i+1==N_IndepntScalarEqns)
             {
              ConvectionFirstTime = FALSE;
              UpdateRhs = FALSE; 
             }
          }

          // add rhs from current sub time step to rhs array B
          Daxpy(N_Active[i], tau*TDatabase::TimeDB->THETA4,  RhsArray[i], B[i]);

          // update rhs by Laplacian and convective term from previous
          // time step
          // scaled by current sub time step length and theta2
          // currently : M := M + gamma A
          // M = M + (- tau*TDatabase::TimeDB->THETA2)
          MatAdd(MatricesM[i], MatricesA[i], - tau*TDatabase::TimeDB->THETA2);
          // set current factor of steady state matrix
          gamma = -tau*TDatabase::TimeDB->THETA2;

          if(Disctypes[i] == SDFEM)
            MatAdd(MatricesM[i], MatricesK[i], 1);

          // defect = M * sol
          // B:= B + defec
          memset(defect, 0, Max_N_Unknowns*SizeOfDouble);
          MatVectActive(MatricesM[i], SolArray[i], defect);
          Daxpy(N_Active[i], 1, defect, B[i]);

          // set Dirichlet values
          // RHSs[0] still available from assembling
          memcpy(B[i]+N_Active[i], RhsArray[i]+N_Active[i], (N_Uarray[i]-N_Active[i])*SizeOfDouble);
          // copy Dirichlet values from rhs into sol
          memcpy(SolArray[i]+N_Active[i], RhsArray[i]+N_Active[i], (N_Uarray[i]-N_Active[i])*SizeOfDouble);

          // system matrix
          MatAdd(MatricesM[i], MatricesA[i], -gamma + tau*TDatabase::TimeDB->THETA1);
          gamma = tau*TDatabase::TimeDB->THETA1;

        //======================================================================
        // solve linear system
        //======================================================================


          if( TDatabase::ParamDB->SC_VERBOSE > 0 )
            t1 = GetTime();

        // solve the system with a parallel solver
        switch(int(TDatabase::ParamDB->SOLVER_TYPE))
        {
          case 0:
            // AMG Solver
            cout << "solver type not implemented yet !!!!!!!!!!!!!" << endl;
            exit(0);
            break;

          case 1:
            // GMG Solver
            cout << "solver type not implemented yet !!!!!!!!!!!!!" << endl;
            exit(0);

            break;

          case 2:
//             DirectSolver(MatricesM[i], B[i], SolArray[i]);
// 	      rb_flag = 3;
               DirectSolver(MatricesM[i], B[i], SolArray[i], MatValues_Solver[i], 
                            KCol_Solver[i], Row_Solver[i], Symbolic_Solver[i], Numeric_Solver[i], rb_flag);
            break; 
#ifdef _OMP
          case 100:

            if(Initialize_ScalarSolver)
            {
              /** Construct SMP Direct Solver for the scalar **/
              // Sort_ForDirectSolver() should be called before this construction direct solver!
              // collect data of the global system

              Scalar_Entries = MatricesM[i]->GetEntries();
              N_Eqn = N_Uarray[i];
              ScalarKCol = sqstructureA[i]->GetKCol();
              ScalarRowPtr = sqstructureA[i]->GetRowPtr();
              N_Entries = ScalarRowPtr[N_Eqn];
              KCol = new int[N_Entries];
              RowPtr = new int[N_Eqn+1];

              RowPtr[0] = 0 + 1;                  // fortran style for ParDiso direct solver
              Eq_pos = 0;
              for(j=0;j<N_Eqn;j++)
              {
                begin = ScalarRowPtr[j];
                end = ScalarRowPtr[j+1];
                for(k=begin;k<end;k++)
                {
                                                  // fortran style
                  KCol[Eq_pos] = ScalarKCol[k] + 1;
                  Eq_pos++;
                }

                RowPtr[j+1] = Eq_pos + 1;         // fortran style for ParDiso direct solver
              }

              Scalar_SMPSolver[i] = new TParDirectSolver(N_Eqn, RowPtr, KCol, N_Entries);
              Scalar_SMPSolver[i]->SymbolicFactorize(Scalar_Entries);

              Initialize_ScalarSolver = FALSE;
              delete [] KCol; delete [] RowPtr;
            }                                     // if(Initialize_ScalarSolver)

            // Scalar_SMPSolver[i]->FactorizeAndSolve(MatricesM[i], B[i], SolArray[i]);
            Scalar_SMPSolver[i]->Solve(MatricesM[i], B[i], SolArray[i]);

            break;
#endif

          default:
            cout << "wrong parallel solver type !!!!!!!!!!!!!" << endl;
            exit(0);
            break;
        }

 
       if( TDatabase::ParamDB->SC_VERBOSE > 0 )
        {
          t2 = GetTime();
          OutPut( "time for solving: " << t2-t1 << endl);
        }
        //======================================================================
        // end solve linear system
        //======================================================================
          // restore matrices
          MatAdd(MatricesM[i], MatricesA[i], -gamma);
          gamma = 0;

          if(Disctypes[i] == SDFEM)
            MatAdd(MatricesM[i], MatricesK[i], -1);
      }     // for(i=0;i<N_IndepntScalarEqns;i++)

     // IndepntScalarEqns matrices are same over the time, so no need to LU decomp evey timestep
     if(rb_flag==0) rb_flag =1;    
      
#ifdef __PBS__
      //===============================================================================
      // PBE system solution -- start
      // assume that the concentration eqn solved
      // first solve w.r.t  w.r.t physical space and then onternal space (L-direction)
      // PBE system X-direction solution --start
      //===============================================================================
      //set the boundary value at L_min and L_max for the current timestep
      // not necessar here, as DIRICHLET BC at L_min is not solved spatially 
//         if(!dGDisc)
//          SetLMinLMaxBoundValue(N_PBE_IntlLevels, ScalarPbeFunctions, N_U, RhsArray_Pbe);

        //set the boundary value at Dirichlet nodes for all levels in the current timestep
        for(i=start_pbe; i<end_pbe; i++)
        {
          if(UpdatePBEConvection || UpdatePBERhs ||  PBEConvectionFirstTime )
          {
            TDatabase::ParamDB->REACTOR_P29=IntlPosL[i];

            // set parameters
            N_Rhs = 1;
            N_FESpaces = 1;
            fesp[0] = Scalar_Spaces[PBE_INDEX];
            ferhs[0] = Scalar_Spaces[PBE_INDEX];

            // reset matrices
            N_SquareMatrices = 1;
            SQMATRICES[0] = MatricesA[PBE_INDEX];
            SQMATRICES[0]->Reset();
            DiscreteForm = DiscreteFormMatrixARhs[PBE_INDEX];

            if(Disctypes[PBE_INDEX] == SDFEM)
            {
              N_SquareMatrices = 2;
              SQMATRICES[1] = MatricesK[PBE_INDEX];
              SQMATRICES[1]->Reset();
              DiscreteForm = DiscreteFormMatrixARhs_SUPG[PBE_INDEX];
            }

            BDCond[0] = BoundaryConditions[PBE_INDEX];
            BDValue[0] = BoundValues[PBE_INDEX];

            // assemble A and rhs
            if(VeloFunction)
            {
              N_FESpaces = 2;
              fesp[1] = Velocity_Spaces;
              fefct[0] = velo1;
              fefct[1] = velo2;

              // defined in TimeConvDiff2D.h
              aux =  new TAuxParam2D(TimeCDParamsVeloFieldN_FESpaces,
                TimeCDParamsVeloFieldN_Fct,
                TimeCDParamsVeloFieldN_ParamFct,
                TimeCDParamsVeloFieldN_FEValues,
                fesp+1, fefct,
                TimeCDParamsVeloFieldFct,
                TimeCDParamsVeloFieldFEFctIndex,
                TimeCDParamsVeloFieldFEMultiIndex,
                TimeCDParamsVeloFieldN_Params,
                TimeCDParamsVeloFieldBeginParam);

            }
            else
              { aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 1, NULL); }

              memset(RhsArray_Pbe+i*N_U, 0, N_U*SizeOfDouble);
            RHSs[0] = RhsArray_Pbe+i*N_U;

            Assemble2D(N_FESpaces, fesp,
              N_SquareMatrices, SQMATRICES,
              0, NULL,
              N_Rhs, RHSs, ferhs,
              DiscreteForm,
              BDCond,
              BDValue,
              aux);

            delete aux;

            PBEConvectionFirstTime = FALSE;
            rb_flag_pbe = 0;
          }  // if(UpdatePBEConvection || UpdatePBERhs ||  PBEConvectionFirstTime )

          // copy Dirichlet values from rhs into sol
          SolPbe_OneLevel = SolPbe+i*N_U;
          memcpy(SolPbe_OneLevel+N_Active[PBE_INDEX], RhsArray_Pbe+(i*N_U +N_Active[PBE_INDEX]),
            (N_Uarray[PBE_INDEX]-N_Active[PBE_INDEX])*SizeOfDouble);
        }  // for(i=start_pbe; i<end_pbe; i++)

        memcpy(SolPbe_Lmin, SolPbe, N_U*SizeOfDouble);
        memcpy(SolPbe_Lmax, SolPbe+(N_U*(N_PBE_IntlLevels-1)), N_U*SizeOfDouble);
        //====================================================================================

        gamma=0.;
        // working array for rhs is B, initialize B
        // same lhs for alll levels
        memset(B_Pbe, 0, N_U*N_PBE_IntlLevels*SizeOfDouble);

        for(i=0; i<N_PBE_IntlLevels; i++)
        {
          // add rhs from current sub time step to rhs array B
          Daxpy(N_Active[PBE_INDEX], tau*TDatabase::TimeDB->THETA3, OldRhsArray_Pbe+i*N_U, B_Pbe+i*N_U);

          Daxpy(N_Active[PBE_INDEX], tau*TDatabase::TimeDB->THETA4, RhsArray_Pbe+i*N_U, B_Pbe+i*N_U);

          if(gamma==0.)
          {
            MatAdd(MatricesM[PBE_INDEX], MatricesA[PBE_INDEX], -tau*TDatabase::TimeDB->THETA2);
            // set current factor of steady state matrix
            gamma = -tau*TDatabase::TimeDB->THETA2;
          }

          SolPbe_OneLevel = SolPbe+i*N_U;

          // copy Dirichlet values from rhs into sol
          memcpy(SolPbe_OneLevel+N_Active[PBE_INDEX], RhsArray_Pbe+(i*N_U +N_Active[PBE_INDEX]),
            (N_Uarray[PBE_INDEX]-N_Active[PBE_INDEX])*SizeOfDouble);

          memset(defect, 0, Max_N_Unknowns*SizeOfDouble);
          MatVectActive(MatricesM[PBE_INDEX], SolPbe_OneLevel, defect);
          Daxpy(N_Active[PBE_INDEX], 1, defect, B_Pbe+i*N_U);

          // set Dirichlet values
          memcpy(B_Pbe+(i*N_U+N_Active[PBE_INDEX]), RhsArray_Pbe+(i*N_U + N_Active[PBE_INDEX]),
            (N_Uarray[PBE_INDEX]-N_Active[PBE_INDEX])*SizeOfDouble);
        }                                         // for(i=0; i<N_PBE_IntlLevels; i++)

        memcpy(OldRhsArray_Pbe, RhsArray_Pbe, N_U*N_PBE_IntlLevels*SizeOfDouble);

        // system matrix
        MatAdd(MatricesM[PBE_INDEX], MatricesA[PBE_INDEX], -gamma + tau*TDatabase::TimeDB->THETA1);
        gamma = tau*TDatabase::TimeDB->THETA1;

      //======================================================================
      // solve linear PB system with multiple rhs
      //======================================================================
      if( TDatabase::ParamDB->SC_VERBOSE > 0 )
        t1 = GetTime();

      // solve the system with a parallel solver
      switch(int(TDatabase::ParamDB->SOLVER_TYPE))
      {
        case 0:
          // AMG Solver
          cout << "solver type not implemented yet !!!!!!!!!!!!!" << endl;
          exit(0);
          break;

        case 1:
          // GMG Solver
          cout << "solver type not implemented yet !!!!!!!!!!!!!" << endl;
          exit(0);
          break;

        case 2:
          // solve linear PB system with multiple rhs
          SQMATRICESSOLVER[0] = MatricesM[PBE_INDEX];
//           DirectSolver(SQMATRICESSOLVER[0], B_Pbe, SolPbe, end_pbe, start_pbe);
          DirectSolver(SQMATRICESSOLVER[0], B_Pbe, SolPbe, end_pbe, start_pbe, MatValues_Solver[PBE_INDEX], 
                       KCol_Solver[PBE_INDEX], Row_Solver[PBE_INDEX], Symbolic_Solver[PBE_INDEX], Numeric_Solver[PBE_INDEX], rb_flag_pbe);
          break;
#ifdef _OMP
        case 100:
          cout << "solver type not implemented yet for multiple rhs !!!!!!!!!!!!!" << endl;
          exit(0);
          break;
#endif

        default:
          cout << "wrong parallel solver type !!!!!!!!!!!!!" << endl;
          exit(0);
          break;
      }

     if( TDatabase::ParamDB->SC_VERBOSE > 0 )
      {
        t2 = GetTime();
        OutPut( "time for solving: " << t2-t1 << endl);
      }
      //======================================================================
      // end solve linear system
      //======================================================================
        // restore matrices
        MatAdd(MatricesM[PBE_INDEX], MatricesA[PBE_INDEX], -gamma);
        gamma = 0;

     // x-dir PBE matrix is same over the time, so no need to LU decomp evey timestep
     if(rb_flag_pbe==0) rb_flag_pbe =1;    	

     pbs_t2 = GetTime(); 
     pbs_step1 += pbs_t2 - pbs_t1;
      //===============================================================================
      //PBE system X-direction solution -- end
      //PBE system L-direction solution -- start  
      //===============================================================================
     if(m%100==0 || m==1)
      {     
       t1 = GetTime();  
      }     
     
    pbs_t1 = GetTime();   
     
     
   if(TDatabase::ParamDB->PBE_P0!=0)
    {
      /** Restrict sol from high fespace to low fespace */
       if(PBE_ORDER != ORDER)
       {
        RestrictFunction(Scalar_Spaces[PBE_INDEX_LOW], Scalar_Spaces[PBE_INDEX], N_PBE_IntlLevels,
               SolPbe_low, SolPbe, PbeAux);

#ifndef __PBSConstT__   
        RestrictFunction(Scalar_Spaces_low[0], Scalar_Spaces[0], SolArray_low[0], SolArray[0], PbeAux);

        RestrictFunction(Scalar_Spaces_low[1], Scalar_Spaces[1], SolArray_low[1], SolArray[1], PbeAux);
#endif
       }


//        if( TDatabase::ParamDB->SC_VERBOSE > 0 )
//           t2 = GetTime();
	
        // get the current timestep concentration
        // T, C, FESpace_PBE, C, C_sat, N_PBE_SpatialPts
        // get the current timestep BD updated PB for all internal levels
        if(TDatabase::ParamDB->REACTOR_P5==0)
         {
          if(N_IndepntScalarEqns)
            GetPtsValuesForIntl(ScalarFunctions_low[0], ScalarFunctions_low[1], Scalar_Spaces[PBE_INDEX_LOW],
                                C_IntlPtValues, T_IntlPtValues, C_Sat_IntlPtValues, N_PBE_SpatialPts);

          /** transfer sol from spatial to internal */
          GetPtsValuesForIntl(N_PBE_IntlLevels, ScalarPbeFunctions_low, SolPbe_low, N_U_low, PBE_IntlPtValuesT, N_PBE_SpatialPts);

          /** evaluate nodal values from nodal points */
          L_Nodal2Sol(FeSpace_Intl, PBE_IntlPtValuesT, N_PBE_SpatialPts, N_PBE_IntlLevels, PBE_LSolT);
         }
        else
         {
          if(N_IndepntScalarEqns)
	   {
            GetQuadPtsValuesForIntl(ScalarFunctions[0], ScalarFunctions[1], Scalar_Spaces[PBE_INDEX],  C_IntlPtValues,  C_Sat_IntlPtValues, N_PBE_SpatialPts);
	    cout << " T_IntlPtValues for heat coupling with PBE is not yet implemented for QuadPts" <<endl;
	    exit(0);
	   }
	   /** transfer sol from spatial to internal */
          GetPtsValuesForIntl_Quad(N_PBE_IntlLevels, ScalarPbeFunctions_low,  PBE_IntlPtValuesT, N_PBE_SpatialPts);
         }
 
     
        /* compute rhs */
        if (TDatabase::ParamDB->PBE_P0 == 2 || TDatabase::ParamDB->PBE_P0 == 4 )
        {
          pbs_t3 = GetTime();
          OutPut("===== PSD with Aggregation ======="<< endl);

       /** stored old Aggregation for rhs */
//        if(TDatabase::TimeDB->THETA3)
//         Daxpy(N_LDof*N_PBE_SpatialPts, tau*TDatabase::TimeDB->THETA3, PBE_IntlPtRhsValuesT, OldPBE_IntlPtRhsValuesT); 



         Trans_XLtoLX(PBE_IntlPtValuesT, PBE_IntlPtValues, N_PBE_SpatialPts, N_PBE_IntlLevels);

         if(UpdateConvection)   
          {
           GetVeloGradForIntlPts(Velocity_Spaces, velo1, velo2, Scalar_Spaces[PBE_INDEX], Velo_PBE_SpatialPts, GradVelo_PBE_SpatialPts, N_PBE_SpatialPts);


// #ifdef __AXIAL3D__      
//            Dscal(3*N_PBE_SpatialPts, u_infty, Velo_PBE_SpatialPts);  
//            Dscal(9*N_PBE_SpatialPts, u_infty/l_infty, GradVelo_PBE_SpatialPts); //  u_infty/l_infty =1, so no need       
// #else
           Dscal(2*N_PBE_SpatialPts, u_infty, Velo_PBE_SpatialPts);  
           Dscal(4*N_PBE_SpatialPts, u_infty/l_infty, GradVelo_PBE_SpatialPts); //  u_infty/l_infty =1, so no need
//            memset(GradVelo_PBE_SpatialPts, 0, 4*N_PBE_SpatialPts*SizeOfDouble);
// #endif
          }
          
          Dscal(N_PBE_SpatialPts*N_LDof, f_infty, SolPbe_low);
  
//           OutPut("velo " << Ddot(2*N_PBE_SpatialPts, Velo_PBE_SpatialPts, Velo_PBE_SpatialPts)<<endl);
//           OutPut("GradVelo " << Ddot(4*N_PBE_SpatialPts, GradVelo_PBE_SpatialPts, GradVelo_PBE_SpatialPts)<<endl);
//           OutPut("temp " << Ddot(N_PBE_SpatialPts, T_IntlPtValues, T_IntlPtValues)<<endl);
//           OutPut("Conc " << Ddot(N_PBE_SpatialPts, C_IntlPtValues, C_IntlPtValues)<<endl);
//           OutPut("f " << Ddot(N_PBE_SpatialPts*N_LDof, SolPbe_low, SolPbe_low)<<endl);

         memset(PBE_IntlPtRhsValues, 0, N_PBE_SpatialPts*N_PBE_IntlLevels*SizeOfDouble);

         apply_integral_operator(N_PBE_SpatialPts, 1, 1, N_PBE_IntlLevels, SolPbe_low,
                                 Velo_PBE_SpatialPts, GradVelo_PBE_SpatialPts, T_IntlPtValues,
                                 PBE_IntlPtRhsValues, layers_mass,
                                 TDatabase::ParamDB->UREA_D_P_MAX,
                                 TDatabase::ParamDB->UREA_f_infty,
                                 params
                                );
  
          Dscal(N_PBE_SpatialPts*N_LDof, 1./f_infty, SolPbe_low);

         /**if aggr needs to be considered implicitly, then Aggr function has to be called within nonlinear loop */
         /**so, aggregation is treated explicitly, as the aggr routies are expensive, so dt is mult. here  */
          Dscal(N_PBE_SpatialPts*N_LDof, tau*l_infty/(u_infty*f_infty), PBE_IntlPtRhsValues);  
	  
         if(TDatabase::ParamDB->REACTOR_P5==0)
          {
           /** transfer rhs from spatial to internal */
           GetPtsValuesForIntl(N_PBE_IntlLevels, ScalarPbeFunctions_low, PBE_IntlPtRhsValues, N_U_low, PBE_IntlPtRhsValuesT, N_PBE_SpatialPts);

          // not necessary for P1 element in internal space, otherwise Aggregation has to be calculated accordingly for each internal quad. pt.
          /** evaluate nodal values from nodal points */
          // L_Nodal2Sol(FeSpace_Intl, OldPBE_IntlPtRhsValuesT, N_PBE_SpatialPts, N_PBE_IntlLevels, PBE_IntlPtRhsValuesT);
          pbs_t4 = GetTime();
          aggr_time += pbs_t4 - pbs_t3; 
         }
//           OutPut("Aggregation " << Ddot(N_PBE_SpatialPts*N_LDof, PBE_IntlPtRhsValues, PBE_IntlPtRhsValues)<<endl);
//           OutPut("===== PSD with Aggregation end======="<<endl);
//           exit(0);
        }
//      cout << "test " <<endl;
//      exit(0);
//     OutPut(N_PBE_SpatialPts << " Solution f_x " << Ddot(N_PBE_SpatialPts*N_LDof, PBE_LSolT, PBE_LSolT)<<endl);
// //  #pragma omp parallel for

//        t1 =  Ddot(N_PBE_SpatialPts*N_LDof, PBE_LSolT, PBE_LSolT);
        for(i=0;i<N_PBE_SpatialPts;i++)
        {
          if(DirichletBDPt[i])
            continue;

//      i =20;
//           BDUpdatedSolPbe_Intl_Loc = PBE_IntlPtValuesT + i*N_PBE_IntlLevels;
          RhsPbe_Intl_Loc = PBE_IntlPtRhsValuesT + i*N_LDof;
          BDUpdatedSolPbe_Intl_Loc = PBE_LSolT + (i*N_LDof);

          // old rhs, Agg + [A+ (dt\theta_3)M]f^{n-1}, when the meatrix part needs to return from the SolveAllLevels   
//           RhsPbe_Intl_Loc = OldPBE_IntlPtRhsValuesT + (i*N_LDof);  // Mat part will be stored after each SolveAllLevels, and then Agg

//           memset(Oldrhs_IntlLoc, 0, N_LDof*SizeOfDouble);
//           Daxpy(N_LDof, 1, RhsPbe_Intl_Loc, Oldrhs_IntlLoc);

//           if(TDatabase::TimeDB->THETA4)
//            Daxpy(N_LDof, tau*TDatabase::TimeDB->THETA4, PBE_IntlPtRhsValuesT + (i*N_LDof), Oldrhs_IntlLoc);

          if(N_IndepntScalarEqns)
          {
            C = C_IntlPtValues[i];
            C_Sat = C_Sat_IntlPtValues[i];
	    T = T_IntlPtValues[i];
            //  cout<< " C " << C << " C_Sat " << C_Sat;
          }
          else
          {
            T=0;    
            C = 0.;
            C_Sat = 0.;
          }

//           /** set Dirichlet BC at L_min */
//           if(!dGDisc)  done in internal level solver
//            { BDUpdatedSolPbe_Intl_Loc[0] = LMinBDValue(C, T); }

//           cout<< i << " C_Sat_IntlPtValues " << C_IntlPtValues[i] << endl;
          ADI_System[i]->SolveAllLevels(MatValues, BDUpdatedSolPbe_Intl_Loc, RhsPbe_Intl_Loc, 
                                        C, C_Sat, T, BilinearCoeffs_Psd_Intl, tau, 
                                        cond_Lmin, cond_Lmax, NULL); //Exact_Psd_Intl
	  
// 	  t2 = Ddot(N_PBE_SpatialPts*N_LDof, PBE_LSolT, PBE_LSolT);
// 	  cout<< endl;
//      if(i==1)
//      {
//       OutPut(i << " Solution f_x " << Ddot(N_PBE_SpatialPts*N_LDof, PBE_LSolT, PBE_LSolT)<<endl);
// 
// exit(0);
//      } 
//         for(ii=0;ii<N_LDof;ii++)
//           RhsPbe_Intl_Loc[ii] = Oldrhs_IntlLoc[ii];

// // ===============================================================================================
// //  test 2D_1D_ConstT_UnitSqr_IntlOnly.h
// // ===============================================================================================
//           L_Sol2Nodal(FeSpace_Intl, PBE_LSolT+i*N_LDof, N_PBE_IntlLevels, PBE_IntlPtValuesT, 1);
// 
//              for(j=0;j<N_PBE_IntlLevels;j++)
//               {
//                L = IntlPosL[j];
//                x = TDatabase::TimeDB->CURRENTTIME;
//                cout << " L " << setw(10)<<   L   
//                     << setw(14)<< "  Sol " <<  PBE_IntlPtValuesT[j]    
//                     << setw(20)<<"  Err " <<  (1. - L)*exp(-x*L*L) - PBE_IntlPtValuesT[j] << endl;
//                }
//            exit(0);
//  
// ===============================================================================================
          // update due to the concentration coupling with PBE
          if(N_IndepntScalarEqns)
	   {
            C_IntlPtValues[i] = C;
	    T_IntlPtValues[i] = T;
	   }
        }   //  for(i=0;i<N_

//         OutPut("Solution f_x_L " << Ddot(N_PBE_SpatialPts*N_LDof, PBE_LSolT, PBE_LSolT)<<endl);

// if(m==1)
// exit(0);


        // transpose the solution, since we need it on rhs of the X-direction Stage 2
        if(TDatabase::ParamDB->REACTOR_P5==0)
        {
          if(N_IndepntScalarEqns)
	   {
            GetSolFromNodalPtVales(ScalarFunctions_low[1], C_IntlPtValues, Scalar_Spaces[PBE_INDEX_LOW]);
            GetSolFromNodalPtVales(ScalarFunctions_low[0], T_IntlPtValues, Scalar_Spaces[PBE_INDEX_LOW]);
	   }

          L_Sol2Nodal(FeSpace_Intl, PBE_LSolT, N_PBE_IntlLevels, PBE_IntlPtValuesT, N_PBE_SpatialPts);

          GetSolFromNodalPtVales(N_PBE_IntlLevels, MaxN_PtsForNodal, ScalarPbeFunctions_low, SolPbe_low, N_U_low, PBE_IntlPtValuesT, N_PBE_SpatialPts);
        }
        else
        {
          if(N_IndepntScalarEqns)
	   {
            GetSolFromQuadPtVales(ScalarFunctions_low[1], C_IntlPtValues, Scalar_Spaces[PBE_INDEX_LOW]);
	   }

          GetSolFromQuadPtVales(N_PBE_IntlLevels, MaxN_PtsForNodal, ScalarPbeFunctions_low, SolPbe_low, N_U_low, PBE_IntlPtValuesT, N_PBE_SpatialPts);
        }

     /** Prolongate sol from low fespace to high fespace for next timestep*/
     if(PBE_ORDER != ORDER)
     {
      Prolongate(Scalar_Spaces[PBE_INDEX_LOW], Scalar_Spaces[PBE_INDEX], N_PBE_IntlLevels,
               SolPbe_low, SolPbe, PbeAux); 
#ifndef __PBSConstT__        
      Prolongate(Scalar_Spaces_low[0], Scalar_Spaces[0],  SolArray_low[0], SolArray[0], PbeAux); 

      Prolongate(Scalar_Spaces_low[1], Scalar_Spaces[1],  SolArray_low[1], SolArray[1], PbeAux);  
#endif
     }

//        if( TDatabase::ParamDB->SC_VERBOSE > 0 )
//         {
//           t1 = GetTime();
//           OutPut("Time for solving Configuration space system : " <<  t1- t2 << "s"<<endl;);
//         }

//         // no change in L if it is DIRICHLET BC, correction due to interpolation/Transpose
// //        BD is setted in internal level so no update here
//         if(cond_Lmin==DIRICHLET && !dGDisc)
//           memcpy(SolPbe, SolPbe_Lmin, N_U*SizeOfDouble);
// 
//         if(cond_Lmax==DIRICHLET && !dGDisc)
//           memcpy(SolPbe+(N_U*(N_PBE_IntlLevels-1)), SolPbe_Lmax, N_U*SizeOfDouble);

     } //IF PBE_P0!=0
// cout<< " L_Sol2Nodal C_Sat " << C_Sat << endl;


     pbs_t2 = GetTime(); 
     pbs_step2 += (pbs_t2 - pbs_t1) - (pbs_t4 - pbs_t3);

     if(m%100==0 || m==1)
      {
       t2 = GetTime();
       OutPut( "Time for direct solving PBE_L: " << t2-t1 << endl);    
      }
//         cout << " N_PBE_SpatialPts " << N_PBE_SpatialPts << endl;
//  
// exit(0);

      //===============================================================================
      //PBE system L-direction solution --end
      //PBE system X-direction solution --start
      //===============================================================================

      //  //===============================================================================
      //  // test the errors after L-direction -end
      //  //===============================================================================
      //
      //   #ifdef _MPI
      //   if(ActiveProcess)
      //   #endif
      //    GetOSError(N_PBE_IntlLevels, ScalarPbeFunctions,  FeSpace_Intl, SolPbe, N_U, errors, Exact_Psd_Intl);
      //
      //
      //   #ifdef _MPI
      //   MPI_Reduce(errors, &l2, 1, MPI_DOUBLE, MPI_SUM, out_rank, MPI_COMM_WORLD);
      //   MPI_Reduce(errors+1, &H1, 1, MPI_DOUBLE, MPI_SUM, out_rank, MPI_COMM_WORLD);
      //   #else
      //   l2 = errors[0];
      //   H1 = errors[1];
      //   #endif // ifdef _MPI
      //
      //
      //   #ifdef _MPI
      //   if(rank==out_rank)
      //   #endif
      //    {
      //     olderror = l2 = sqrt(l2);
      //     olderror1 =  H1 = sqrt(H1);
      //     OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
      //     OutPut(" L2: " << l2);
      //     OutPut(" H1-semi: " << H1 << endl);
      //    }
      //  //===============================================================================
      //  // test the errors after L-direction -end
      //  //===============================================================================

#endif  //  #ifdef __PBS__

    }  //for(l=0;l<N_SubSteps;l++)


//         OutPut("Solution f_x_L " << Ddot(N_PBE_SpatialPts*N_LDof, SolPbe, SolPbe)<<endl);
// exit(0);

 
 
 
 
    /*****************************************************************************/
    /* postprocessing of current time step                                       */
    /*****************************************************************************/
    // measure errors to known solution
    if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
      for(i=0;i<N_IndepntScalarEqns;i++)
      {
        //  ScalarFunctions[i]->Interpolate(Exact);
          aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 1, NULL);
          ScalarFunctions[i]->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors,
            Coefficients[i], aux, 1, fesp, errors);
          delete aux;

        OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
        OutPut(" L2: " << errors[0]);
        OutPut(" H1-semi: " << errors[1] << endl);

        errors[3] += (errors[0]*errors[0] +olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
        olderror = errors[0];
        OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(errors[3]) << " ");

        errors[4] += (errors[1]*errors[1] +olderror1 * olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
        OutPut( "L2(0,T;H1) " << sqrt(errors[4]) << endl);
        olderror1 = errors[1];
      }

#ifdef __PBSConstT__
//       // interpolate exact solution - begin
//        memset(PBE_IntlPtValuesT, 0, N_PBE_SpatialPts*N_PBE_IntlLevels*SizeOfDouble);
//         for(i=0;i<N_PBE_SpatialPts;i++)
//           {
//            BDUpdatedSolPbe_Intl_Loc = PBE_IntlPtValuesT + i*N_PBE_IntlLevels;
//            ADI_System[i]->Interpolate(BDUpdatedSolPbe_Intl_Loc, Exact_Psd_Intl);
//  
//           }
//          if(TDatabase::ParamDB->REACTOR_P5==0)
//           {
//            GetSolFromNodalPtVales(N_PBE_IntlLevels, MaxN_PtsForNodal, ScalarPbeFunctions, SolPbe, N_U, PBE_IntlPtValuesT, N_PBE_SpatialPts);
// 
//           }
//          else
//           {
//            GetSolFromQuadPtVales(N_PBE_IntlLevels, MaxN_PtsForNodal, ScalarPbeFunctions, SolPbe, N_U, PBE_IntlPtValuesT, N_PBE_SpatialPts);
// 
//           }
//        // interpolate exact solution - end

     GetOSError(N_PBE_IntlLevels, ScalarPbeFunctions,  FeSpace_Intl, SolPbe, N_U, errors, Exact_Psd_Intl);

     l2 = errors[0];
     H1 = errors[1];

      //small l2
        errors[5] += (TDatabase::TimeDB->TIMESTEPLENGTH)*l2;

//         l2 = sqrt(l2);
//         H1 = sqrt(H1);
        OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
        OutPut(" L2: " << sqrt(l2));
        OutPut(" H1-semi: " << sqrt(H1) << endl);



        if(L2error_Max< sqrt(l2))
        {
          L2error_Max= sqrt(l2);
          L2error_Max_t=TDatabase::TimeDB->CURRENTTIME;
        }

        errors[3] += (l2*l2 +olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
        olderror = l2;
        OutPut(TDatabase::TimeDB->CURRENTTIME <<  " l2(0,T;L2) " << sqrt(errors[5]) << " ");
  // ifdef __PBSConstT__
        OutPut(" L2(0,T;L2) " << sqrt(errors[3]) << " ");

        OutPut(L2error_Max_t <<  " L2error_Max " << L2error_Max << " ");

        errors[4] += (H1*H1 +olderror1 * olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
        OutPut( "L2(0,T;H1) " << sqrt(errors[4]) << endl);
        olderror1 = H1;

#endif  // ifdef __PBSConstT__
    }
    //  printf( "Rank : %d,   %e\n", rank,  errors[0]  );// endif MEASURE_ERRORS
    // printf("Rank: %d time step %f  \n", rank, tau1);
    // MPI_Finalize();
//     exit(0);


#ifndef __PBSConstT__
#ifdef __PBS__
   if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0 || m==1)
    {   

     OutPut( "Time, Total time taken for OS_Step1, aggregation,   OS_Step2, N_Spatial Pts, N_TimeSteps: " << TDatabase::TimeDB->CURRENTTIME << "  " << pbs_step1 << "  " << aggr_time <<"  " << pbs_step2 << " " << N_PBE_SpatialPts << " " << m << endl);
     
     t4 =  GetTime();
     total_time += t4 - t3;
     OutPut("Time, Total running time, N_TimeSteps: " << TDatabase::TimeDB->CURRENTTIME << "  " << total_time <<  " " << m << endl);
     t3 =  GetTime();
     
     // calculating the psd mass
     for(i=0; i<N_PBE_SpatialPts; i++)
      {
       BDUpdatedSolPbe_Intl_Loc = PBE_IntlPtValuesT+(i*N_LDof);
       q3_max=ADI_System[i]->GetQ3Max(BDUpdatedSolPbe_Intl_Loc);
       pbeq3massAtNodalPts[i] = q3_max;
      }
      
     if(TDatabase::ParamDB->REACTOR_P5==0)
      {
       GetSolFromNodalPtVales(PSD_mass,  pbeq3massAtNodalPts, Scalar_Spaces[PBE_INDEX_LOW]);
      }
     else
      {
       GetSolFromQuadPtVales(PSD_mass, pbeq3massAtNodalPts, Scalar_Spaces[PBE_INDEX_LOW]);
      }
        
//      Get_FeFunction2DMass(PSD_mass, param);   
     PSD_mass->GetMassAndMean(param);
     OutPut( "Time, Vol, Mass of PSD_q3, T, C: " <<TDatabase::TimeDB->CURRENTTIME<< " " <<param[1]<<" " <<param[0]<< " "); 
     ScalarFunctions[0]->GetMassAndMean(param);
     OutPut(param[0]<< " ");  
     ScalarFunctions[1]->GetMassAndMean(param);
     OutPut(param[0]<< " "<<endl);     
     
     
     // WriteAvgQ3AtOutlet
     // first calculate the average for eqch i:  avg_L_i = \sum_j^N{ f(t_n, x_j, l_i) }/N
     memset(Oldrhs_IntlLoc, 0, N_LDof*SizeOfDouble);

     for(j=0; j<N_f_OutletIndex; j++)
      {
       BDUpdatedSolPbe_Intl_Loc = PBE_IntlPtValuesT+(f_OutletIndex[j]*N_LDof);
       
       for(i=0; i<N_LDof; i++)
        {
          Oldrhs_IntlLoc[i] += BDUpdatedSolPbe_Intl_Loc[i];    
         } //  for(i=0; i<N_LDof	     
      } //for(j=0; j<N_

      for(i=0; i<N_LDof; i++)
       { Oldrhs_IntlLoc[i] /=(double)N_f_OutletIndex; }
     
       // can be calculate using any pt, since cal. is done for give avg value of f
       q3_max=ADI_System[f_OutletIndex[0]]->GetQ3Max(Oldrhs_IntlLoc);
     
         os.seekp(std::ios::beg);
        if (N_Q3Data<10) { os << "Q3Data/OutletAvgQ3_0000"<<N_Q3Data<<".data" << ends; }
         else if(N_Q3Data<100) { os << "Q3Data/OutletAvgQ3_000"<<N_Q3Data<<".data" << ends; }
         else { os << "Q3Data/OutletAvgQ3_00"<<N_Q3Data<<".data" << ends; }
        std::ofstream dat(os.str().c_str());
        dat.setf(std::ios::fixed);
        dat << setprecision(12);

        if (!dat)
         {
          cerr << "cannot open file for output" << endl;
          return -1;
         }
        dat << "%% Outlet Avg. Q3 data created for droplet by ParMooN" << endl;
        dat << "%% Current Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
	dat << "%% x=200 line at outlet"<< endl;
        dat << "%% q3_max : " <<  setprecision(12) << q3_max << endl;
        dat << "%% L , Avg.f, Avg. q3: " << endl;	 



        if( fabs(q3_max)>1.0e-100)
         {
           // N_LDof should be equal to N_PBE_IntlLevels !!
          for(j=0; j<N_PBE_IntlLevels; j++)
           {
            L=IntlPosL[j];
            dat <<L <<  " " <<   setprecision(12) << Oldrhs_IntlLoc[j] << " " <<
                        setprecision(12) <<  L*L*L*Oldrhs_IntlLoc[j]/q3_max<<endl;
           }
         }
        else
         {
          for(j=0; j<N_PBE_IntlLevels; j++)
           {
            L=IntlPosL[j];
            dat <<L << " " <<  Oldrhs_IntlLoc[j] << " " << 0. <<endl;
           }
         }

        dat.close();
//         cout << endl;
        OutPut( "OutletAvgQ3 data wrote into file Q3Data/OutletAvgQ3*" << N_Q3Data << ".data" <<endl);
        N_Q3Data++;

   
//        // calculating the mass  
//       for(i=0;i<N_IndepntScalarEqns;i++)
//        {
//         OutPut("Solution " << i << " Ddot " << Ddot(N_Uarray[i], SolArray[i], SolArray[i])<<endl);
//        }

//        OutPut("Solution f_x_L " << Ddot(N_PBE_SpatialPts*N_LDof, SolPbe_low, SolPbe_low)<<endl);


       BDUpdatedSolPbe_Intl_Loc = PBE_IntlPtValuesT + f_Integ_NodalPt[0]*N_LDof;

       q3_max=ADI_System[f_Integ_NodalPt[0]]->GetQ3Max(BDUpdatedSolPbe_Intl_Loc);
//        OutPut( endl);
       OutPut(TDatabase::TimeDB->CURRENTTIME << " X " <<IntlX[f_Integ_NodalPt[0]]<< " Y " <<IntlY[f_Integ_NodalPt[0]] << " q3_max: " <<  setprecision(8) << q3_max << endl);
       
        if(m % 5000 == 0)
         {
          os.seekp(std::ios::beg);
          if (N_OutData<10) { os << "Q3Data/Outlet_X_"<<IntlX[f_Integ_NodalPt[0]]<< "_Y_" <<IntlY[f_Integ_NodalPt[0]]<< "_00"<<N_OutData<<".data" << ends; }
           else if(N_OutData<100) { os << "Q3Data/Outlet_X_"<<IntlX[f_Integ_NodalPt[0]]<< "_Y_" <<IntlY[f_Integ_NodalPt[0]]<< "_0"<<N_OutData<<".data" << ends; }
           else { os << "Q3Data/Outlet_X_"<<IntlX[f_Integ_NodalPt[0]]<< "_Y_" <<IntlY[f_Integ_NodalPt[0]]<< "_"<<N_OutData<<".data" << ends; }
           
          std::ofstream dat(os.str().c_str());
          dat.setf(std::ios::fixed);
          dat << setprecision(12);   
	  
         if (!dat)
          {
           cerr << "cannot open file for output" << endl;
           return -1;
          }
         dat << "%% Outlet Q3 data created by ParMooN" << endl;
         dat << "%% Current Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
	 dat << "%% x "<<IntlX[f_Integ_NodalPt[0]] << " Y " <<IntlY[f_Integ_NodalPt[0]]<<   endl;
         dat << "%% q3_max : " <<  setprecision(12) << q3_max << endl;
         dat << "%% L , f,  q3: " << endl;	 

        if( fabs(q3_max)>1.0e-100)
         {
           // N_LDof should be equal to N_PBE_IntlLevels !!
          for(j=0; j<N_PBE_IntlLevels; j++)
           {
            L=IntlPosL[j];
            dat <<L <<  " " <<   setprecision(12) << BDUpdatedSolPbe_Intl_Loc[j] << " " << setprecision(12) <<  L*L*L*BDUpdatedSolPbe_Intl_Loc[j]/q3_max<<endl;
           }
         }
        else
         {
          for(j=0; j<N_PBE_IntlLevels; j++)
           {
            L=IntlPosL[j];
            dat <<L << " " <<  BDUpdatedSolPbe_Intl_Loc[j] << " " << 0. <<endl;
           }
         }

        dat.close();
        OutPut( "Outlet point Q3 data wrote " <<endl);
    } // if(m % 5000 == 0)
       
 
       
       BDUpdatedSolPbe_Intl_Loc = PBE_IntlPtValuesT + f_Integ_NodalPt[1]*N_LDof; 
       q3_max=ADI_System[f_Integ_NodalPt[1]]->GetQ3Max(BDUpdatedSolPbe_Intl_Loc);
//        OutPut( endl);
       OutPut(TDatabase::TimeDB->CURRENTTIME << " X " <<IntlX[f_Integ_NodalPt[1]]<< " Y " <<IntlY[f_Integ_NodalPt[1]] << " q3_max: " <<  setprecision(8) << q3_max << endl);
        if(m % 5000 == 0)
         {
          os.seekp(std::ios::beg);
          if (N_OutData<10) { os << "Q3Data/Outlet_X_"<<IntlX[f_Integ_NodalPt[1]]<< "_Y_" <<IntlY[f_Integ_NodalPt[1]]<< "_00"<<N_OutData<<".data" << ends; }
           else if(N_OutData<100) { os << "Q3Data/Outlet_X_"<<IntlX[f_Integ_NodalPt[1]]<< "_Y_" <<IntlY[f_Integ_NodalPt[1]]<< "_0"<<N_OutData<<".data" << ends; }
           else { os << "Q3Data/Outlet_X_"<<IntlX[f_Integ_NodalPt[1]]<< "_Y_" <<IntlY[f_Integ_NodalPt[1]]<< "_"<<N_OutData<<".data" << ends; }
           
          std::ofstream dat(os.str().c_str());
          dat.setf(std::ios::fixed);
          dat << setprecision(12);   
	  
         if (!dat)
          {
           cerr << "cannot open file for output" << endl;
           return -1;
          }
         dat << "%% Outlet Q3 data created by ParMooN" << endl;
         dat << "%% Current Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
	 dat << "%% x "<<IntlX[f_Integ_NodalPt[1]] << " Y " <<IntlY[f_Integ_NodalPt[1]]<<   endl;
         dat << "%% q3_max : " <<  setprecision(12) << q3_max << endl;
         dat << "%% L , f,  q3: " << endl;	 

        if( fabs(q3_max)>1.0e-100)
         {
           // N_LDof should be equal to N_PBE_IntlLevels !!
          for(j=0; j<N_PBE_IntlLevels; j++)
           {
            L=IntlPosL[j];
            dat <<L <<  " " <<   setprecision(12) << BDUpdatedSolPbe_Intl_Loc[j] << " " << setprecision(12) <<  L*L*L*BDUpdatedSolPbe_Intl_Loc[j]/q3_max<<endl;
           }
         }
        else
         {
          for(j=0; j<N_PBE_IntlLevels; j++)
           {
            L=IntlPosL[j];
            dat <<L << " " <<  BDUpdatedSolPbe_Intl_Loc[j] << " " << 0. <<endl;
           }
         }

        dat.close();
        OutPut( "Outlet point Q3 data wrote " <<endl);
    } // if(m % 5000 == 0)       

       BDUpdatedSolPbe_Intl_Loc = PBE_IntlPtValuesT + f_Integ_NodalPt[2]*N_LDof;        
       q3_max=ADI_System[f_Integ_NodalPt[2]]->GetQ3Max(BDUpdatedSolPbe_Intl_Loc);
//        OutPut( endl);
       OutPut(TDatabase::TimeDB->CURRENTTIME << " X " <<IntlX[f_Integ_NodalPt[2]]<< " Y " <<IntlY[f_Integ_NodalPt[2]] << " q3_max: " <<  setprecision(8) << q3_max << endl);
        if(m % 5000 == 0)
         {
          os.seekp(std::ios::beg);
          if (N_OutData<10) { os << "Q3Data/Outlet_X_"<<IntlX[f_Integ_NodalPt[2]]<< "_Y_" <<IntlY[f_Integ_NodalPt[2]]<< "_00"<<N_OutData<<".data" << ends; }
           else if(N_OutData<100) { os << "Q3Data/Outlet_X_"<<IntlX[f_Integ_NodalPt[2]]<< "_Y_" <<IntlY[f_Integ_NodalPt[2]]<< "_0"<<N_OutData<<".data" << ends; }
           else { os << "Q3Data/Outlet_X_"<<IntlX[f_Integ_NodalPt[2]]<< "_Y_" <<IntlY[f_Integ_NodalPt[2]]<< "_"<<N_OutData<<".data" << ends; }
           
          std::ofstream dat(os.str().c_str());
          dat.setf(std::ios::fixed);
          dat << setprecision(12);   
	  
         if (!dat)
          {
           cerr << "cannot open file for output" << endl;
           return -1;
          }
         dat << "%% Outlet Q3 data created by ParMooN" << endl;
         dat << "%% Current Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
	 dat << "%% x "<<IntlX[f_Integ_NodalPt[2]] << " Y " <<IntlY[f_Integ_NodalPt[2]]<<   endl;
         dat << "%% q3_max : " <<  setprecision(12) << q3_max << endl;
         dat << "%% L , f,  q3: " << endl;	 

        if( fabs(q3_max)>1.0e-100)
         {
           // N_LDof should be equal to N_PBE_IntlLevels !!
          for(j=0; j<N_PBE_IntlLevels; j++)
           {
            L=IntlPosL[j];
            dat <<L <<  " " <<   setprecision(12) << BDUpdatedSolPbe_Intl_Loc[j] << " " << setprecision(12) <<  L*L*L*BDUpdatedSolPbe_Intl_Loc[j]/q3_max<<endl;
           }
         }
        else
         {
          for(j=0; j<N_PBE_IntlLevels; j++)
           {
            L=IntlPosL[j];
            dat <<L << " " <<  BDUpdatedSolPbe_Intl_Loc[j] << " " << 0. <<endl;
           }
         }

        dat.close();
        OutPut( "Outlet point Q3 data wrote " <<endl);
        N_OutData++;	  
     } // if(m % 5000 == 0)
    
    
    
//  ==============================================================================   
    
       BDUpdatedSolPbe_Intl_Loc = PBE_IntlPtValuesT + f_Integ_NodalPt[3]*N_LDof;        
       q3_max=ADI_System[f_Integ_NodalPt[3]]->GetQ3Max(BDUpdatedSolPbe_Intl_Loc);
//        OutPut( endl);
       OutPut(TDatabase::TimeDB->CURRENTTIME << " X " <<IntlX[f_Integ_NodalPt[3]]<< " Y " <<IntlY[f_Integ_NodalPt[3]] << " q3_max: " <<  setprecision(8) << q3_max << endl);
        if(m % 5000 == 0)
         {
          os.seekp(std::ios::beg);
          if (N_OutData<10) { os << "Q3Data/Outlet_X_"<<IntlX[f_Integ_NodalPt[3]]<< "_Y_" <<IntlY[f_Integ_NodalPt[3]]<< "_00"<<N_OutData<<".data" << ends; }
           else if(N_OutData<100) { os << "Q3Data/Outlet_X_"<<IntlX[f_Integ_NodalPt[3]]<< "_Y_" <<IntlY[f_Integ_NodalPt[3]]<< "_0"<<N_OutData<<".data" << ends; }
           else { os << "Q3Data/Outlet_X_"<<IntlX[f_Integ_NodalPt[3]]<< "_Y_" <<IntlY[f_Integ_NodalPt[3]]<< "_"<<N_OutData<<".data" << ends; }
           
          std::ofstream dat(os.str().c_str());
          dat.setf(std::ios::fixed);
          dat << setprecision(12);   
	  
         if (!dat)
          {
           cerr << "cannot open file for output" << endl;
           return -1;
          }
         dat << "%% Outlet Q3 data created by ParMooN" << endl;
         dat << "%% Current Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
	 dat << "%% x "<<IntlX[f_Integ_NodalPt[3]] << " Y " <<IntlY[f_Integ_NodalPt[3]]<<   endl;
         dat << "%% q3_max : " <<  setprecision(12) << q3_max << endl;
         dat << "%% L , f,  q3: " << endl;	 

        if( fabs(q3_max)>1.0e-100)
         {
           // N_LDof should be equal to N_PBE_IntlLevels !!
          for(j=0; j<N_PBE_IntlLevels; j++)
           {
            L=IntlPosL[j];
            dat <<L <<  " " <<   setprecision(12) << BDUpdatedSolPbe_Intl_Loc[j] << " " << setprecision(12) <<  L*L*L*BDUpdatedSolPbe_Intl_Loc[j]/q3_max<<endl;
           }
         }
        else
         {
          for(j=0; j<N_PBE_IntlLevels; j++)
           {
            L=IntlPosL[j];
            dat <<L << " " <<  BDUpdatedSolPbe_Intl_Loc[j] << " " << 0. <<endl;
           }
         }

        dat.close();
        OutPut( "Outlet point Q3 data wrote " <<endl);
        N_OutData++;	  
     } // if(m % 5000 == 0)
        
    
//  ==============================================================================   
    
       BDUpdatedSolPbe_Intl_Loc = PBE_IntlPtValuesT + f_Integ_NodalPt[4]*N_LDof;        
       q3_max=ADI_System[f_Integ_NodalPt[4]]->GetQ3Max(BDUpdatedSolPbe_Intl_Loc);
//        OutPut( endl);
       OutPut(TDatabase::TimeDB->CURRENTTIME << " X " <<IntlX[f_Integ_NodalPt[4]]<< " Y " <<IntlY[f_Integ_NodalPt[4]] << " q3_max: " <<  setprecision(8) << q3_max << endl);
        if(m % 5000 == 0)
         {
          os.seekp(std::ios::beg);
          if (N_OutData<10) { os << "Q3Data/Outlet_X_"<<IntlX[f_Integ_NodalPt[4]]<< "_Y_" <<IntlY[f_Integ_NodalPt[4]]<< "_00"<<N_OutData<<".data" << ends; }
           else if(N_OutData<100) { os << "Q3Data/Outlet_X_"<<IntlX[f_Integ_NodalPt[4]]<< "_Y_" <<IntlY[f_Integ_NodalPt[4]]<< "_0"<<N_OutData<<".data" << ends; }
           else { os << "Q3Data/Outlet_X_"<<IntlX[f_Integ_NodalPt[4]]<< "_Y_" <<IntlY[f_Integ_NodalPt[4]]<< "_"<<N_OutData<<".data" << ends; }
           
          std::ofstream dat(os.str().c_str());
          dat.setf(std::ios::fixed);
          dat << setprecision(12);   
	  
         if (!dat)
          {
           cerr << "cannot open file for output" << endl;
           return -1;
          }
         dat << "%% Outlet Q3 data created by ParMooN" << endl;
         dat << "%% Current Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
	 dat << "%% x "<<IntlX[f_Integ_NodalPt[4]] << " Y " <<IntlY[f_Integ_NodalPt[4]]<<   endl;
         dat << "%% q3_max : " <<  setprecision(12) << q3_max << endl;
         dat << "%% L , f,  q3: " << endl;	 

        if( fabs(q3_max)>1.0e-100)
         {
           // N_LDof should be equal to N_PBE_IntlLevels !!
          for(j=0; j<N_PBE_IntlLevels; j++)
           {
            L=IntlPosL[j];
            dat <<L <<  " " <<   setprecision(12) << BDUpdatedSolPbe_Intl_Loc[j] << " " << setprecision(12) <<  L*L*L*BDUpdatedSolPbe_Intl_Loc[j]/q3_max<<endl;
           }
         }
        else
         {
          for(j=0; j<N_PBE_IntlLevels; j++)
           {
            L=IntlPosL[j];
            dat <<L << " " <<  BDUpdatedSolPbe_Intl_Loc[j] << " " << 0. <<endl;
           }
         }

        dat.close();
        OutPut( "Outlet point Q3 data wrote " <<endl);
        N_OutData++;	  
     } // if(m % 5000 == 0)
        
    
//  ==============================================================================       

  } //  if(m % 200 == 0)
    
//     exit(0);
    
//    if(m % 100 == 0)
//     {
//        if(TDatabase::ParamDB->REACTOR_P5==0)
//          {
// 
//           /** transfer sol from spatial to internal */
//           GetPtsValuesForIntl(N_PBE_IntlLevels, ScalarPbeFunctions_low, SolPbe_low, N_U_low, PBE_IntlPtValuesT, N_PBE_SpatialPts);
// 
//           /** evaluate nodal values from nodal points */
//           L_Nodal2Sol(FeSpace_Intl, PBE_IntlPtValuesT, N_PBE_SpatialPts, N_PBE_IntlLevels, PBE_LSolT);
// 	  
// 	  
//          }
//         else
//          {
// //           GetPtsValuesForIntl_Quad(N_PBE_IntlLevels, ScalarPbeFunctions,  PBE_IntlPtValuesT, N_PBE_SpatialPts);
// 	   	    cout << " GetPtsValuesForIntl_Quad not yet implemented for QuadPts" <<endl;
// 	    exit(0);
//          }


//        BDUpdatedSolPbe_Intl_Loc = PBE_LSolT + f_Integ_NodalPt*N_LDof;

//        q3_max=ADI_System[f_Integ_NodalPt]->GetQ3Max(BDUpdatedSolPbe_Intl_Loc);
//        OutPut( endl);
//        OutPut( "q3_max " <<  setprecision(8) << q3_max << endl);
// 
//         os.seekp(std::ios::beg);
//         if (N_Q3Data<10) os << "Q3Data_0000"<<N_Q3Data<<".data" << ends;
//          else os << "Q3Data_000"<<N_Q3Data<<".data" << ends;
//         std::ofstream dat(os.str().c_str());
//         dat.setf(std::ios::fixed);
//         dat << setprecision(8);
// 
//         if (!dat)
//          {
//           cerr << "cannot open file for output" << endl;
//           return -1;
//          }
//         dat << "%% Q3 data created for droplet by MooNMD" << endl;
//         dat << "%% Current Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
// 	dat << "%% x=200 and y=0.5"<< endl;
//         dat << "%% q3_max : " <<  setprecision(8) << q3_max << endl;
//         dat << "%% L , f, q3/q3_max: " << endl;	 
// 
//        BDUpdatedSolPbe_Intl_Loc = PBE_IntlPtValuesT + f_Integ_NodalPt*N_PBE_IntlLevels;
// 
//         if( fabs(q3_max)>1.0e-100)
//          {
//           for(j=0; j<N_PBE_IntlLevels; j++)
//            {
//             L=IntlPosL[j];
//             dat <<L <<  " " <<   setprecision(8) << BDUpdatedSolPbe_Intl_Loc[j] << " " <<
//                         setprecision(8) <<  L*L*L*BDUpdatedSolPbe_Intl_Loc[j]/q3_max<<endl;
//            }
//          }
//         else
//          {
//           for(j=0; j<N_PBE_IntlLevels; j++)
//            {
//             L=IntlPosL[j];
//             dat <<L << " " <<  BDUpdatedSolPbe_Intl_Loc[j] << " " << L*L*L*BDUpdatedSolPbe_Intl_Loc[j] <<endl;
//            }
//          }
// 
//         dat.close();
//         cout << endl;
//         OutPut( "Q3 data wrote into file " << N_Q3Data <<endl);
//         N_Q3Data++;
// 
//         // int at the level (j==6 ||j==13 || j== 18 )
//        CompOutletAvgFunct(SolPbe, N_U, PbeFunctions, Out_CellIndex, Out_EdgeIndex, N_OutEdges);

//      } //if(m % TDatabase::TimeDB->10 == 0)
#endif
#endif // #ifndef __PBSConstT__

    if(m==1 || m % 200 == 0)
    {
      if(TDatabase::ParamDB->WRITE_VTK)
      {
        //      # ifdef _MPI
        //      if(N_IndepntScalarEqns)
        //       Output->Write_ParVTK(MPI_COMM_WORLD, img, PhySol);
        //
        //      #ifdef __PBS__
        //      if(N_PBEqns)
        //       {
        //        if(ActiveProcess)
        //         GetSolFromNodalPtVales(N_PBE_IntlLevels, MaxN_PtsForNodal, ScalarPbeFunctions, SolPbe, N_U, PBE_IntlPtValuesT, N_PBE_SpatialPts);
        //
        //        for(j=0; j<5; j++)
        //         {
        //          if(j)
        //          {  memcpy(SolPbe_Output, SolPbe+(j*Out_Level -1)*N_U, N_U*SizeOfDouble); }
        //          else
        //           { memcpy(SolPbe_Output, SolPbe+j*N_U, N_U*SizeOfDouble); }
        //
        //          strcpy (SubID, IntLevel);
        //          sprintf (RankBaseName, "%d", j);
        //          strcat(SubID, RankBaseName);
        //          Output_PBS->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
        //         } //for(j=0; j<5; j+
        //
        //        memcpy(SolPbe_Output, SolPbe+N_U, N_U*SizeOfDouble);
        //        strcpy (SubID, IntLevel);
        //        sprintf (RankBaseName, "%d", j);
        //        strcat(SubID, RankBaseName);
        //        Output_PBS->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
        // 
        //        } // if(N_PBEqns)
        //       #endif
        //      img++;
        //     #else

          if(N_IndepntScalarEqns)
          {
            os.seekp(std::ios::beg);
            if(img<10) os <<"VTK/"<< VtkBaseName<<".0000"<<img<<".vtk" << ends;
            else if(img<100) os <<"VTK/"<< VtkBaseName<<".000"<<img<<".vtk" << ends;
            else if(img<1000) os <<"VTK/"<< VtkBaseName<<".00"<<img<<".vtk" << ends;
            else if(img<10000) os <<"VTK/"<< VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<"VTK/"<< VtkBaseName<<"."<<img<<".vtk" << ends;
            Output->WriteVtk(os.str().c_str());
          }

#ifdef __PBS__
          for(j=0; j<N_PBE_IntlLevels; j++)
          {
//             TDatabase::ParamDB->REACTOR_P29=IntlPosL[j];
// //                       cout << j << " IntlPosL[j] " << IntlPosL[j] <<endl;
//             //          ScalarPbeFunctions[j]->Interpolate(Exact);
//             memcpy(SolPbe_Output, SolPbe+j*N_U, N_U*SizeOfDouble);
// 
//             os.seekp(std::ios::beg);
//             if(img<10) os <<"VTK/"<<  GnuBaseName<<j<<".0000"<<img<<".vtk" << ends;
//             else if(img<100) os <<"VTK/"<<  GnuBaseName<<j<<".000"<<img<<".vtk" << ends;
//             else if(img<1000) os <<"VTK/"<<  GnuBaseName<<j<<".00"<<img<<".vtk" << ends;
//             else if(img<10000) os <<"VTK/"<<  GnuBaseName<<j<<".0"<<img<<".vtk" << ends;
//             else  os <<"VTK/"<<  GnuBaseName<<j<<"."<<img<<".vtk" << ends;
//             Output_PBS->WriteVtk(os.str().c_str());

//             //testing
//             //          os.seekp(std::ios::beg);
//             //          os <<"VTK/"<<  GnuBaseName<<".0000"<<j+1<<".vtk" << ends;
//             //          else if(img<100) os <<"VTK/"<<  GnuBaseName<<".000"<<j+1<<".vtk" << ends;
//             //          else if(img<1000) os <<"VTK/"<<  GnuBaseName<<".00"<<j+1<<".vtk" << ends;
//             //          else if(img<10000) os <<"VTK/"<<  GnuBaseName<<".0"<<j+1<<".vtk" << ends;
//             //          else  os <<"VTK/"<<  GnuBaseName<<"."<<j+1<<".vtk" << ends;
//             //          Output_PBS->WriteVtk(os.str().c_str());
          }
#endif
          img++;
      }                                           // if(TDatabase::ParamDB->WRITE_VTK)

//         MPI_Finalize();
//         exit(0);

    }   
    



#ifdef __PBS__  

 



#endif    
// if(m % TDatabase::TimeDB->S
//            cout<< " test main time" <<endl;
//            exit(0);

//     if(TDatabase::TimeDB->CURRENTTIME>1500)
//       TDatabase::TimeDB->STEPS_PER_IMAGE = 10000;
//         OutPut("Solution f_x_L " << Ddot(N_PBE_SpatialPts*N_LDof, SolPbe, SolPbe)<<endl);
// exit(0);
  }                                               //while(TDatabase::TimeDB->CURRENTTIME< end_time)

  //    #ifdef _MPI
  //    for(i=0;i<N_IndepntScalarEqns;i++)
  //       ParSolVect[i]->AssembleAtRoot();
  //
  //    ParSolPbe->AssembleAtRoot();
  //
  //    if(!ActiveProcess)
  //      #endif
  //       {
  //
  //        if(N_IndepntScalarEqns)
  //         {
  //          os.seekp(std::ios::beg);
  //          if(img<10) os <<"VTK/"<< VtkBaseName<<".0000"<<img<<".vtk" << ends;
  //          else if(img<100) os <<"VTK/"<< VtkBaseName<<".000"<<img<<".vtk" << ends;
  //          else if(img<1000) os <<"VTK/"<< VtkBaseName<<".00"<<img<<".vtk" << ends;
  //          else if(img<10000) os <<"VTK/"<< VtkBaseName<<".0"<<img<<".vtk" << ends;
  //          else  os <<"VTK/"<< VtkBaseName<<"."<<img<<".vtk" << ends;
  //          Output->WriteVtk(os.str().c_str());
  //         }
  //
  //      #ifdef __PBS__
  //        for(j=0; j<N_PBE_IntlLevels; j++)
  //         {
  //          TDatabase::ParamDB->REACTOR_P29=IntlPosL[j];
  // //           cout << " IntlPosL[j] " << IntlPosL[j] <<endl;
  //          ScalarPbeFunctions[j]->Interpolate(Exact);
  //          memcpy(SolPbe_Output, SolPbe+j*N_U, N_U*SizeOfDouble);
  //
  //          os.seekp(std::ios::beg);
  //          if(img<10) os <<"VTK/"<<  GnuBaseName<<j<<".0000"<<img<<".vtk" << ends;
  //          else if(img<100) os <<"VTK/"<<  GnuBaseName<<j<<".000"<<img<<".vtk" << ends;
  //          else if(img<1000) os <<"VTK/"<<  GnuBaseName<<j<<".00"<<img<<".vtk" << ends;
  //          else if(img<10000) os <<"VTK/"<<  GnuBaseName<<j<<".0"<<img<<".vtk" << ends;
  //          else  os <<"VTK/"<<  GnuBaseName<<j<<"."<<img<<".vtk" << ends;
  //          Output_PBS->WriteVtk(os.str().c_str());
  //
  //
  //        }
  //       #endif
  //    img++;
  //     }

  t4 =  GetTime();
  total_time += t4 - t3;
  OutPut("Total running time: " << total_time << endl);
  CloseFiles();
  return 0;
}
