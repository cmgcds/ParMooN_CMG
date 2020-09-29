// =======================================================================
//
// Purpose:     main program ADI scheme
//               1D + 1D  
// 
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 06.11.2009
// 
// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <DirectSolver.h>
#include <SquareStructure1D.h>
#include <SquareMatrix1D.h>
#include <FEFunction1D.h>
#include <LineAffin.h>
#include <LinAlg.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include <MainUtilities.h>
#include <TimeUtilities.h>

#include <ADICell1D.h>

#include "../Examples/TCD_2D/Time1_ADI.h"
// #include "../Examples/TCD_2D/Gaussian-BenchMark_ADI.h"
// ======================================================================
// utilities for main program
// ======================================================================
void GetQuadPtsValuesForIntl(int N, TBaseCell *Cell, TFEFunction1D **u, int N_InternalLevel, double *G,
                             double *QuadPtsRhsT)
{
  int  i, j, k, l, N_Points, N_BaseFunct, N_Sets=1;
  int *GlobalNumbers, *BeginIndex, *DOF;

  double *Weights, *zeta, X[MaxN_QuadPoints_1D], AbsDetjk[MaxN_QuadPoints_1D];
  double **origvaluesD0, *sol, val, *value;
  double *orgD0;

  bool Needs2ndDer[1];

  TFESpace1D *FeSpace;
  FE1D FEId;
  TFE1D *Element;
  TBaseFunct1D *bf;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TRefTrans1D *F_K;
  BaseFunct1D BaseFunct_ID, BaseFunct[1];

  Needs2ndDer[0] = FALSE;

  // assume that all internal levels use same FE Space
  FeSpace = u[0]->GetFESpace1D();
  GlobalNumbers = FeSpace->GetGlobalNumbers();
  BeginIndex = FeSpace->GetBeginIndex();

  FEId = FeSpace->GetFE1D(N, Cell);
  Element = TFEDatabase2D::GetFE1D(FEId);
  N_BaseFunct = Element->GetN_DOF();
  BaseFunct_ID = Element->GetBaseFunct1D_ID();

  bf = Element->GetBaseFunct1D();
  l = bf->GetPolynomialDegree();
  LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
  qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
  qf1->GetFormulaData(N_Points, Weights, zeta);

  F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
  ((TLineAffin *)F_K)->SetCell(Cell);
  //((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, X, AbsDetjk);

  //OutPut("N_Points  : "<< setw(4) << N_Points << endl);

  BaseFunct[0] = BaseFunct_ID;
  ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta,  LineQuadFormula,  Needs2ndDer);

  origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
//   origvaluesD1=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D1);

  DOF = GlobalNumbers + BeginIndex[N];

  for(i=0; i<N_InternalLevel; i++)
   {
    sol = u[i]->GetValues();

    for(j=0; j<N_Points; j++)
     {
      orgD0 = origvaluesD0[j];
      val = 0.;
      for(k=0;k<N_BaseFunct;k++)
       {
        val += orgD0[k]*sol[ DOF[k] ]; // for M
       }
      QuadPtsRhsT[j*N_InternalLevel + i] = val;
//         if(max < val) max = val;
     } //  for(j=0; 
   } // for(i=0*/

//   BilinearCoeffs();
  for(i=0; i<N_Points; i++)
   {
    G[i] = 0.1;
   } //  for(i=0; 



}


void AssembleM(TFESpace1D *FeSpace, TSquareMatrix1D *M, BoundCondFunct2D *BoundaryCondition, BoundValueFunct2D *BoundValue)
{
 int i, j, k, l, N_Cells, N_BaseFunct, N_U;
 int N_Points, N_Sets=1, *GlobalNumbers, *BeginIndex, *DOF;
 int TestDOF, begin, end, *RowPtr, *KCol;

 double *Weights, *zeta, X[20], AbsDetjk[20];
 double LocMatrixM[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double **origvaluesD0, **origvaluesD1, Mult;
 double *orgD0, *orgD1, test0, test1, ansatz0, ansatz1, *ValuesM;
 double x=0.;

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
    //origvaluesD1=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D1);

    memset(LocMatrixM, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);

    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Points;j++)
     {
      Mult = Weights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];
      //orgD1 = origvaluesD1[j];

      for(k=0;k<N_BaseFunct;k++)
       {
        test0  = orgD0[k];
        //test1  = orgD1[k];

        // cout<< " uref " << origvaluesD1[j][k] << " orgD0 " << orgD1[k] <<endl;
        for(l=0;l<N_BaseFunct;l++)
         {
          ansatz0  = orgD0[l];
          //ansatz1  = orgD1[l];
         LocMatrixM[k*N_BaseFunct + l] += Mult*ansatz0*test0;
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
           //ValuesA[k] +=LocMatrixA[j*N_BaseFunct + l];
           ValuesM[k] +=LocMatrixM[j*N_BaseFunct + l];
           break;
          }
        } // for(m=0;m<N_BaseFunct_low
      } // for(n=begin;n<end;n++)
     } // for(l=0;l<N_BaseFunct_low
   }// for(i=0; i<N_Cells; i++)

 //update boundary data
 // starting point
  BoundaryCondition(0, x, BDType);

  if(BDType==DIRICHLET)
   {
    begin = RowPtr[0];
    end = RowPtr[1];

    for(k=begin;k<end;k++)
     {
      if(KCol[k] == 0 )
       ValuesM[k] = 1.;
      else
       ValuesM[k] = 0.;
    }
   }

 // end point
  BoundaryCondition(1, x, BDType);
  if(BDType==DIRICHLET)
   {
    begin = RowPtr[N_U-1];
    end = RowPtr[N_U];

    for(k=begin;k<end;k++)
     {
      if(KCol[k] == N_U-1 )
       ValuesM[k] = 1.;
      else
       ValuesM[k] = 0.;
     }
   }

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

// exit(0);
}


void AssembleMatARhs(TFESpace1D *FeSpace, TSquareMatrix1D *A, double *Rhs, double *OldSol, CoeffFct2D *BilinearCoeffs,
                     BoundCondFunct2D *BoundaryCondition, BoundValueFunct2D *BoundValue, double y,
                     int N_InternalLevel, int *N_QuadPts, int TotalN_QuadPoints, double *QuadPtsSol)
{
 int i, j, k, l, N_Cells, N_BaseFunct;
 int N_Points, N_Sets=1, *GlobalNumbers, *BeginIndex, *DOF;
 int TestDOF, begin, end, *RowPtr, *KCol, disp, N_U;

 double *Weights, *zeta, X[MaxN_QuadPoints_1D], AbsDetjk[MaxN_QuadPoints_1D];
 double **aux, **coeff, *Coeff;
 double Space_X[MaxN_QuadPoints_1D];
 double LocMatrixA[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double oldsol[MaxN_BaseFunctions1D];
 double **origvaluesD0, **origvaluesD1, Mult;
 double *orgD0, *orgD1, test0, test1, ansatz0, ansatz1, *ValuesA;
 double c0, c1, g0, rhsval, val, BDValue, *rhsvalues;
 bool Needs2ndDer[1];

 TCollection *Coll;
 TBaseCell *Cell;
 FE1D FEId;
 TFE1D *Element;
 TBaseFunct1D *bf;
 QuadFormula1D LineQuadFormula;
 TQuadFormula1D *qf1;
 TRefTrans1D *F_K;
 BaseFunct1D BaseFunct_ID, BaseFunct[1];
 BoundCond BDType;

  GlobalNumbers = FeSpace->GetGlobalNumbers();
  BeginIndex = FeSpace->GetBeginIndex();
  Coll = FeSpace->GetCollection();
  N_U = FeSpace->GetN_DegreesOfFreedom();

  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  ValuesA = A->GetEntries(); 

  N_Cells = Coll->GetN_Cells();
  Needs2ndDer[0] = FALSE;

  aux = new double *[MaxN_QuadPoints_1D];
  coeff = new double *[MaxN_QuadPoints_1D];

  for(i=0; i<MaxN_QuadPoints_1D; i++)
   {
    Space_X[i] = y;
    aux[i] = new double[6];
    coeff[i] = new double[6];
   }


  disp = 0;
  for(i=0; i<N_Cells; i++)
   {
    rhsvalues = QuadPtsSol+disp;

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

    if(N_QuadPts[i] != N_Points)
     {
      OutPut("assemble space N_QuadPts[i] != N_Points : "<<  N_QuadPts[i] <<" " << N_Points << endl);
      exit(0);
     }
    disp +=N_Points;

    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(Cell);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, X, AbsDetjk);

    BilinearCoeffs(N_Points, Space_X, X, aux, coeff);

    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta,  LineQuadFormula,  Needs2ndDer);

    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    origvaluesD1=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D1);

    memset(LocMatrixA, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(oldsol, 0, N_BaseFunct*SizeOfDouble);

    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Points;j++)
     {
      Coeff = coeff[j];
      c0 = Coeff[0];
      g0 = Coeff[1];
      c1 = Coeff[2];

      // internal coordinate value at this quad point
      rhsval = rhsvalues[j];

      Mult = Weights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];
      orgD1 = origvaluesD1[j];

      for(k=0;k<N_BaseFunct;k++)
       {
        test0  = orgD0[k];
        test1  = orgD1[k];

        oldsol[k] +=  rhsval*Mult*test0;

        for(l=0;l<N_BaseFunct;l++)
         {
          ansatz0  = orgD0[l];
          ansatz1  = orgD1[l];

          val  = c0*ansatz1*test1;
          val += g0*ansatz1*test0;
          val += c1*ansatz0*test0;

         // cout << c0<< " Mult*val " <<  Mult*val << endl;
          LocMatrixA[k*N_BaseFunct + l] += Mult*val;
        }
       }
     }

    // add to global matrices
    for(j=0;j<N_BaseFunct;j++)
     {
      TestDOF = DOF[j];
      OldSol[TestDOF] += oldsol[j]; // soluble surfactant relation

      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(k=begin;k<end;k++)
       {
       for(l=0;l<N_BaseFunct;l++)
        {
         if(KCol[k] == DOF[l])
          {
           ValuesA[k] +=LocMatrixA[j*N_BaseFunct + l];
           break;
          }
        } // for(m=0;m<N_BaseFunct_low
      } // for(n=begin;n<end;n++)
     } // for(l=0;l<N_BaseFunct_low

  }// for(i=0; i<N_Cells; i++)


//   memcpy(OldSol, Rhs, N_U*SizeOfDouble);

 //update boundary data
 // starting point
  BoundaryCondition(0, y, BDType);
  BoundValue(0, y, BDValue);

  if(BDType==DIRICHLET)
   {
    Rhs[0] = BDValue;

    begin = RowPtr[0];
    end = RowPtr[1];

    for(k=begin;k<end;k++)
     {
      if(KCol[k] == 0 )
       ValuesA[k] = 1.;
      else
       ValuesA[k] = 0.;
    }
   }
  else if(BDType== NEUMANN)
   {
    Rhs[0] += BDValue;
   }
  else
   {
    OutPut("Unknown internal boundary condition !"<< endl); 
    exit(0);
   }

 // end point
  BoundaryCondition(1, y, BDType);
  BoundValue(1, y, BDValue);

  if(BDType==DIRICHLET)
   {
    Rhs[N_U-1] = BDValue;

    begin = RowPtr[N_U-1];
    end = RowPtr[N_U];

    for(k=begin;k<end;k++)
     {
      if(KCol[k] == N_U-1 )
       ValuesA[k] = 1.;
      else
       ValuesA[k] = 0.;
     }
   }
  else if(BDType== NEUMANN)
   {
    Rhs[N_U-1] += BDValue;
   }
  else
   {
    OutPut("Unknown internal boundary condition !"<< endl); 
    exit(0);
   }


// //print matrix
//   for(j=0;j<N_U;j++)
//    {
//     begin = RowPtr[j];
//     end = RowPtr[j+1];
//     for(k=begin;k<end;k++)
//      {
//       cout << "A(" << j << ", "<< KCol[k] << ") = " << ValuesA[k] <<endl;
//      }
//     cout<<endl;
//    }
// exit(0);

  for(i=0; i<MaxN_QuadPoints_1D; i++)
   {
    delete [] aux[i];
    delete [] coeff[i];
   }

  delete [] aux;
  delete [] coeff;



}

// ======================================================================
// main program
// ======================================================================

int main(int argc, char* argv[])
{
  // ======================================================================
  // variable declaration
  // ======================================================================
  int i, j, k, l, m, ret, N_Cells, N_Cells_Internal;
  int N_U, N_V, N_Points, MaxN_QdPoints, *N_InternalLevels, *N_QuadPts;
  int min, max, N_InternalLevel, N_SubSteps, time_discs, very_first_time=0;
  int TotalN_QuadPoints, disp;

  double t1, t2, t3, total_time, x, y, gamma, oldtau, end_time;
  double **Sol, **Rhs, *Weights, *zeta, X[20], absdetjk[20];
  double *OldSol, minval, maxval;
  double *QuadPtsRhsT, *Internal_Xpos, *GridX, tau;
  double *QuadPtsSol, *QuadPtsSol_Loc, G[20], *values;
  double *B, *defect, BDValue, exactvalues[5];
  
  TDomain *Domain = new TDomain();
  TDomain *Domain_Internal = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TSquareStructure1D *Struct, *SqStruct1D_Intl;
  TSquareMatrix1D *M, *A, *SQMATRICES[2], *M_Intl, *A_Intl;
  TFEFunction1D **u, *GridFEFunction1D;
  BoundCond BDType;
  FE1D *FE1D_List, *FE1D_List_Intl;
  TCollection *Coll, *Coll_Intl;
  TBaseCell *Cell;
  TFESpace1D *FeSpace, *FeSpace_Intl, *FeSpaces[2];
  FE1D FEId;
  TFE1D *FE;
  TBaseFunct1D *bf;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  TRefTrans1D *F_K;
  TADICell1D **ADIColl;
  std::ostringstream os;
  os << " ";

  char *PRM, *GEO;
  char *PsBaseName, *GrapeBaseName, *GnuBaseName, *ReadGrapeBaseName;
  char *VtkBaseName;
  char ReadinDat[] = "readin.dat";
  char UString[] = "U";
  char GString[] = "Grid";
  char VString[] = "V";
//======================================================================
// read parameter file
//======================================================================
  total_time = GetTime();
  if(argc>=2)
    ret=Domain->ReadParam(argv[1]);
  else
    ret=Domain->ReadParam(ReadinDat);


  if(ret==-1)
  {
    exit(-1);
  }

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();
  i = int(TDatabase::ParamDB->REACTOR_P10);
  Generate1DUniformMesh(Domain, 0., 1., i);

  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;

  t3 = GetTime();
  total_time = t3 - total_time;
// *****************************************************************************
// FeSpace and memory for all matrices
// assume that the convection and the reaction coefficient terms are independent 
// of internal coordinates, so lhs matrices are same for all lelvels of internal coordinate
// *****************************************************************************
  Coll = Domain->GetCollection(It_Finest, 0);
  N_Cells= Coll->GetN_Cells();
  cout<< "N_Cells " << N_Cells <<endl;

  ADIColl = new TADICell1D*[N_Cells];

  FE1D_List = new FE1D[N_Cells];
  for(j=0;j<N_Cells;j++)
   FE1D_List[j] = FE1D(TDatabase::ParamDB->ANSATZ_ORDER);

  FeSpace = new TFESpace1D(Coll , UString, UString, FE1D_List);
  N_U = FeSpace->GetN_DegreesOfFreedom();

  Struct = new TSquareStructure1D(FeSpace);
  Struct->Sort();

  M = new TSquareMatrix1D(Struct);
  A = new TSquareMatrix1D(Struct);

  OutPut("dof space : "<< setw(10) << N_U << endl);

// *****************************************************************************
// Construct FeSpace and all data for the internal domain
// FESpace is same in internal coordinate for all  QuadPts
// *****************************************************************************
  i = int(TDatabase::ParamDB->REACTOR_P11);
  Generate1DUniformMesh(Domain_Internal, 0., 1., i);
  Coll_Intl = Domain_Internal->GetCollection(It_Finest, 0);
  N_Cells_Internal= Coll_Intl->GetN_Cells();
  cout<< "N_Cells_Internal " << N_Cells_Internal <<endl;

  FE1D_List_Intl = new FE1D[N_Cells_Internal];
  for(i=0;i<N_Cells_Internal;i++)
   FE1D_List_Intl[i] = FE1D(TDatabase::ParamDB->TEST_ORDER);

  FeSpace_Intl = new TFESpace1D(Coll_Intl, VString, VString, FE1D_List_Intl);

  OutPut("dof space : "<< setw(10) << FeSpace_Intl->GetN_DegreesOfFreedom() << endl);


  SqStruct1D_Intl = new TSquareStructure1D(FeSpace_Intl);
  SqStruct1D_Intl->Sort();
  M_Intl = new TSquareMatrix1D(SqStruct1D_Intl);
  A_Intl = new TSquareMatrix1D(SqStruct1D_Intl);

  N_InternalLevels = new int [N_Cells];
  N_QuadPts = new int [N_Cells];

  MaxN_QdPoints = 0;
  TotalN_QuadPoints = 0;
  for(i=0; i<N_Cells; i++)
   {
    Cell = Coll->GetCell(i);
    FEId = FeSpace->GetFE1D(i, Cell);
    FE = TFEDatabase2D::GetFE1D(FEId);
    bf = FE->GetBaseFunct1D();
    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_Points, Weights, zeta);

    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(Cell);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, X, absdetjk);

    ADIColl[i] = new TADICell1D(Cell, FeSpace_Intl, M_Intl, A_Intl, N_Points, X,
                                Initial, Exact);

    // assume that the internal domain is uniformaly discretised in a cell
    // i.e., N_InternalLevels is same for all QuadPt in this cells
    N_V = ADIColl[i]->GetN_InternalLevels();
    N_InternalLevels[i] = N_V;
    N_QuadPts[i] = N_Points;
    TotalN_QuadPoints += N_Points;
    if(MaxN_QdPoints<N_Points) MaxN_QdPoints=N_Points;
   }

// *****************************************************************************

   TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

   min = 1e8;
   max = -100;
   for(i=0; i<N_Cells; i++)
    {
     if(min>N_InternalLevels[i]) min =  N_InternalLevels[i];
     if(max<N_InternalLevels[i]) max = N_InternalLevels[i];
    }

  if(min!=max)
   {
    OutPut(" unstructured internal discretisation not yet implemented "<<  min << "  "  << max << endl);
    exit(0);
   }
  N_InternalLevel = max;
  OutPut("MaxN_QdPoints : "<< setw(4) << MaxN_QdPoints << endl);
  OutPut("N_InternalLevel : "<< setw(4) << N_InternalLevel << endl);

  QuadPtsRhsT = new double[MaxN_QdPoints*N_InternalLevel];
  QuadPtsSol_Loc = new double[MaxN_QdPoints*N_InternalLevel];
  QuadPtsSol = new double[N_Cells*MaxN_QdPoints*N_InternalLevel];

  Sol = new double*[N_InternalLevel];
  Rhs = new double*[N_InternalLevel];
  u = new TFEFunction1D*[N_InternalLevel];


  // assume that all grids in internal domain has uniform structure
  //
  Internal_Xpos = ADIColl[0]->Get_Xpos();

  for(i=0; i<N_InternalLevel; i++)
   {
//     cout<< i << " y " << Internal_Xpos[i] << endl;
    Sol[i] = new double[N_U];
    Rhs[i] = new double[N_U];

    memset(Sol[i], 0, N_U*SizeOfDouble);
    memset(Rhs[i], 0, N_U*SizeOfDouble);

    u[i] = new TFEFunction1D(FeSpace, UString, UString, Sol[i], N_U);
    y = Internal_Xpos[i];
    u[i]->Interpolate(1, y, Initial);
   }


  B = new double[N_U];
  defect = new double[N_U];
  OldSol = new double[N_U];
  GridX = new double[N_U];
  GridFEFunction1D = new TFEFunction1D(FeSpace, GString, GString, GridX, N_U);
  GridFEFunction1D->GridToData();

/*
  for(i=0; i<N_InternalLevel; i++)
   for(j=0; j<N_U; j++)
    cout<<GridX[j] << " " << Internal_Xpos[i] <<  " " << Sol[i][j] <<endl;
*/
 

  // parameters for time stepping scheme
  gamma = 0;
  m = 0;
  N_SubSteps = GetN_SubSteps();
  oldtau = 1;
  end_time = TDatabase::TimeDB->ENDTIME;

  // not active : TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL = 0
  if (TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL)
    time_discs = 2;
  else
    time_discs = 1;

  // Assemble the mass matrices(Space & Intl) once, it will not change in time loop
  M->Reset();
  AssembleM(FeSpace, M, BoundCondition_Space, BoundValue_Space);
  M_Intl->Reset();
  AssembleM(FeSpace_Intl, M_Intl, BoundCondition_Internal, BoundValue_Internal);

  maxval = -1.e8;
  
  
  
  
  
  //======================================================================
  // start of time cycle
  //======================================================================
  while(TDatabase::TimeDB->CURRENTTIME< end_time)
   {                                               // time cycle
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

    for(l=0;l<N_SubSteps;l++)                     // sub steps of fractional step theta
     {
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


      OutPut(endl << "CURRENT TIME: ");
      OutPut(TDatabase::TimeDB->CURRENTTIME << endl);
        // OutPut("memory: " << setw(10) << GetMemory() << endl);

     //********************************************************************************************
     // solve the system in internal direction --- begin
     // so that all values will be available in all internal levels
     // then solve in all internal level
     //********************************************************************************************
     disp = 0;
     for(i=0; i<N_Cells; i++)
      {
       Cell = Coll->GetCell(i);
       // find the convection term G(X) in all Quad ptsd for the internal direction
       // assumed that the internal convection is independent of internal coordinates
       // i.e, internal convection is constant for internal direction
       // so internal convection at internal level "0" is used for all levels
       GetQuadPtsValuesForIntl(i, Cell, u, N_InternalLevel, G, QuadPtsRhsT);

       // mass matrix is already assembled, restore it at end since needed for next cell
       // only assemble stiffness Mat and solve
       ADIColl[i]->SolveAllQdPts(G, QuadPtsRhsT, BilinearCoeffs_Internal, 
                                 BoundCondition_Internal, BoundValue_Internal,
                                 tau, QuadPtsSol_Loc);

      // copy internal sol in each cell to all levels
      for(j=0; j<N_QuadPts[i]; j++)
       {
        values = QuadPtsSol_Loc + j*N_InternalLevel;
//         values = QuadPtsRhsT + j*N_InternalLevel; // testing
        for(k=0; k<N_InternalLevel; k++)
         {
          QuadPtsSol[k*TotalN_QuadPoints + disp + j ] = values[k];
         }
       }
      disp +=N_QuadPts[i];
     }

     //********************************************************************************************
     // solve the system in internal direction --- end
     // solve the system in space direction for all internal levels --- begin
     //********************************************************************************************

//         OutPut("Writing data"<< TDatabase::TimeDB->CURRENTTIME<<".data"<< endl);
//         os.seekp(std::ios::beg);
//         os << "data"<< TDatabase::TimeDB->CURRENTTIME<<".data" << ends;
//         std::ofstream dat(os.str().c_str());
// 
//        if (!dat)
//         {
//          cerr << "cannot open file for output" << endl;
//          exit(0);
//         }
// 
//         if (!dat)
//         {
// 	  cerr << "cannot open file for output" << endl;
// 	  exit(0);
//         }



     for(i=0; i<N_InternalLevel; i++)
      {
       // working array for rhs is B, initialize B
       memset(B, 0, N_U*SizeOfDouble);
//        Daxpy(N_U, tau*TDatabase::TimeDB->THETA3, Rhs[i], B);

       //assemble the stiffness matrix and rhs 
       A->Reset();
       memset(Rhs[i], 0, N_U*SizeOfDouble);
       memset(OldSol, 0, N_U*SizeOfDouble);
       AssembleMatARhs(FeSpace, A, Rhs[i], OldSol, BilinearCoeffs_Space, BoundCondition_Space, BoundValue_Space, Internal_Xpos[i],
                       N_InternalLevel, N_QuadPts, TotalN_QuadPoints, QuadPtsSol+i*TotalN_QuadPoints);

       //memcpy(Sol[i], OldSol, N_U*SizeOfDouble);

//        Daxpy(N_U, tau*TDatabase::TimeDB->THETA4, Rhs[i], B);

       // M = M -(tau*TDatabase::TimeDB->THETA2) A
//        MatAdd(M, A, -tau*TDatabase::TimeDB->THETA2);
       // set current factor of steady state matrix
//        gamma = -tau*TDatabase::TimeDB->THETA2;

       // defect = M * sol
       // B:= B + defect
//        memset(defect, 0, N_U*SizeOfDouble);
//        MatVectActive(M, OldSol, defect);
//        Daxpy(N_U, 1, defect, B);
       Daxpy(N_U, 1, OldSol, B);
       //for(j=0; j<N_U; j++)
        // cout<< OldSol[j] << " defect " << i << " " <<  j<<  " " << defect[j] <<endl;
       //***************************************************************************************
       // set Dirichlet values
       // starting point
       BoundCondition_Space(0, Internal_Xpos[i], BDType);
       BoundValue_Space(0, Internal_Xpos[i], BDValue);

       if(BDType==DIRICHLET)
        {
         B[0] = BDValue;
         Sol[i][0] = BDValue;
        }

       //endpoint
       BoundCondition_Space(1, Internal_Xpos[i], BDType);
       BoundValue_Space(1, Internal_Xpos[i], BDValue);

       if(BDType==DIRICHLET)
        {
         B[N_U-1] = BDValue;
         Sol[i][N_U-1] = BDValue;
          }
       //***************************************************************************************

       // system matrix
       MatAdd(M, A, -gamma + tau*TDatabase::TimeDB->THETA1);
       // set current factor of steady state matrix
       gamma = tau*TDatabase::TimeDB->THETA1;

      DirectSolver(M, B, Sol[i]);

       MatAdd(M, A, -gamma);
       // set current factor of steady state matrix
       gamma = 0;

//       u[i]->Interpolate(1, Internal_Xpos[i], Initial);
//      for(j=0; j<N_U; j++)
//       {
//       dat <<GridX[j] << " " << Internal_Xpos[i] <<  " " << Sol[i][j] <<endl;
//        dat <<GridX[j] << " " << Internal_Xpos[i] <<  " " << OldSol[j] <<endl;

//               if(maxval < OldSol[j]) maxval = OldSol[j];
//       }


      } // for(i=0; i<N_InternalLevel; i++)

//        dat.close();
// 
// cout << " max main value " << maxval <<endl;
//        cout << "Min Value " << minval  << " Max Value " << maxval << endl;
// 
// exit(0);


     if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      {
        OutPut("Writing data"<< TDatabase::TimeDB->CURRENTTIME<<".data"<< endl);
        os.seekp(std::ios::beg);
        os << "data"<< TDatabase::TimeDB->CURRENTTIME<<".data" << ends;
        std::ofstream dat(os.str().c_str());

       if (!dat)
        {
         cerr << "cannot open file for output" << endl;
         exit(0);
        }

        if (!dat)
        {
	  cerr << "cannot open file for output" << endl;
	  exit(0);
        }

       minval= 1.e8;
       maxval= -1.e8;
       for(i=0; i<N_InternalLevel; i++)
        {
         for(j=0; j<N_U; j++)
          {
//           Exact(GridX[j], Internal_Xpos[i], exactvalues);
          dat <<GridX[j] << " " << Internal_Xpos[i] <<  " " <<  Sol[i][j] <<endl;
//           dat <<GridX[j] << " " << Internal_Xpos[i] <<  " " << Sol[i][j] - exactvalues[0] <<endl;
          if( minval>  Sol[i][j] ) minval=Sol[i][j];
           if( maxval<  Sol[i][j] ) maxval=Sol[i][j];
//            if( minval> (Sol[i][j] - exactvalues[0]) ) minval=Sol[i][j] - exactvalues[0];
 //          if( maxval< fabs( exactvalues[0] - Sol[i][j]) ) maxval=fabs( exactvalues[0] - Sol[i][j]);
         }
        dat << endl;
        }

       dat.close();

       cout << "Min Value " << minval  << " Max Value " << maxval << endl;
      }
     } // for(l=0;l<N_SubSteps;l++)   

   } // while(TDatabase::TimeDB->CURRENTTIME< end_time

//   cout << "test ADI final" <<endl;
  CloseFiles();
  return 0;
}





