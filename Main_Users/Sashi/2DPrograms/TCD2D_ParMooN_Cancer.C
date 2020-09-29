// =======================================================================
//
// Purpose:     main program for scalar equations with new kernels of ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 08.08.2014

// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <SystemMatTimeScalar2D.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
#include <TimeUtilities.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

// =======================================================================
// include current example
// =======================================================================
// #include "../Examples/TCD_2D/exp.h"
// #include "../Examples/TCD_2D/SinCos1.h"
// #include "../Examples_All/TCD_2D/cancer.h"
#include "../Examples_All/TCD_2D/cancerlsg.h"
// #include "../Examples/TCD_2D_Cancer/test.h"
// #include "../Examples/TCD_2D_Cancer/testexp.h"
// #include "../Examples/TCD_2D_Cancer/test1.h"
// #include "../Examples/TCD_2D_Cancer/test2.h"
// #include "../Examples/TCD_2D_Cancer/testw.h"
//    #include "../Examples/TCD_2D_Cancer/testv.h"

void AssembleCancerDensity(TSquareMatrix2D *A, TFEFunction2D *u_FeFunction, TFEFunction2D *v_FeFunction)
{
  int i, j, k, l, m, n, begin, end, N_Cells, N_LocalUsedElements;
  int N_Points, N_U_LocDof, N_V_LocDof, TestDOF;
  int *N_BaseFunct, *BeginIndex, *GlobalNumbers, *RowPtr, *KCol, *DOF, *VBeginIndex, *VGlobalNumbers, *VDOF;
    
  double v, *weights, *xi, *eta, c0, lambda, *ValuesA;  
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double **uorig, **uxorig,**uyorig, *Orig, *x_Orig, *y_Orig, U_value;
  double **vorig, **vxorig,**vyorig, *V_Orig, *Vx_Orig, *Vy_Orig, V_value, Vx_value, Vy_value ;
  double Mult, mass, *U_Sol, *V_Sol;
  double test00, test10, test01, ansatz00;
  double uvfact, vgrad_test;
  double LocMatrixA[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];  
  
  TFESpace2D *fespace, *V_fespace;  
  TCollection *Coll;
  TBaseCell *Me;
  FE2D FEId, V_FEId;
  TFE2D *ele;
  FE2D LocalUsedElements[2];
  BaseFunct2D BaseFunct, V_BaseFunct, *BaseFuncts;
   
  bool *SecondDer;
  
  SecondDer = new bool[2];
  SecondDer[0] = FALSE;
  SecondDer[1] = FALSE;
  
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D(); 
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  
  // U space
  fespace = A->GetFESpace();    // outer phase
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  U_Sol = u_FeFunction->GetValues();
  
  BeginIndex =  fespace->GetBeginIndex();
  GlobalNumbers =  fespace->GetGlobalNumbers();

  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  ValuesA = A->GetEntries();

  // V space
  V_Sol = v_FeFunction->GetValues();
  V_fespace = v_FeFunction->GetFESpace2D();   
  VBeginIndex =  V_fespace->GetBeginIndex();
  VGlobalNumbers =  V_fespace->GetGlobalNumbers();
  
  
  N_LocalUsedElements = 2;
//   mass = 0;
  c0 = TDatabase::ParamDB->REACTOR_P2; // chi value
  lambda = TDatabase::ParamDB->REACTOR_P3;
  
  for(i=0;i<N_Cells;i++)
   { 
    Me = Coll->GetCell(i);
    FEId = fespace->GetFE2D(i, Me);
    V_FEId = V_fespace->GetFE2D(i, Me);

    // ======================================================
    // calculate values on original element
    // ======================================================
    LocalUsedElements[0] = FEId;
    LocalUsedElements[1] = V_FEId;
    
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                           Coll, Me, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);
  
    BaseFunct = BaseFuncts[FEId];
    N_U_LocDof = N_BaseFunct[FEId];
    
    V_BaseFunct = BaseFuncts[V_FEId];
    N_V_LocDof  = N_BaseFunct[V_FEId];
    
    uorig = TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);  
    uxorig = TFEDatabase2D::GetOrigElementValues(BaseFunct, D10); 
    uyorig = TFEDatabase2D::GetOrigElementValues(BaseFunct, D01);    
 
    vorig = TFEDatabase2D::GetOrigElementValues(V_BaseFunct, D00);  
    vxorig = TFEDatabase2D::GetOrigElementValues(V_BaseFunct, D10); 
    vyorig = TFEDatabase2D::GetOrigElementValues(V_BaseFunct, D01);    
    
    memset(LocMatrixA, 0, N_U_LocDof*N_U_LocDof*SizeOfDouble);

    DOF = GlobalNumbers + BeginIndex[i];    
    VDOF = VGlobalNumbers + VBeginIndex[i];         
    
    for(j=0;j<N_Points;j++)
     { 
      // calculate old u at this quad point 
      Orig = uorig[j];
      x_Orig = uxorig[j];     
      y_Orig = uyorig[j];   
       
      U_value = 0.; 
      for(l=0;l<N_U_LocDof;l++)
       {
        U_value += U_Sol[DOF[l]] * Orig[l];
       } 
       
      // calculate old v at this quad point  
      V_Orig = vorig[j];
      Vx_Orig = vorig[j];
      Vy_Orig = vorig[j];      
      V_value = 0.; Vx_value = 0.; Vy_value = 0.;
      for(l=0;l<N_V_LocDof;l++)
       {
        v = V_Sol[VDOF[l]];
        V_value += v * V_Orig[l];
        Vx_value += v * Vx_Orig[l];
        Vy_value+= v * Vy_Orig[l];
       } 
       
      Mult = weights[j]*AbsDetjk[j];
      
      // assemble the local matrix
      for(k=0;k<N_U_LocDof;k++)
      {
       test00  = Orig[k];
       test10  = x_Orig[k];
       test01  = y_Orig[k];

       vgrad_test =   Mult*c0*(Vx_value*test10 + Vy_value*test01);
       uvfact = Mult*lambda*(1. - U_value - V_value)*test00;
       
        //cout<< " uref " << test0  <<endl;
        for(l=0;l<N_U_LocDof;l++)
        {
          ansatz00  = Orig[l];
          LocMatrixA[k*N_U_LocDof + l] -= (vgrad_test + uvfact)*ansatz00;
        }
      } //    for(k=0;k<N_U_LocDof;k++    
      
     // checking 
//      mass += Mult;
     } //    for(j=0;j<N_Points   

    //add local to global matrix     
     for(l=0;l<N_U_LocDof;l++)
     {
      TestDOF = DOF[l];
//       RHS[TestDOF] += LocRhs[l]; 

      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      
      for(n=begin;n<end;n++)
       {
       for(m=0;m<N_U_LocDof;m++)
        {
         if(KCol[n] == DOF[m])
          {
           ValuesA[n] +=LocMatrixA[l*N_U_LocDof+m];
//            cout << TestDOF  << ", " << DOF_LOW[m] <<" A: "<< LocMatrixA[l*N_U_LocDof+m] <<endl;   
           break;
          }
        } // for(m=0;m<N_U_LocDof
      } // for(n=begin;n<end;n++)
     } // for(l=0;l<N_U_LocDof
     
   } //for(i=0;i<N_Cells;i
  
 
//   cout << " AssembleCancerDensity " << mass <<endl;
//   exit(0);
  
} // AssembleCancerDensity



void AssembleWArhs(TSquareMatrix2D *A, TFEFunction2D *u_FeFunction, double *RHS)
{
  int i, j, k, l, m, n, begin, end, N_Cells, N_LocalUsedElements;
  int N_Points, N_U_LocDof, N_V_LocDof, TestDOF;
  int *N_BaseFunct, *BeginIndex, *GlobalNumbers, *RowPtr, *KCol, *DOF, *UBeginIndex, *UGlobalNumbers, *UDOF;
    
  double *weights, *xi, *eta, alpha, beta, *ValuesA;  
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double **uorig,  *Orig, *x_Orig, *y_Orig, U_value;
  double **vorig,  *V_Orig, V_value;
  double Mult, *U_Sol;
  double test00, ansatz00;
  double ufact, alphaU, alphaUbeta;
  double LocMatrixA[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D], LocRhs[MaxN_BaseFunctions2D];  
  
  TFESpace2D *fespace, *U_fespace;  
  TCollection *Coll;
  TBaseCell *Me;
  FE2D FEId, U_FEId;
  TFE2D *ele;
  FE2D LocalUsedElements[2];
  BaseFunct2D BaseFunct, U_BaseFunct, *BaseFuncts;
   
  bool *SecondDer;
  
  SecondDer = new bool[2];
  SecondDer[0] = FALSE;
  SecondDer[1] = FALSE;
  
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D(); 
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  
  // W space
  fespace = A->GetFESpace();   
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  
  BeginIndex =  fespace->GetBeginIndex();
  GlobalNumbers =  fespace->GetGlobalNumbers();

  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  ValuesA = A->GetEntries();

  // U space
  U_Sol = u_FeFunction->GetValues();
  U_fespace = u_FeFunction->GetFESpace2D();
  UBeginIndex =  U_fespace->GetBeginIndex();
  UGlobalNumbers =  U_fespace->GetGlobalNumbers();
  
  
  N_LocalUsedElements = 2;
  alpha = TDatabase::ParamDB->REACTOR_P5;
  beta = TDatabase::ParamDB->REACTOR_P6;
  
  for(i=0;i<N_Cells;i++)
   { 
    Me = Coll->GetCell(i);
    
    FEId = fespace->GetFE2D(i, Me);
    U_FEId = U_fespace->GetFE2D(i, Me);

    // ======================================================
    // calculate values on original element
    // ======================================================
    LocalUsedElements[0] = FEId;
    LocalUsedElements[1] = U_FEId;
    
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                           Coll, Me, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);
  
    BaseFunct = BaseFuncts[FEId];
    N_V_LocDof = N_BaseFunct[FEId];
    
    U_BaseFunct = BaseFuncts[U_FEId];
    N_U_LocDof  = N_BaseFunct[U_FEId];
    
    vorig = TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);     
    uorig = TFEDatabase2D::GetOrigElementValues(U_BaseFunct, D00);  
 
    memset(LocMatrixA, 0, N_V_LocDof*N_V_LocDof*SizeOfDouble);
    memset(LocRhs, 0, N_V_LocDof*SizeOfDouble);

    DOF = GlobalNumbers + BeginIndex[i];       
    UDOF = UGlobalNumbers + UBeginIndex[i];    
     
    for(j=0;j<N_Points;j++)
     {
      // calculate old u at this quad point 
      Orig = uorig[j];  
       
      U_value = 0.; 
      for(l=0;l<N_U_LocDof;l++)
       {
        U_value += U_Sol[UDOF[l]] * Orig[l];
       } 
       
      // weights 
      Mult = weights[j]*AbsDetjk[j];
      alphaU     = Mult*alpha*U_value;
      alphaUbeta = Mult*(alpha*U_value + beta);

      // assemble the local matrix and rhs
      for(k=0;k<N_V_LocDof;k++)
      {
       test00  = Orig[k];
       LocRhs[k] += alphaU*test00; 
       
       ufact = alphaUbeta*test00;
       
       for(l=0;l<N_V_LocDof;l++)
        {
          ansatz00  = Orig[l];
          LocMatrixA[k*N_V_LocDof + l] += ufact*ansatz00;
        }
      } //    for(k=0;k<N_V_LocDof;k++    

     } //    for(j=0;j<N_Points   

    
    //add local to global matrix     
     for(l=0;l<N_V_LocDof;l++)
     {
      TestDOF = DOF[l];
      RHS[TestDOF] += LocRhs[l]; 

      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      
      for(n=begin;n<end;n++)
       {
       for(m=0;m<N_V_LocDof;m++)
        {
         if(KCol[n] == DOF[m])
          {
           ValuesA[n] +=LocMatrixA[l*N_V_LocDof+m];
//            cout << TestDOF  << ", " << DOF_LOW[m] <<" A: "<< LocMatrixA[l*N_V_LocDof+m] <<endl;   
           break;
          }
        } // for(m=0;m<N_V_LocDof
      } // for(n=begin;n<end;n++)
     } // for(l=0;l<N_V_LocDof
     
   } //for(i=0;i<N_Cells;i
  
 
//   cout << " AssembleVrhs "  <<endl;
//   exit(0);
  
} // AssembleWArhs



void AssembleVMat(TSquareMatrix2D *A, TFEFunction2D *w_FeFunction, TFEFunction2D *u_FeFunction, TFEFunction2D *v_FeFunction)
{
  int i, j, k, l, m, n, begin, end, N_Cells, N_LocalUsedElements;
  int N_Points, N_W_LocDof, N_V_LocDof, TestDOF, N_U_LocDof;
  int *N_BaseFunct, *BeginIndex, *GlobalNumbers, *RowPtr, *KCol, *DOF, *WBeginIndex, *WGlobalNumbers, *WDOF;
  int *UBeginIndex, *UGlobalNumbers, *UDOF;
  
  double *weights, *xi, *eta, coeff, *ValuesA, *U_Sol, mu_a;  
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double **uorig, **worig, *Orig, *x_Orig, *y_Orig, W_value, U_value;
  double **vorig,  *V_Orig, V_value;
  double Mult, *W_Sol, *V_Sol;
  double test00, ansatz00;
  double wfact;
  double LocMatrixA[MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];  
  
  TFESpace2D *fespace, *W_fespace, *U_fespace;  
  TCollection *Coll;
  TBaseCell *Me;
  FE2D FEId, W_FEId, U_FEId;
  TFE2D *ele;
  FE2D LocalUsedElements[3];
  BaseFunct2D BaseFunct, W_BaseFunct, U_BaseFunct, *BaseFuncts;
   
  bool *SecondDer;
  
  SecondDer = new bool[2];
  SecondDer[0] = FALSE;
  SecondDer[1] = FALSE;
  
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D(); 
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  
  // V space
  V_Sol = v_FeFunction->GetValues();
  fespace = A->GetFESpace();   
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  
  BeginIndex =  fespace->GetBeginIndex();
  GlobalNumbers =  fespace->GetGlobalNumbers();

  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  ValuesA = A->GetEntries();

  // W space
  W_Sol = w_FeFunction->GetValues();
  W_fespace = w_FeFunction->GetFESpace2D();
  WBeginIndex =  W_fespace->GetBeginIndex();
  WGlobalNumbers =  W_fespace->GetGlobalNumbers();
  
  //U space  
  U_Sol = u_FeFunction->GetValues();
  U_fespace = u_FeFunction->GetFESpace2D();
  UBeginIndex =  U_fespace->GetBeginIndex();
  UGlobalNumbers =  U_fespace->GetGlobalNumbers();
  
 
  N_LocalUsedElements = 3;
  coeff = TDatabase::ParamDB->REACTOR_P4;
  mu_a = TDatabase::ParamDB->REACTOR_P7;
   
  for(i=0;i<N_Cells;i++)
   { 
    Me = Coll->GetCell(i);
    
    FEId = fespace->GetFE2D(i, Me);
    W_FEId = W_fespace->GetFE2D(i, Me);
    U_FEId = U_fespace->GetFE2D(i, Me);
    
    // ======================================================
    // calculate values on original element
    // ======================================================
    LocalUsedElements[0] = FEId;
    LocalUsedElements[1] = W_FEId;
    LocalUsedElements[2] = U_FEId;
    
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                           Coll, Me, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);
  
    BaseFunct = BaseFuncts[FEId];
    N_V_LocDof = N_BaseFunct[FEId];
    
    W_BaseFunct = BaseFuncts[W_FEId];
    N_W_LocDof  = N_BaseFunct[W_FEId];
    
    U_BaseFunct = BaseFuncts[U_FEId];
    N_U_LocDof  = N_BaseFunct[U_FEId];   
    
    vorig = TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);     
    uorig = TFEDatabase2D::GetOrigElementValues(U_BaseFunct, D00);  
    worig = TFEDatabase2D::GetOrigElementValues(W_BaseFunct, D00);  
    
    memset(LocMatrixA, 0, N_V_LocDof*N_V_LocDof*SizeOfDouble);
    
    DOF = GlobalNumbers + BeginIndex[i];    
    WDOF = WGlobalNumbers + WBeginIndex[i];    
    UDOF = UGlobalNumbers + UBeginIndex[i];  
    
    for(j=0;j<N_Points;j++)
     {
      // calculate old w at this quad point 
      Orig = worig[j];     
      W_value = 0.; 
      for(l=0;l<N_W_LocDof;l++)
       {
        W_value += W_Sol[WDOF[l]] * Orig[l];
       } 
       
      // calculate old w at this quad point 
      Orig = uorig[j];       
      U_value = 0.; 
      for(l=0;l<N_U_LocDof;l++)
       {
        U_value += U_Sol[UDOF[l]] * Orig[l];
       } 
       
     // calculate old v at this quad point 
      Orig = vorig[j];       
      V_value = 0.; 
      for(l=0;l<N_V_LocDof;l++)
       {
        V_value += V_Sol[DOF[l]] * Orig[l];
       }     
          
      // weights and coeff
      Mult = weights[j]*AbsDetjk[j]*( coeff*W_value - mu_a*(1.-U_value-V_value) );
      
      // assemble the local matrix and rhs
      for(k=0;k<N_V_LocDof;k++)
      {
       test00  = Orig[k];
       wfact = Mult*test00;
       
       for(l=0;l<N_V_LocDof;l++)
        {
          ansatz00  = Orig[l];
          LocMatrixA[k*N_V_LocDof + l] += wfact*ansatz00;
        }
      } //    for(k=0;k<N_V_LocDof;k++    
     } //    for(j=0;j<N_Points   
     
    //add local to global matrix     
     for(l=0;l<N_V_LocDof;l++)
     {
      TestDOF = DOF[l];
      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      
      for(n=begin;n<end;n++)
       {
       for(m=0;m<N_V_LocDof;m++)
        {
         if(KCol[n] == DOF[m])
          {
           ValuesA[n] +=LocMatrixA[l*N_V_LocDof+m];
//            cout << TestDOF  << ", " << DOF_LOW[m] <<" A: "<< LocMatrixA[l*N_V_LocDof+m] <<endl;   
           break;
          }
        } // for(m=0;m<N_V_LocDof
      } // for(n=begin;n<end;n++)
     } // for(l=0;l<N_V_LocDof
     
   } //for(i=0;i<N_Cells;i
//   cout << " AssembleVrhs "  <<endl;
//   exit(0);
  
} // AssembleVMat


int main(int argc, char* argv[])
{
// ======================================================================
// variable declaration
// ======================================================================
  int i, j, l, m, N_SubSteps, ORDER, N_Cells, N_DOF,  N_V_DOF, N_W_DOF, img=1;
  int N_Active, N_V_Active, N_W_Active, Max_It_scalar, Min_It_scalar;
  
  double *sol, *oldsol, *rhs, *oldrhs, t1, t2, errors[5];
  double tau, end_time, *defect, olderror, olderror1;
  double *V_sol, *V_oldsol, *V_rhs, *V_oldrhs, *V_defect;
  double *W_sol, *W_oldsol, *W_rhs, *W_oldrhs, *W_defect, residual_scalar, oldresidual_scalar, limit_scalar;
  
  double Linfty;
  double Parameters[10], hmin, hmax;
  
  bool UpdateStiffnessMat, UpdateRhs,  ConvectionFirstTime;
  
  TDomain *Domain;
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D(); 
  TCollection *coll;
  TFESpace2D *Scalar_FeSpace, *Scalar_V_FeSpace, *Scalar_W_FeSpace, *fesp[1];
  TFEFunction2D *Scalar_FeFunction, *Scalar_V_FeFunction, *Scalar_W_FeFunction;
  TOutput2D *Output;
  TSystemMatTimeScalar2D *SystemMatrix, *V_SystemMatrix, *W_SystemMatrix;
  TAuxParam2D *aux;
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
  TSquareMatrix2D *V_A_Matrix;
  
  const char vtkdir[] = "VTK"; 
  char *PsBaseName, *VtkBaseName, *GEO;
  char Name[] = "name";
  char Description[] = "description";
  char UString[] = "u";
  char VString[] = "v";
  char WString[] = "w";
    
  std::ostringstream os;
  os << " ";   
    
  mkdir(vtkdir, 0777);
  
// ======================================================================
// set the database values and generate mesh
// ======================================================================    
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  Domain = new TDomain(argv[1]);  

  /* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
  // standard mesh  
  GEO = TDatabase::ParamDB->GEOFILE;
  Domain->ReadGeo(GEO);

  
  Domain->ConvertQuadToTri(TDatabase::ParamDB->CONVERT_QUAD_TO_TRI);
    
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
  
//=========================================================================
// construct all finite element spaces
//=========================================================================
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  
  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  OutPut("N_Cells (space) : " << N_Cells <<endl);
  
  // u fespaces for scalar equation 
   Scalar_FeSpace  =  new TFESpace2D(coll, Name, Description, BoundCondition, ORDER, NULL); 
   N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
   N_Active =  Scalar_FeSpace->GetActiveBound();
   OutPut("dof all      : "<< setw(10) << N_DOF  << endl);
   
  // v fespaces for scalar equation 
   Scalar_V_FeSpace  =  new TFESpace2D(coll, Name, Description, V_BoundCondition, ORDER, NULL); 
   N_V_DOF = Scalar_V_FeSpace->GetN_DegreesOfFreedom();
   N_V_Active =  Scalar_V_FeSpace->GetActiveBound();
   OutPut("dof all      : "<< setw(10) << N_V_DOF  << endl);
      
  // w fespaces for scalar equation 
   Scalar_W_FeSpace  =  new TFESpace2D(coll, Name, Description, W_BoundCondition, ORDER, NULL);   
   N_W_DOF = Scalar_W_FeSpace->GetN_DegreesOfFreedom();
   N_W_Active =  Scalar_W_FeSpace->GetActiveBound();
   OutPut("dof all      : "<< setw(10) << N_W_DOF  << endl);   
   
//======================================================================
// construct all finite element functions
//======================================================================
   
   // u function
    sol = new double[N_DOF];
    oldsol = new double[N_DOF];   
    rhs = new double[N_DOF];
    oldrhs   = new double[N_DOF];
    defect = new double[N_DOF];
    
    memset(sol, 0, N_DOF*SizeOfDouble);
    memset(rhs, 0, N_DOF*SizeOfDouble);

    Scalar_FeFunction = new TFEFunction2D(Scalar_FeSpace, UString, UString, sol, N_DOF); 
    
    //interpolate the initial value
    Scalar_FeFunction->Interpolate(InitialCondition);
    
    
   // v function
    V_sol = new double[N_V_DOF];
    V_oldsol = new double[N_V_DOF];    
    V_rhs = new double[N_V_DOF];
    V_oldrhs = new double[N_V_DOF];
    V_defect = new double[N_V_DOF];
    
    memset(V_sol, 0, N_V_DOF*SizeOfDouble);
    memset(V_rhs, 0, N_V_DOF*SizeOfDouble);

    Scalar_V_FeFunction = new TFEFunction2D(Scalar_V_FeSpace, VString, VString, V_sol, N_V_DOF); 
    
    //interpolate the initial value
    Scalar_V_FeFunction->Interpolate(V_InitialCondition);
    
   // w function
    W_sol = new double[N_W_DOF];
    W_oldsol = new double[N_W_DOF];    
    W_rhs = new double[N_W_DOF];
    W_oldrhs = new double[N_W_DOF];
    W_defect = new double[N_W_DOF];
    
    memset(W_sol, 0, N_W_DOF*SizeOfDouble);
    memset(W_rhs, 0, N_W_DOF*SizeOfDouble);

    Scalar_W_FeFunction = new TFEFunction2D(Scalar_W_FeSpace, WString, WString, W_sol, N_W_DOF); 
    
    //interpolate the initial value
    Scalar_W_FeFunction->Interpolate(W_InitialCondition);

//======================================================================
// SystemMatrix construction and solution
//======================================================================  
    //U  SYSTEM
    // Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT 
    SystemMatrix = new TSystemMatTimeScalar2D(Scalar_FeSpace, GALERKIN, DIRECT);
    
    // initilize the system matrix with the functions defined in Example file
    SystemMatrix->Init(BilinearCoeffs, BoundCondition, BoundValue);
       
    // assemble the system matrix with given aux, sol and rhs 
    // aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
    // otherwise, just pass with NULL 
    SystemMatrix->AssembleMRhs(NULL, sol, rhs); 
   
// 		cout << Ddot( N_DOF, rhs, rhs) << endl;
// 		exit(0);    
    
    // V SYSTEM 
     V_SystemMatrix = new TSystemMatTimeScalar2D(Scalar_V_FeSpace, GALERKIN, DIRECT);
     V_SystemMatrix->Init(V_BilinearCoeffs, V_BoundCondition, V_BoundValue);    
     V_SystemMatrix->AssembleMRhs(NULL, V_sol, V_rhs);     
     
    //W SYSTEM
     W_SystemMatrix = new TSystemMatTimeScalar2D(Scalar_W_FeSpace, GALERKIN, DIRECT);
     W_SystemMatrix->Init(W_BilinearCoeffs, W_BoundCondition, W_BoundValue);    
     W_SystemMatrix->AssembleMRhs(NULL, W_sol, W_rhs); 

//======================================================================
// produce outout at t=0
//======================================================================
    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
    Output = new TOutput2D(3, 3, 1, 1, Domain);

    Output->AddFEFunction(Scalar_FeFunction);
    Output->AddFEFunction(Scalar_V_FeFunction);
    Output->AddFEFunction(Scalar_W_FeFunction);
    
//   Scalar_FeFunction->Interpolate(Exact);   
    if(TDatabase::ParamDB->WRITE_VTK)
     {
      os.seekp(std::ios::beg);
       if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      img++;
     }   

//     exit(0);
     
     
// measure errors to known solution u
    /*if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      fesp[0] = Scalar_FeSpace;       
      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      
     for(j=0;j<5;j++)
       errors[j] = 0;
     
      Scalar_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, BilinearCoeffs, aux, 1, fesp, errors);  
      olderror = errors[0];
      olderror1 = errors[1]; 
     
      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
      OutPut(" L2: " << errors[0]);
      OutPut(" H1-semi: " << errors[1] << endl);     
     Linfty=errors[0];
     } *//*if(TDatabase::ParamDB->MEASURE_ERRORS)  
         
     exit(0);*/
     
 // measure errors to known solution w
//      if(TDatabase::ParamDB->MEASURE_ERRORS)
//      {
//       fesp[0] =Scalar_W_FeSpace;       
//       aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
//       
//      for(j=0;j<5;j++)
//        errors[j] = 0;
//      
//       Scalar_W_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, BilinearCoeffs, aux, 1, fesp, errors);  
//       olderror = errors[0];
//       olderror1 = errors[1]; 
//      
//       OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
//       OutPut(" L2: " << errors[0]);
//       OutPut(" H1-semi: " << errors[1] << endl);     
//      Linfty=errors[0];
//      } 
 
 
 
 // measure errors to known solution v
     if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      fesp[0] =Scalar_V_FeSpace;       
      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      
     for(j=0;j<5;j++)
       errors[j] = 0;
     
      Scalar_V_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, BilinearCoeffs, aux, 1, fesp, errors);  
      olderror = errors[0];
      olderror1 = errors[1]; 
     
      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
      OutPut(" L2: " << errors[0]);
      OutPut(" H1-semi: " << errors[1] << endl);     
     Linfty=errors[0];
     } 
     
     
//======================================================================
// calculating h value 
//======================================================================    
//   coll->GetHminHmax(&hmin,&hmax);
//   OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
//   
//   TDatabase::TimeDB->TIMESTEPLENGTH =  hmax;
//   
//  exit(0);



//======================================================================
// time disc loop
//======================================================================    
   // parameters for time stepping scheme
   m = 0;
   N_SubSteps = GetN_SubSteps();
   end_time = TDatabase::TimeDB->ENDTIME; 
   Max_It_scalar = (int)TDatabase::ParamDB->REACTOR_P10;
   Min_It_scalar = TDatabase::ParamDB->REACTOR_P11;
   limit_scalar = TDatabase::ParamDB->REACTOR_P12;
   
   UpdateStiffnessMat = TRUE; //check BilinearCoeffs in example file
   UpdateRhs = FALSE; //check BilinearCoeffs in example file
   ConvectionFirstTime=TRUE;
   
   // time loop starts
   while(TDatabase::TimeDB->CURRENTTIME< end_time)
    {                                               // time cycle
     m++;
     TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

     for(l=0;l<N_SubSteps;l++) // sub steps of fractional step theta
      {
       SetTimeDiscParameters();

      if(m==1)
       {
        OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
        OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
        OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
        OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
       }

      tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;
   
      OutPut(endl << "CURRENT TIME: ");
      OutPut(TDatabase::TimeDB->CURRENTTIME << endl);   
     
      //copy rhs to oldrhs
      memcpy(oldrhs, rhs, N_DOF*SizeOfDouble); 
      memcpy(V_oldrhs, V_rhs, N_V_DOF*SizeOfDouble);   
      memcpy(W_oldrhs, W_rhs, N_W_DOF*SizeOfDouble);       
      
      //copy sol to oldsol   
      memcpy(oldsol, sol, N_DOF*SizeOfDouble); 
      memcpy(V_oldsol, V_sol, N_V_DOF*SizeOfDouble);   
      memcpy(W_oldsol, W_sol, N_W_DOF*SizeOfDouble);       
      
      // GAUSS-SEIDEL TYPE LINEARIZIZATION
      for(j=0;j<Max_It_scalar;j++)
       {   
//=======================================================================       
//U SYSTEM
//======================================================================= 
        // unless the stiffness matrix or rhs change in time, it is enough to assemble only once at the begning   
        if(UpdateStiffnessMat || UpdateRhs ||  ConvectionFirstTime)
         {    
          SystemMatrix->AssembleARhs(NULL, sol, rhs);

          AssembleCancerDensity(SystemMatrix->GetAMatrix(), Scalar_FeFunction, Scalar_V_FeFunction);

//     M:= M + (tau*TDatabase::TimeDB->THETA1)*A
//     rhs: =(tau*TDatabase::TimeDB->THETA4)*rhs +(tau*TDatabase::TimeDB->THETA3)*oldrhs + [ M - (tau*TDatabase::TimeDB->THETA2)A]*oldsol
//     note! sol contains only the previous time step value, so just pass sol for oldsol
          SystemMatrix->AssembleSystMat(oldrhs, oldsol, rhs, sol);
         }
              
        residual_scalar = SystemMatrix->GetResidual(sol);
        OutPut("Scalar nonlinear step " << setw(3) << j);
        OutPut(setw(14) << residual_scalar); // sqrt of residual_scalar is alread done in ScalarDefect

        if (j>0)
         { 
          if(fabs(oldresidual_scalar)>0)
          OutPut(setw(14) <<  residual_scalar/oldresidual_scalar ); 
          OutPut(endl);
         }
        else
         { OutPut(endl); }
   

        oldresidual_scalar = residual_scalar;
   
        if( ((residual_scalar<=limit_scalar)||(j==Max_It_scalar-1))  && (j>=Min_It_scalar) )
         { 
          if(UpdateStiffnessMat || UpdateRhs)
           {         
            SystemMatrix->RestoreMassMat();
           }   
          break;
         }        
   
        // solve the system matrix 
        SystemMatrix->Solve(sol, rhs);

        // restore the mass matrix for the next time step    
        // unless the stiffness matrix or rhs change in time, it is not necessary to assemble the system matrix in every time step
        if(UpdateStiffnessMat || UpdateRhs)
         {         
          SystemMatrix->RestoreMassMat();
         }
         
//=======================================================================       
//W SYSTEM
//======================================================================= 
        if(UpdateStiffnessMat || UpdateRhs ||  ConvectionFirstTime)
         {    
          W_SystemMatrix->AssembleARhs(NULL, W_sol, W_rhs);
          AssembleWArhs(W_SystemMatrix->GetAMatrix(), Scalar_FeFunction, W_rhs);
          W_SystemMatrix->AssembleSystMat(W_oldrhs,  W_oldsol, W_rhs, W_sol);
       
          ConvectionFirstTime = FALSE;
         }
     
        W_SystemMatrix->Solve(W_sol, W_rhs);

        if(UpdateStiffnessMat || UpdateRhs)
         {         
          W_SystemMatrix->RestoreMassMat();
         }     
//=======================================================================       
// V SYSTEM
//======================================================================= 
//         no need a=to assemble A matrix
         
         V_SystemMatrix->AssembleARhs(NULL, V_sol, V_rhs);
         
        V_A_Matrix = V_SystemMatrix->GetAMatrix();
        V_A_Matrix->Reset();
        AssembleVMat(V_A_Matrix, Scalar_W_FeFunction, Scalar_FeFunction, Scalar_V_FeFunction);
        V_SystemMatrix->AssembleSystMat(V_oldrhs, V_oldsol, V_rhs, V_sol);
	
// 	exit(0);

        V_SystemMatrix->Solve(V_sol, V_rhs);

        if(UpdateStiffnessMat || UpdateRhs)
         {         
          V_SystemMatrix->RestoreMassMat();
         }  
//=======================================================================      
     } //for(j=0;j<Max_It_scalar;j++)
//       exit(0);
    } // for(l=0;l<N_SubSteps;l++) 
//======================================================================
// produce outout
//======================================================================
    if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
     if(TDatabase::ParamDB->WRITE_VTK)
      {
       os.seekp(std::ios::beg);
        if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());
       img++;
      }
      
//======================================================================
// measure errors to known solution
//======================================================================    
    /*if(TDatabase::ParamDB->MEASURE_ERRORS)
     {      
      Scalar_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, BilinearCoeffs, aux, 1, fesp, errors);

      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
      OutPut(" L2: " << errors[0]);
      OutPut(" H1-semi: " << errors[1] << endl);

      errors[3] += (errors[0]*errors[0] + olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      olderror = errors[0];
      OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(errors[3]) << " ");

      errors[4] += (errors[1]*errors[1] +olderror1 * olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      OutPut( "L2(0,T;H1) " << sqrt(errors[4]) << endl);
      olderror1 = errors[1];     
      
      if(Linfty<errors[0])
	Linfty=errors[0];
      
     } *///  if(TDatabase::ParamDB->MEASURE_ERRORS)
     
/*     if(m==2)
     exit(0); */ 
 // } // while(TDatabase::TimeDB->CURRENTTIME< end_time)
  
     
//======================================================================
// measure errors to known solution w
//======================================================================  
//       if(TDatabase::ParamDB->MEASURE_ERRORS)
//      {      
//       Scalar_W_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, BilinearCoeffs, aux, 1, fesp, errors);
// 
//       OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
//       OutPut(" L2: " << errors[0]);
//       OutPut(" H1-semi: " << errors[1] << endl);
// 
//       errors[3] += (errors[0]*errors[0] + olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
//       olderror = errors[0];
//       OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(errors[3]) << " ");
// 
//       errors[4] += (errors[1]*errors[1] +olderror1 * olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
//       OutPut( "L2(0,T;H1) " << sqrt(errors[4]) << endl);
//       olderror1 = errors[1];     
//       
//       if(Linfty<errors[0])
// 	Linfty=errors[0];
//       
//      } 

    //}

//======================================================================
// measure errors to known solution v
//======================================================================  
      if(TDatabase::ParamDB->MEASURE_ERRORS)
     {      
      Scalar_V_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, BilinearCoeffs, aux, 1, fesp, errors);

      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
      OutPut(" L2: " << errors[0]);
      OutPut(" H1-semi: " << errors[1] << endl);

      errors[3] += (errors[0]*errors[0] + olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      olderror = errors[0];
      OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(errors[3]) << " ");

      errors[4] += (errors[1]*errors[1] +olderror1 * olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      OutPut( "L2(0,T;H1) " << sqrt(errors[4]) << endl);
      olderror1 = errors[1];     
      
      if(Linfty<errors[0])
	Linfty=errors[0];
      
     } 
     
//                if(m==2)
//       exit(0);
    }
//   OutPut( "Linfty " << Linfty << endl);
//======================================================================
// produce final outout
//======================================================================
     
     if(TDatabase::ParamDB->WRITE_VTK)
      {
       os.seekp(std::ios::beg);
        if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());
       img++;
      }
      
      
      
  CloseFiles();
  
  return 0;
} // end main









