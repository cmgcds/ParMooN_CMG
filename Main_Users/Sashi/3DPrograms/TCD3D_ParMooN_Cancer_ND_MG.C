// =======================================================================
//
// Purpose:     main program for time-dependent scalar cancer equation with new kernels of ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 24.01.15

// =======================================================================
 
#include <Domain.h>
#include <Database.h>
#include <TSystemTCD3D.h>
#include <FEDatabase3D.h>
#include <FESpace3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble3D.h>
#include <Output3D.h>
#include <LinAlg.h>
// #include <TCD3DErrorEstimator.h>
#include <MainUtilities.h>
// #include <TimeUtilities.h>
#include <TimeDiscRout.h>
#include <Solver.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>


#include <tetgen.h>
#include <GridCell.h>
#include <MacroCell.h>
#include <BdPlane.h>
#include <BdSphere.h>
#include <IsoBoundFace.h>
#include <InterfaceJoint3D.h>
#include <IsoInterfaceJoint3D.h>

#ifdef _MPI
#include "mpi.h"
#include <MeshPartition.h>
//#include <MeshPartition2D.h>
// #include <ParFECommunicator3D.h>
// #include <MumpsSolver.h>
// #include <ParVector3D.h>
// #include <ParVectorNSE3D.h>
// #include <Scalar_ParSolver.h>
#endif

double bound = 0;
double timeC = 0;
// =======================================================================
// include current example
// =======================================================================
// #include "../Examples/TCD_3D/Sin4.h"
// #include "../Examples/TCD_3D/ConstTSmooth.h"
// #include "../Examples/TCD_3D/ConstT.h"
// #include "../Examples/TCD_3D/amc.h"
// #include "../Examples_All/TCD_3D/cancer3D.h"
// #include "../Examples/TCD_3D/cancer3D.h"
// #include "../Examples/TCD_3D/cancerconv.h"
// #include "../Examples/TCD_3D/cancertest.h"
   #include "../../Examples/TCD_3D/cancer_exp_3d_nd.h"
// #include "../../Examples/TCD_3D/cancer_exp_3d.h"
// #include "../Examples/TCD_3D/cancerconv1.h"
// #include "../Examples/TCD_3D/cancer_nd.h"

void AssembleCancerDensity(int Levels, TSquareMatrix3D **A, TFEFunction3D **u_FeFunction, TFEFunction3D **v_FeFunction, TFEFunction3D **w_FeFunction)
{
  int ii, i, j, k, l, m, n, begin, end, N_Cells, N_LocalUsedElements;
  int N_Points, N_U_LocalDof, N_V_LocalDof, N_W_LocDof, TestDOF;
  int *N_BaseFunct, *BeginIndex, *GlobalNumbers, *RowPtr, *KCol, *DOF, *VBeginIndex, *VGlobalNumbers, *VDOF, *WBeginIndex, *WGlobalNumbers, *WDOF;
  
  double v, *weights, *xi, *eta, *zeta, c0, lambda, *ValuesA, d1, d2, num, a;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  double **uorig, **uxorig, **uyorig, **uzorig, *Orig, *x_Orig, *y_Orig, *z_Orig, U_value, Ux_value, Uy_value, Uz_value;
  double **vorig, **vxorig, **vyorig, **vzorig, *V_Orig, *Vx_Orig, *Vy_Orig, *Vz_Orig, V_value, Vx_value, Vy_value, Vz_value;
  double **worig, **wxorig, **wyorig, **wzorig, *W_Orig, *Wx_Orig, *Wy_Orig, *Wz_Orig, W_value, Wx_value, Wy_value, Wz_value;
  double Mult, mass, *U_Sol, *V_Sol, *W_Sol;
  double test000, test100, test010, test001, ansatz000, ansatz001, ansatz010, ansatz100, val1, val2;
  double uvfact, vgrad_test, diff_fact, cancerdiffusion;
  double LocMatrixA[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  
  TFESpace3D *fespace, *V_fespace, *W_fespace;
  TCollection *Coll;
  TBaseCell *Me;
  FE3D FEId, V_FEId, W_FEId;
  TFE3D *ele;
  FE3D LocalUsedElements[2];
  BaseFunct3D BaseFunct, V_BaseFunct, W_BaseFunct, *BaseFuncts;
  
  bool *SecondDer;
  
  SecondDer = new bool[2];
  SecondDer[0] = FALSE;
  SecondDer[1] = FALSE;
  
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();
  
//   U solution Space
  for(ii=0;ii<Levels; ii++)
   {
    fespace = A[ii]->GetFESpace();
    Coll = fespace->GetCollection();
    N_Cells = Coll->GetN_Cells();
    U_Sol = u_FeFunction[ii]->GetValues();
  
    BeginIndex = fespace->GetBeginIndex();
    GlobalNumbers = fespace->GetGlobalNumbers();
  
    RowPtr = A[ii]->GetRowPtr();
    KCol = A[ii]->GetKCol();
    ValuesA = A[ii]->GetEntries();
  
    //V solution Space
    V_Sol = v_FeFunction[ii]->GetValues();
    V_fespace = v_FeFunction[ii]->GetFESpace3D();
    VBeginIndex = V_fespace->GetBeginIndex();
    VGlobalNumbers = V_fespace->GetGlobalNumbers();
    
    //W solution space
    W_Sol = w_FeFunction[ii]->GetValues();
    W_fespace = w_FeFunction[ii]->GetFESpace3D();
    WBeginIndex = W_fespace->GetBeginIndex();
    WGlobalNumbers = W_fespace->GetGlobalNumbers();
  
    N_LocalUsedElements = 2;
    c0 = TDatabase::ParamDB->REACTOR_P2; //chi
    lambda = TDatabase::ParamDB->REACTOR_P3;
    d1= TDatabase::ParamDB->REACTOR_P13; //diffusion coeff
    
    d2= TDatabase::ParamDB->REACTOR_P14; //stepin diff 
    num=TDatabase::ParamDB->REACTOR_P15; //stepin diff
    a=TDatabase::ParamDB->REACTOR_P16; //stepin diff 
  
    for(i=0;i<N_Cells;i++)
     {
      Me = Coll->GetCell(i);
      FEId = fespace->GetFE3D(i,Me);
      V_FEId = V_fespace->GetFE3D(i,Me);
      W_FEId = W_fespace->GetFE3D(i,Me);
    
      //Calculate values on original element
      LocalUsedElements[0]=FEId;
      LocalUsedElements[1]=V_FEId;
      LocalUsedElements[2]=W_FEId;
    
      TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                             Coll, Me, SecondDer, 
                             N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);
    
      BaseFunct = BaseFuncts[FEId];
      N_U_LocalDof = N_BaseFunct[FEId];
    
      V_BaseFunct = BaseFuncts[V_FEId];
      N_V_LocalDof = N_BaseFunct[V_FEId];
      
      W_BaseFunct = BaseFuncts[W_FEId];
      N_W_LocDof = N_BaseFunct[W_FEId];
    
      uorig = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
      uxorig = TFEDatabase3D::GetOrigElementValues(BaseFunct, D100);
      uyorig = TFEDatabase3D::GetOrigElementValues(BaseFunct, D010);
      uzorig = TFEDatabase3D::GetOrigElementValues(BaseFunct, D001);
    
      vorig = TFEDatabase3D::GetOrigElementValues(V_BaseFunct, D000);
      vxorig = TFEDatabase3D::GetOrigElementValues(V_BaseFunct, D100);
      vyorig = TFEDatabase3D::GetOrigElementValues(V_BaseFunct, D010);
      vzorig = TFEDatabase3D::GetOrigElementValues(V_BaseFunct, D001);
      
      worig = TFEDatabase3D::GetOrigElementValues(W_BaseFunct, D000);
      wxorig = TFEDatabase3D::GetOrigElementValues(W_BaseFunct, D100);
      wyorig = TFEDatabase3D::GetOrigElementValues(W_BaseFunct, D010);
      wzorig = TFEDatabase3D::GetOrigElementValues(W_BaseFunct, D001);
    
      memset(LocMatrixA, 0, N_U_LocalDof*N_U_LocalDof*SizeOfDouble);
    
      DOF = GlobalNumbers + BeginIndex[i];
      VDOF = VGlobalNumbers + VBeginIndex[i];
      WDOF = WGlobalNumbers + WBeginIndex[i];
    
      for(j=0;j<N_Points;j++)
       {
        Orig = uorig[j];
        x_Orig = uxorig[j];
        y_Orig = uyorig[j];
        z_Orig = uzorig[j];
      
        U_value = 0.; Ux_value = 0.; Uy_value = 0.; Uz_value = 0.;
        for(l=0; l<N_U_LocalDof;l++)
         {
          U_value += U_Sol[DOF[l]]* Orig[l];
         }
      
        V_Orig = vorig[j];
        Vx_Orig = vxorig[j];
        Vy_Orig = vyorig[j];
        Vz_Orig = vzorig[j];
      
        V_value = 0.; Vx_value = 0.; Vy_value = 0.; Vz_value= 0.;
        for(l=0; l<N_V_LocalDof;l++)
         {
          v = V_Sol[VDOF[l]];
          V_value += v * V_Orig[l];
          Vx_value += v * Vx_Orig[l];
          Vy_value += v * Vy_Orig[l];
          Vz_value += v * Vz_Orig[l];
         }
      
       
       W_Orig = worig[j];
       
       W_value = 0.;
       for(l=0; l<N_W_LocDof; l++)
       {
	 W_value += W_Sol[WDOF[l]] * Orig[l];
       }
       
        Mult = weights[j] * AbsDetjk[j];
      
       //assembling the local matrix
       for(k=0;k<N_U_LocalDof;k++)
        {
         test000 = Orig[k];
         test100 = x_Orig[k];
         test010 = y_Orig[k];
         test001 = z_Orig[k];
//=================================================Diffusion terms===============================================
//Constant Diffusion
//            diff_fact = Mult * d1;	 
//Meral Diffusion
// 	 diff_fact = Mult * d1 * (1. / ( 1. + U_value * V_value));
//Anderson Diffusion
         diff_fact = Mult * d1 * W_value;	
//Stepin Diffusion
//           diff_fact = Mult * d1 + ( (d2 * pow(U_value, num))/ (a + pow(U_value,num)));	
//Chaplain-Lolas Diffusion
//           diff_fact = Mult * d1 * (1. + U_value * V_value);	  
//=================================================================================================================
	   
//================================================Haptotactic terms================================================	   
	   
//Direct Measurement	   
// 	 vgrad_test = Mult * c0 * (Vx_value*test100 + Vy_value*test010 + Vz_value*test001);
//Logarithmic	 
// 	   vgrad_test = Mult * (c0/(1 + V_value)) * (Vx_value*test100 + Vy_value*test010 + Vz_value*test001);
//A receptor kinetic law
//         vgrad_test = Mult * (c0/((1 + V_value ) * (1 + V_value))) * (Vx_value*test100 + Vy_value*test010 + Vz_value*test001); 
//Cooperative binding
//         vgrad_test = Mult * (c0 * 2* V_value  /((1 + V_value * V_value)*(1 + V_value * V_value))) * (Vx_value*test100 + Vy_value*test010 + Vz_value*test001); 
//Michaelis-Menten receptor kinetics	 
	 vgrad_test = Mult * c0 * (V_value /(1 + V_value)) * (Vx_value*test100 + Vy_value*test010 + Vz_value*test001); 
//=====================================================================================================================	   
	 
         uvfact = Mult * lambda * (1.-U_value-V_value) * test000;

         for(l=0;l<N_U_LocalDof;l++)
          { 
           ansatz000 = Orig[l];
	   ansatz100 = x_Orig[l];
	   ansatz010 = y_Orig[l];
	   ansatz001 = z_Orig[l];
	   
	   val1 = diff_fact * (ansatz100 * test100 + ansatz010 * test010 + ansatz001 * test001);
	   val2 = (vgrad_test + uvfact)*ansatz000;
	   
           LocMatrixA[k*N_U_LocalDof + l] +=  val1-val2;
          } //  for(l=0;l<N_U_LocalDof;l++
        } //  for(k=0;k<N_U_LocalDof;k++)
      } // for(j=0;j<N_Points;j++)
      
      //local to global matrix
     for(l=0;l<N_U_LocalDof;l++)
      {
       TestDOF = DOF[l];

       begin = RowPtr[TestDOF];
       end = RowPtr[TestDOF+1];

       for(n=begin;n<end;n++)
        {
         for(m=0;m<N_U_LocalDof;m++)
          {
           if(KCol[n] == DOF[m])
             {
              ValuesA[n] += LocMatrixA[l*N_U_LocalDof+m];
              break;
             }
	  }
         }
       } //  for(l=0;l<N_U_LocalDof;l++)
      } // for(i=0;i<N_Cells;i++)
    
//    cout << ii << " test assemble " << endl;

  } //  for(ii=0;ii<Levels; ii+
//      exit(0);
}

 void AssembleWArhs(int Levels, TSquareMatrix3D **A, TFEFunction3D **u_FeFunction, double **RHS)
 {
 int ii, i, j, k, l, m, n, begin, end, N_Cells, N_LocalUsedElements;
 int N_Points, N_U_LocalDof, N_V_LocalDof, TestDOF;
 int *N_BaseFunct, *BeginIndex, *GlobalNumbers, *RowPtr, *KCol, *DOF, *UBeginIndex, *UGlobalNumbers, *UDOF;
 
 double *weights, *xi, *eta, *zeta, alpha, beta, *ValuesA, *VRHS;
 double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
 double AbsDetjk[MaxN_QuadPoints_3D];
 double **uorig, *Orig, *x_Orig, *y_Orig, *z_Orig, U_value;
 double **vorig, *V_Orig, V_value;
 double Mult, *U_Sol;
 double test000, ansatz00;
 double ufact, alphaU, alphaUbeta;
 double LocMatrixA[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D], LocRhs[MaxN_BaseFunctions3D];

 TFESpace3D *fespace, *U_fespace;
 TCollection *Coll;
 TBaseCell *Me;
 FE3D FEId, U_FEId;
 TFE3D *ele;
 FE3D LocalUsedElements[2];
 BaseFunct3D BaseFunct, U_BaseFunct, *BaseFuncts;

 bool *SecondDer;
 
 SecondDer = new bool[2];
 SecondDer[0] = FALSE;
 SecondDer[1] = FALSE;
 
 BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
 N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();
 
//  //W solution spaces
for(ii=0;ii<Levels; ii++)
  {
    fespace = A[ii]->GetFESpace();
    Coll = fespace->GetCollection();
    N_Cells = Coll->GetN_Cells();

    BeginIndex = fespace->GetBeginIndex();
    GlobalNumbers = fespace->GetGlobalNumbers();
    
    RowPtr = A[ii]->GetRowPtr();
    KCol = A[ii]->GetKCol();
    ValuesA = A[ii]->GetEntries();
    
  //U solution spaces
    U_Sol = u_FeFunction[ii]->GetValues();
    U_fespace = u_FeFunction[ii]->GetFESpace3D();
    UBeginIndex = U_fespace->GetBeginIndex();
    UGlobalNumbers = U_fespace->GetGlobalNumbers();
    
    N_LocalUsedElements = 2;
    alpha = TDatabase::ParamDB->REACTOR_P5;
    beta = TDatabase::ParamDB->REACTOR_P6;
    
    for(i=0;i<N_Cells;i++)
    {
      Me = Coll->GetCell(i);
      
      FEId = fespace->GetFE3D(i, Me);
      U_FEId = U_fespace->GetFE3D(i, Me);
      
      LocalUsedElements[0]=FEId;
      LocalUsedElements[1]=U_FEId;
      
      TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
			      Coll, Me, SecondDer, 
			      N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);
      
      BaseFunct = BaseFuncts[FEId];
      N_V_LocalDof = N_BaseFunct[FEId];
      
      U_BaseFunct = BaseFuncts[U_FEId];
      N_U_LocalDof = N_BaseFunct[FEId];
      
      vorig = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
      uorig = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
      
      memset(LocMatrixA, 0, N_V_LocalDof*N_V_LocalDof*SizeOfDouble);
      memset(LocRhs, 0, N_V_LocalDof*SizeOfDouble);
      
      DOF = GlobalNumbers + BeginIndex[i];
      UDOF = UGlobalNumbers + UBeginIndex[i];
    
      for(j=0;j<N_Points;j++)
      {
	Orig = uorig[j];
	
	U_value = 0;
	
	for(l=0; l<N_U_LocalDof; l++)
	{
	  U_value += U_Sol[UDOF[l]] * Orig[l];
	}
	
	Mult = weights[j]*AbsDetjk[j];
	alphaU = Mult * alpha * U_value;
	alphaUbeta = Mult *(alpha* U_value + beta);
	
	//Assembling local matrix and rhs
	
	for(k=0;k<N_V_LocalDof;k++)
	{
	  test000 = Orig[k];
	  LocRhs[k] += alphaU*test000;
	  
	  ufact = alphaUbeta*test000;
	  
	  for(l=0;l<N_V_LocalDof;l++)
	  {
	    ansatz00 = Orig[l];
	    LocMatrixA[k*N_V_LocalDof + l] += ufact*ansatz00;
	    
	  }//for(l=0;l<N_V_LocalDof;l++)
	    
	}//for(k=0;k<N_V_LocalDof;k++)
	
    }//for(j=0;j<N_Points;j++)
    
    //local to global matrix
   
   for(l=0;l<N_V_LocalDof;l++)
   {
     TestDOF = DOF[l];
    // RHS[TestDOF] += LocRhs[l];
     VRHS=RHS[ii];
     VRHS[TestDOF]+=LocRhs[l];
     
     begin = RowPtr[TestDOF];
     end = RowPtr[TestDOF+1];
     
     for(n=begin;n<end;n++)
     {
       for(m=0;m<N_V_LocalDof;m++)
       {
	 if(KCol[n] == DOF[m])
	 {
	   ValuesA[n] += LocMatrixA[l*N_V_LocalDof+m];
	   
	   break;
	 } //if(KCol[n] == DOF[m])

       }//for(m=0;m<N_V_LocalDof;m++)

     }//for(n=begin;n<end;n++)

   } //for(l=0;l<N_V_LocalDof;l++)
   
    }//for(i=0;i<N_Cells;i++)
//  cout << ii << " test assemble for w " << endl;
   } // for(ii=0;ii<Levels; ii++)
 }
 
void AssembleVMat(int Levels, TSquareMatrix3D **A, TFEFunction3D **w_FeFunction, TFEFunction3D **u_FeFunction, TFEFunction3D **v_FeFunction)
{
  int ii, i, j, k, l, m, n, begin, end, N_Cells, N_LocalUsedElements;
  int N_Points, N_W_LocDof, N_V_LocDof, TestDOF, N_U_LocDof;
  int *N_BaseFunct, *BeginIndex, *GlobalNumbers, *RowPtr, *KCol, *DOF, *WBeginIndex, *WGlobalNumbers, *WDOF;
  int *UBeginIndex, *UGlobalNumbers, *UDOF;
  
  double *weights, *xi, *eta, *zeta, coeff, *ValuesA, *U_Sol, mu_a;  
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];;
  double AbsDetjk[MaxN_QuadPoints_3D];
  double **uorig, **worig, *Orig, *x_Orig, *y_Orig, W_value, U_value;
  double **vorig,  *V_Orig, V_value;
  double Mult, *W_Sol, *V_Sol;
  double test000, ansatz000;
  double wfact;
  double LocMatrixA[MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
  
  TFESpace3D *fespace, *W_fespace, *U_fespace;  
  TCollection *Coll;
  TBaseCell *Me;
  FE3D FEId, W_FEId, U_FEId;
  TFE3D *ele;
  FE3D LocalUsedElements[3];
  BaseFunct3D BaseFunct, W_BaseFunct, U_BaseFunct, *BaseFuncts;
   
  bool *SecondDer;
  
  SecondDer = new bool[2];
  SecondDer[0] = FALSE;
  SecondDer[1] = FALSE;
  
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();
  
  
    for(ii=0;ii<Levels;ii++)
      {
    //v solutions space
    V_Sol = v_FeFunction[ii]->GetValues();
    fespace = A[ii]->GetFESpace();
    Coll = fespace->GetCollection();
    N_Cells = Coll->GetN_Cells();
    
    BeginIndex = fespace->GetBeginIndex();
    GlobalNumbers = fespace->GetGlobalNumbers();
    
    RowPtr = A[ii]->GetRowPtr();
    KCol = A[ii]->GetKCol();
    ValuesA = A[ii]->GetEntries();
    
    //w solution space
     W_Sol = w_FeFunction[ii]->GetValues();
     W_fespace = w_FeFunction[ii]->GetFESpace3D();
     WBeginIndex = W_fespace->GetBeginIndex();
     WGlobalNumbers = W_fespace->GetGlobalNumbers();
     
     //u solution space
      U_Sol = u_FeFunction[ii]->GetValues();
      U_fespace = u_FeFunction[ii]->GetFESpace3D();
      UBeginIndex = U_fespace->GetBeginIndex();
      UGlobalNumbers = U_fespace->GetGlobalNumbers();
      
      N_LocalUsedElements = 3;
      coeff = TDatabase::ParamDB->REACTOR_P4;
      mu_a = TDatabase::ParamDB->REACTOR_P7;
    
	for(i=0;i<N_Cells;i++)
	{
	Me = Coll->GetCell(i);
	
 	FEId = fespace->GetFE3D(i, Me);
	W_FEId = W_fespace->GetFE3D(i, Me);
	U_FEId = U_fespace->GetFE3D(i, Me);
	
	LocalUsedElements[0] = FEId;
	LocalUsedElements[1] = W_FEId;
	LocalUsedElements[2] = U_FEId;
	
	TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
			       Coll, Me, SecondDer, 
			       N_Points, xi, eta, zeta, weights, X, Y, Z, AbsDetjk);
	
	BaseFunct = BaseFuncts[FEId];
	N_V_LocDof = N_BaseFunct[FEId];
	
	W_BaseFunct = BaseFuncts[W_FEId];
	N_W_LocDof  = N_BaseFunct[W_FEId];
	
	U_BaseFunct = BaseFuncts[U_FEId];
	N_U_LocDof = N_BaseFunct[U_FEId];
	
	vorig = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
	uorig = TFEDatabase3D::GetOrigElementValues(U_BaseFunct, D000);
	worig = TFEDatabase3D::GetOrigElementValues(W_BaseFunct, D000);
	
	memset(LocMatrixA, 0, N_V_LocDof*N_V_LocDof*SizeOfDouble);
	
	DOF = GlobalNumbers + BeginIndex[i];
	WDOF = WGlobalNumbers + WBeginIndex[i];
	UDOF = UGlobalNumbers + UBeginIndex[i];
	
	for(j=0;j<N_Points;j++)
	{
	  
	  //calculation for old w
	  Orig = worig[j];
	  W_value =0.;
	  for(l=0;l<N_W_LocDof;l++)
	  {
	    W_value += W_Sol[WDOF[l]] * Orig[l];
	  }
	  
	  //calculation for old u
	  Orig = uorig[j];
	  U_value = 0.;
	  for(l=0; l<N_U_LocDof;l++)
	  {
	    U_value += U_Sol[UDOF[l]]*Orig[l];
	  }
	  
	  //calculation for old v
	  Orig = vorig[j];
	  V_value = 0.;
	  for(l=0;l<N_V_LocDof;l++)
	  {
	    V_value += V_Sol[DOF[l]]*Orig[l];
	  }
	  
	  //weights and other terms
	  Mult = weights[j] * AbsDetjk[j] * (coeff * W_value - mu_a *(1.-U_value-V_value) );
	 
	  //assemble the local matrix
	  
	  for(k=0;k<N_V_LocDof;k++)
	    {
	      test000 = Orig[k];
	      wfact = Mult * test000;
	      
	      for(l=0;l<N_V_LocDof;l++)
	      {
		ansatz000 = Orig[l];
		LocMatrixA[k*N_V_LocDof + l] += wfact * ansatz000;
	      }//for(l=0;l<N_V_LocDof;l++)
	      
	    }//for{k=0;k<N_V_LocDof;k++)
	      
	      
	}//(j=0;j<N_Points;j++)
	
	//add local to global
	for(l=0;l<N_V_LocDof;l++)
	{
	  TestDOF = DOF[l];
	  begin = RowPtr[TestDOF];
	  end = RowPtr[TestDOF+1];
	  
	  for(n=begin;n<end;n++)
	  {
	    for(m=0;m<N_V_LocDof;m++)
	    {
	      if(KCol[n]==DOF[m])
	      {
		ValuesA[n] += LocMatrixA[l*N_V_LocDof+m];
		break;
	      }
	    }//for(m=0;m<N_V_LocDof;m++)
	  }//for(n=begin;n<end;n++)
	}//for(l=0;l<N_V_LocDof;l++)
	
	
      }//for(i=0;i<N_Cells;i++)
   }//for(ii=0;ii<Levels;ii++)
    
}//assemblevmat



int main(int argc, char* argv[])
{
// ======================================================================
// variable declaration
// ======================================================================
  int i, j, l, m, N_SubSteps, ORDER, LEVELS, mg_level, N_Cells, N_DOF, N_V_DOF, N_W_DOF, img=1;
  int N_Active, N_V_Active, N_W_Active, Max_It_scalar, Min_It_scalar, mg_type;
  
  double *sol, *oldsol, *rhs, *oldrhs, t1, t2, errors[5],  v_errors[5], w_errors[5], **Sol_array, **Rhs_array, end_time;
  double tau, *defect, olderror, v_olderror, w_olderror, olderror1, v_olderror1, w_olderror1, hmin, hmax, temp, reduced_errors[4];
  double *V_sol, *V_oldsol, *V_rhs, *V_oldrhs, *V_defect, **VSol_array, **VRhs_array;
  double *W_sol, *W_oldsol, *W_rhs, *W_oldrhs, *W_defect, **WSol_array, **WRhs_array;
  double residual_scalar, oldresidual_scalar, limit_scalar;
  double Parameters[10], *oldrhs_itr, *V_oldrhs_itr, *W_oldrhs_itr;
  
  double start_time, stop_time, 
	 start_vtk=0, end_vtk=0, total_vtk=0,
	 start_assembling=0, end_assembling=0, total_assembling=0,
	 start_solve=0, end_solve=0, total_solve=0,
	 start_int=0, end_int=0, total_int=0;
	 
  bool UpdateStiffnessMat, UpdateRhs,  ConvectionFirstTime;
  
  TDatabase *Database = new TDatabase();
  
  const char vtkdir[] = "VTK"; 
  char *PsBaseName, *VtkBaseName, *GEO, *PRM;
  char Name[] = "name";
  char Description[] = "description";
  char CString[] = "C";
  char VString[] = "v";
  char WString[] = "w";
  double Linfty=0, v_Linfty=0., w_Linfty=0.;
  char SubID[] = "";
  
  int profiling;
#ifdef _MPI
  int rank, size, out_rank;
  int MaxCpV, MaxSubDomainPerDof;
  
  MPI_Comm Comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  double time1, time2;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
 
  TDatabase::ParamDB->Comm = Comm;
  int Refine;
  int metisType[2] = {0,0};  
#endif
  
  TDomain *Domain;
  TFEDatabase3D *FEDatabase = new TFEDatabase3D(); 
  TCollection *coll;
  TFESpace3D **Scalar_FeSpaces, **Scalar_V_FeSpaces, **Scalar_W_FeSpaces, *fesp[1];
  TFEFunction3D *Scalar_FeFunction, *Scalar_V_FeFunction, *Scalar_W_FeFunction, 
                **Scalar_FeFunctions, **Scalar_V_FeFunctions, **Scalar_W_FeFunctions;
  TOutput3D *Output;
  TSystemTCD3D *SystemMatrix, *V_SystemMatrix, *W_SystemMatrix;
  TAuxParam3D *aux;
  MultiIndex3D AllDerivatives[4] = {D000, D100, D010, D001};
  TSquareMatrix3D  **V_A_Matrix;
  
//   TSquareMatrix3D *sqmatrixA, *V_sqmatrixA, *W_sqmatrixA, **sqmatricesA, *SQMATRICES[3];;
//   TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  
  std::ostringstream os;
  os << " ";   
    
  mkdir(vtkdir, 0777);

// ======================================================================
// set the database values and generate mesh
// ======================================================================    
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  Domain = new TDomain(argv[1]);  
  
  profiling = TDatabase::ParamDB->timeprofiling;
  if(profiling)
  {
#ifdef _MPI
    start_time = MPI_Wtime();
#else
    start_time = GetTime();
#endif
  }
  
  OpenFiles();
  OutFile.setf(std::ios::scientific);

  #ifdef _MPI
  out_rank=TDatabase::ParamDB->Par_P0;
  //out_rank = 0;
  if(rank == out_rank)
  #endif
   {
    Database->WriteParamDB(argv[0]);
    Database->WriteTimeDB();
    ExampleFile();
   }
  
  /* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
  // standard mesh  
  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
  PsBaseName = TDatabase::ParamDB->BASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  Domain->Init(PRM, GEO);

// ================================================================================================
                            // Creating Realsitc Geometry 
// ================================================================================================

/** with atlas mesh, no tetgen*/
    TetrameshCreate(Domain);
    
/** Using tetgen with smesh mesh */
  // TetrameshGen(Domain);
   
/** same as TetrameshGen but without using the face list info from tetgen */
//     TetraGen(Domain);

/** input includes tetrahedron */
   //TetraGenWithInputCells(Domain);
   
   
//    Domain->ConvertQuadToTri(TDatabase::ParamDB->CONVERT_QUAD_TO_TRI);
  
// ================================================================================================  
  
//   refine grid up to the coarsest level
//   for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
//     Domain->RegRefineAll();
  
  #ifdef _MPI
  Domain->GenerateEdgeInfo();
  
  if(profiling)  t1 = MPI_Wtime();
  
  if(rank == 0)
       {
	  printf("\n----------------------------------------------------------------------------------------\n");
	  printf("metis type set to %d\n",TDatabase::ParamDB->Par_P2);
	  printf("----------------------------------------------------------------------------------------\n\n");
       }
  
  do
  {
    metisType[TDatabase::ParamDB->Par_P2] = 1;
    Refine = Partition_Mesh3D(Comm, Domain, MaxCpV);	//MaxCpV=maximum cell per vertex
    
    if(metisType[0]*metisType[1] == 1 && Refine)
    {
      metisType[0] = 0;      metisType[1] = 0;
      TDatabase::ParamDB->Par_P2 = 0;
      if(rank == 0)
       {
	  printf("\n----------------------------------------------------------------------------------------\n");
	  printf("Warning :: both metisType used. Now refining the mesh by one step \n");
	  printf("metis type set to 0\n");
	  printf("----------------------------------------------------------------------------------------\n\n");
       }
      Domain->RegRefineAll();
      Domain->GenerateEdgeInfo();
      TDatabase::ParamDB->UNIFORM_STEPS +=1;
    }
  }while(Refine);
  
  if(profiling)  t2 = MPI_Wtime(); 
  
  if(profiling){
    time2 = t2-t1;
    MPI_Reduce(&time2, &time1, 1, MPI_DOUBLE, MPI_MAX, out_rank, Comm);
    if(rank == out_rank)
      printf("Time taken for Domain Decomposition is %e\n", time1);
  }

  Domain->GenerateEdgeInfo();
  MaxSubDomainPerDof = MIN(MaxCpV, size);
  TDatabase::ParamDB->WRITE_PS = 0;
  #endif   
  
  if(TDatabase::ParamDB->WRITE_PS)
   {
    // write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
   }
 
//=========================================================================
// set data for multigrid
//=========================================================================  
  LEVELS = TDatabase::ParamDB->LEVELS;
  LEVELS = 2;		//we will use order -1, 1, 2 for developing the hierarchy
  TDatabase::ParamDB->LEVELS = 2;
  int levelOrder[3];
  levelOrder[0]=-1;
  levelOrder[1]=1;
  levelOrder[2]=2;
  // set type of multilevel
  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SCALAR;
 
  if(TDatabase::ParamDB->SOLVER_TYPE==AMG_SOLVE|| TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
   { 
     mg_type=0; 
     TDatabase::ParamDB->SC_MG_TYPE_SCALAR = mg_type;
    }
  
//   if(mg_type)
//    {
//     mg_level =  LEVELS + 1;
//     ORDER = -1;
//    }
//   else
   {
    mg_level = LEVELS;
    ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
   }
   
  if(TDatabase::ParamDB->SOLVER_TYPE==GMG)
   {
#ifdef _MPI  
    if(rank == out_rank)
    {
#endif 
    OutPut("=======================================================" << endl);
    OutPut("======           GEOMETRY  LEVEL ");
    OutPut(LEVELS << "              ======" << endl);
    OutPut("======           MULTIGRID LEVEL ");
    OutPut(mg_level << "              ======" << endl);
    OutPut("=======================================================" << endl);   
#ifdef _MPI 
    }
#endif 
   }
    
  Scalar_FeSpaces = new TFESpace3D*[mg_level];  
  Scalar_V_FeSpaces = new TFESpace3D*[mg_level]; 
  Scalar_W_FeSpaces = new TFESpace3D*[mg_level]; 
  
  Scalar_FeFunctions = new TFEFunction3D*[mg_level]; 
  Scalar_V_FeFunctions = new TFEFunction3D*[mg_level];
  Scalar_W_FeFunctions = new TFEFunction3D*[mg_level];
  
  Sol_array = new double*[mg_level];
  Rhs_array = new double*[mg_level];
  VSol_array = new double*[mg_level];
  VRhs_array = new double*[mg_level];
  WSol_array = new double*[mg_level];
  WRhs_array = new double*[mg_level];
 
//=========================================================================
// construct all finite element spaces
// loop over all levels (not a multigrid level but for convergence study)  
//=========================================================================
  for(i=0;i<LEVELS;i++)
   {   
//     if(i)
//      { Domain->RegRefineAll(); }
     
//      #ifdef _MPI
//      if(rank == out_rank)
//        printf("Level :: %d\n\n",i);
//      if(i)
//      {
//        Domain->GenerateEdgeInfo();
//        Domain_Crop(Comm, Domain);       // removing unwanted cells in the hallo after refinement 
//      }
//      #endif

     ORDER = levelOrder[i];

     coll=Domain->GetCollection(It_Finest, 0);
     //continue;
     // fespaces for scalar equation 
     Scalar_FeSpaces[i] =  new TFESpace3D(coll, Name, Description, BoundCondition, ORDER);     
     Scalar_V_FeSpaces[i] =  new TFESpace3D(coll, Name, Description, V_BoundCondition, ORDER); 
     Scalar_W_FeSpaces[i] =  new TFESpace3D(coll, Name, Description, W_BoundCondition, ORDER);    

     #ifdef _MPI
     Scalar_FeSpaces[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
     Scalar_V_FeSpaces[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
     Scalar_W_FeSpaces[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
     #endif
     
//      //multilevel multigrid disc
//      if(i==LEVELS-1 && mg_type==1) 
//       {
//        ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
//        Scalar_FeSpaces[mg_level-1] =  new TFESpace3D(coll, Name, Description, BoundCondition, ORDER);
//        Scalar_V_FeSpaces[mg_level-1] =  new TFESpace3D(coll, Name, Description, V_BoundCondition, ORDER);
//        Scalar_W_FeSpaces[mg_level-1] =  new TFESpace3D(coll, Name, Description, W_BoundCondition, ORDER);
//        
//        #ifdef _MPI
//        Scalar_FeSpaces[mg_level-1]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
//        Scalar_V_FeSpaces[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
//        Scalar_W_FeSpaces[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
//        #endif
//       } //  if(i==LEVELS-1 && i!=mg_level-1) 
     

     
//======================================================================
// construct all finite element functions
//======================================================================
    //===================================================================
    //u solution
    //===================================================================
    N_DOF = Scalar_FeSpaces[i]->GetN_DegreesOfFreedom();

    sol = new double[N_DOF];
    rhs = new double[N_DOF];
    Sol_array[i] = sol;
    Rhs_array[i] = rhs;   
    
    Scalar_FeFunction  = new TFEFunction3D(Scalar_FeSpaces[i], CString, CString, sol, N_DOF);  
    Scalar_FeFunctions[i] = Scalar_FeFunction;
     
//     if(i==LEVELS-1 && mg_type==1) 
//      {  
//       N_DOF = Scalar_FeSpaces[mg_level-1]->GetN_DegreesOfFreedom();
//       sol = new double[N_DOF];
//       rhs = new double[N_DOF];
//       Sol_array[mg_level-1] = sol;
//       Rhs_array[mg_level-1] = rhs;
// 
//       Scalar_FeFunction = new TFEFunction3D(Scalar_FeSpaces[mg_level-1], CString, CString, sol, N_DOF);   
//       Scalar_FeFunctions[mg_level-1] = Scalar_FeFunction;
//      }//   if(i==LEVELS-1 && mg_type==1) 
    
   //========================================================================
   //v solution
   //========================================================================
    N_V_DOF = Scalar_V_FeSpaces[i]->GetN_DegreesOfFreedom();
    V_sol = new double[N_V_DOF];
    V_rhs = new double[N_V_DOF];
    VSol_array[i] = V_sol;
    VRhs_array[i] = V_rhs;   
    
    Scalar_V_FeFunction  = new TFEFunction3D(Scalar_V_FeSpaces[i], VString, VString, V_sol, N_V_DOF);  
    Scalar_V_FeFunctions[i] = Scalar_V_FeFunction;
     
//     if(i==LEVELS-1 && mg_type==1) 
//      {
//       N_V_DOF = Scalar_V_FeSpaces[mg_level-1]->GetN_DegreesOfFreedom();
//       V_sol = new double[N_V_DOF];
//       V_rhs = new double[N_V_DOF];
//       VSol_array[mg_level-1] = V_sol;
//       VRhs_array[mg_level-1] = V_rhs;
// 
//       Scalar_V_FeFunction = new TFEFunction3D(Scalar_V_FeSpaces[mg_level-1], VString, VString, V_sol, N_V_DOF);   
//       Scalar_V_FeFunctions[mg_level-1] = Scalar_V_FeFunction;
//      }//   if(i==LEVELS-1 && mg_type==1) 
  //  =============================================================================
    
    
   //  ========================================================================
   // w solution
   //========================================================================
    N_W_DOF = Scalar_W_FeSpaces[i]->GetN_DegreesOfFreedom();
    W_sol = new double[N_W_DOF];
    W_rhs = new double[N_W_DOF];
    WSol_array[i] = W_sol;
    WRhs_array[i] = W_rhs;   
    
    Scalar_W_FeFunction  = new TFEFunction3D(Scalar_W_FeSpaces[i], WString, WString, W_sol, N_W_DOF);  
    Scalar_W_FeFunctions[i] = Scalar_W_FeFunction;
     
//     if(i==LEVELS-1 && mg_type==1) 
//      {
//       N_W_DOF = Scalar_W_FeSpaces[mg_level-1]->GetN_DegreesOfFreedom();
//       W_sol = new double[N_W_DOF];
//       W_rhs = new double[N_W_DOF];
//       WSol_array[mg_level-1] = W_sol;
//       WRhs_array[mg_level-1] = W_rhs;
// 
//       Scalar_W_FeFunction = new TFEFunction3D(Scalar_W_FeSpaces[mg_level-1], WString, WString, W_sol, N_W_DOF);   
//       Scalar_W_FeFunctions[mg_level-1] = Scalar_W_FeFunction;
//      }//    // if(i==LEVELS-1 && mg_type==1) 
     //========================================================================
     
   }// for(i=0;i<LEVELS;i++)    
   

   oldrhs = new double[N_DOF];
   oldsol = new double[N_DOF];
   oldrhs_itr = new double[N_DOF];
   
   V_oldrhs = new double[N_V_DOF];
   V_oldrhs_itr = new double[N_V_DOF];
   V_oldsol = new double[N_V_DOF];
   
   W_oldrhs = new double[N_W_DOF];
   W_oldrhs_itr = new double[N_W_DOF];   
   W_oldsol = new double[N_W_DOF];
   
#ifdef _MPI
   N_Cells = coll->GetN_Cells();
   printf("rank=%d\t N_Cells   : %d\t Dof all (U)  :%d\t Dof all (V)  :%d\t Dof all (W)  :%d\n",
	  rank,N_Cells,N_DOF,N_V_DOF,N_W_DOF);
#endif

#ifndef _MPI   
   N_Cells = coll->GetN_Cells();
   OutPut("N_Cells      : " << N_Cells <<endl);
   OutPut("Dof all (U)  : " << N_DOF  << endl);
   OutPut("Dof all (V)  : " << N_V_DOF  << endl);
   OutPut("Dof all (W)  : " << N_W_DOF  << endl);
   OutPut(endl);
#endif

//======================================================================
// SystemMatrix construction and solution
//======================================================================  
    /** interpolate the initial value */
    Scalar_FeFunction->Interpolate(InitialCondition);   
    Scalar_V_FeFunction->Interpolate(V_InitialCondition); 
    Scalar_W_FeFunction->Interpolate(W_InitialCondition);   
    
    /** Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
     *  Solver: AMG_SOLVE (or) GMG  (or) DIRECT */
    if(profiling){
#ifdef _MPI
      start_int = MPI_Wtime();
#else
      start_int = GetTime();
#endif
    }
 
   //u sysytem
   SystemMatrix = new TSystemTCD3D(mg_level, Scalar_FeSpaces, TDatabase::ParamDB->DISCTYPE, 
					      TDatabase::ParamDB->SOLVER_TYPE);
  
/*   cout<<"11111111111111111111111yes"<<endl;
   exit(0);  */  
    /** initilize the system matrix with the functions defined in Example file */
    SystemMatrix->Init(BilinearCoeffs, BoundCondition, BoundValue);

    //v sysytem
    V_SystemMatrix = new TSystemTCD3D(mg_level, Scalar_V_FeSpaces, TDatabase::ParamDB->DISCTYPE, 
						TDatabase::ParamDB->SOLVER_TYPE);
    V_SystemMatrix->Init(V_BilinearCoeffs, V_BoundCondition, V_BoundValue);
  
    //w sysytem
    W_SystemMatrix = new TSystemTCD3D(mg_level, Scalar_W_FeSpaces, TDatabase::ParamDB->DISCTYPE, 
 						TDatabase::ParamDB->SOLVER_TYPE);
    W_SystemMatrix->Init(W_BilinearCoeffs, W_BoundCondition, W_BoundValue);
   
    
    if(profiling){
#ifdef _MPI
      end_int = MPI_Wtime();
#else
      end_int = GetTime();
#endif
      total_int = end_int-start_int;
    }
//exit(0);
    /** assemble the system matrix with given aux, sol and rhs 
     *  aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
     *  otherwise, just pass it with NULL  */
    if(profiling){
#ifdef _MPI
      start_assembling = MPI_Wtime();
#else
      start_assembling = GetTime();
#endif
    }
    
    SystemMatrix->AssembleMRhs(NULL, Sol_array, Rhs_array); 
    V_SystemMatrix->AssembleMRhs(NULL, VSol_array, VRhs_array);
    W_SystemMatrix->AssembleMRhs(NULL, WSol_array, WRhs_array);
    
    if(profiling){
#ifdef _MPI
      end_assembling = MPI_Wtime();
#else
      end_assembling = GetTime();
#endif
    }
  //exit(0);  
   /** copy rhs to oldrhs before calling the solver, as rhs will change in multigrid solver */
    memcpy(oldrhs, rhs, N_DOF*SizeOfDouble);
    memcpy(V_oldrhs, V_rhs, N_V_DOF*SizeOfDouble);
    memcpy(W_oldrhs, W_rhs, N_W_DOF*SizeOfDouble);
        
//======================================================================
// produce outout at t=0
//======================================================================   
    Output = new TOutput3D(2, 2, 1, 1, Domain);

    Output->AddFEFunction(Scalar_FeFunction);
    Output->AddFEFunction(Scalar_V_FeFunction);
    Output->AddFEFunction(Scalar_W_FeFunction);

//     Scalar_FeFunction->Interpolate(Exact);
#ifdef _MPI
    if(profiling)	start_vtk = MPI_Wtime();

    if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
        img++;

    if(profiling)	end_vtk = MPI_Wtime();
#else
    if(profiling)	start_vtk = GetTime();
   
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
    if(profiling)	end_vtk = GetTime();
#endif

// MPI_Finalize();
// exit(0);    
//====================================================================== 
    /** measure errors to known solution u */
//======================================================================     
    if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
      fesp[0] = Scalar_FeSpaces[mg_level-1];       
      aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      
      for(j=0;j<5;j++)
       errors[j] = 0;

      Scalar_FeFunction->GetErrors(Exact, 4, AllDerivatives, 2, L2H1Errors,
                                   BilinearCoeffs, aux, 1, fesp, errors);
    
#ifdef _MPI
      MPI_Allreduce(errors, reduced_errors, 2, MPI_DOUBLE, MPI_SUM, Comm);
      for(i=0;i<2;i++)
	errors[i] = sqrt(reduced_errors[i]);
      if(rank == out_rank){
	OutPut(endl);
	OutPut("time: " << TDatabase::TimeDB->CURRENTTIME<<endl);
	OutPut( "L2 (u): " << sqrt(reduced_errors[0]) << endl);
	OutPut( "H1-semi (u): " << sqrt(reduced_errors[1]) << endl);
      }
#else
      OutPut(endl);
      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME<<endl);
      OutPut( "L2 (u): " << errors[0] << endl);
      OutPut( "H1-semi (u): " << errors[1] << endl);
#endif           
      Linfty=errors[0];
    } //  if(TDatabase::ParamDB->MEASURE_ERRORS)
    
//====================================================================== 
    /** measure errors to known solution v */
//======================================================================     
    if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
      fesp[0] = Scalar_V_FeSpaces[mg_level-1];       
      aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      
      for(j=0;j<5;j++)
       v_errors[j] = 0;

      Scalar_V_FeFunction->GetErrors(V_Exact, 4, AllDerivatives, 2, L2H1Errors,
                                   V_BilinearCoeffs, aux, 1, fesp, v_errors);
    
#ifdef _MPI
      MPI_Allreduce(v_errors, reduced_errors, 2, MPI_DOUBLE, MPI_SUM, Comm);
      for(i=0;i<2;i++)
	v_errors[i] = sqrt(reduced_errors[i]);
      if(rank == out_rank){
	OutPut(endl);
	OutPut("time: " << TDatabase::TimeDB->CURRENTTIME<<endl);
	OutPut( "L2 (v): " << sqrt(reduced_errors[0]) << endl);
	OutPut( "H1-semi (v): " << sqrt(reduced_errors[1]) << endl);
      }
#else
      OutPut(endl);
      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME<<endl);
      OutPut( "L2 (v): " << v_errors[0] << endl);
      OutPut( "H1-semi (v): " << v_errors[1] << endl);
#endif           
      v_Linfty=v_errors[0];
      OutPut( "v-Linfty " << v_Linfty << endl);
    } //  if(TDatabase::ParamDB->MEASURE_ERRORS)  

//====================================================================== 
    /** measure errors to known solution w */
//======================================================================     
    if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
      fesp[0] = Scalar_W_FeSpaces[mg_level-1];       
      aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      
      for(j=0;j<5;j++)
       w_errors[j] = 0;

      Scalar_W_FeFunction->GetErrors(W_Exact, 4, AllDerivatives, 2, L2H1Errors,
                                   W_BilinearCoeffs, aux, 1, fesp, w_errors);
    
#ifdef _MPI
      MPI_Allreduce(w_errors, reduced_errors, 2, MPI_DOUBLE, MPI_SUM, Comm);
      for(i=0;i<2;i++)
	w_errors[i] = sqrt(reduced_errors[i]);
      if(rank == out_rank){
	OutPut(endl);
	OutPut("time: " << TDatabase::TimeDB->CURRENTTIME<<endl);
	OutPut( "L2 (w): " << sqrt(reduced_errors[0]) << endl);
	OutPut( "H1-semi (w): " << sqrt(reduced_errors[1]) << endl);
      }
#else
      OutPut(endl);
      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME<<endl);
      OutPut( "L2 (w): " << w_errors[0] << endl);
      OutPut( "H1-semi (w): " << w_errors[1] << endl);
#endif           
      w_Linfty=w_errors[0];
      OutPut( "w-Linfty " << w_Linfty << endl);
    } //  if(TDatabase::ParamDB->MEASURE_ERRORS)        

//======================================================================
// time disc Parameters
//======================================================================       
#ifdef _MPI
     if(rank == out_rank)
#endif     
     cout << "time " << TDatabase::TimeDB->CURRENTTIME << endl;
     
     coll->GetHminHmax(&hmin,&hmax);
     
#ifdef _MPI
     if(rank == out_rank)
#endif
      OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
     
#ifdef _MPI
     MPI_Allreduce(&hmax, &temp, 1, MPI_DOUBLE, MPI_MAX, Comm);
     hmax = temp;
     temp = 0;
#endif
   
//    if(TDatabase::TimeDB->TIME_DISC == 0 || TDatabase::TimeDB->TIME_DISC == 1)
//      TDatabase::TimeDB->TIMESTEPLENGTH = hmax * hmax;
//    else if(TDatabase::TimeDB->TIME_DISC == 2)
//      TDatabase::TimeDB->TIMESTEPLENGTH = hmax;
   
   // parameters for time stepping scheme
   m = 0;
   N_SubSteps = GetN_SubSteps();
   end_time = TDatabase::TimeDB->ENDTIME; 
   
   Max_It_scalar = (int)TDatabase::ParamDB->REACTOR_P10;
   Min_It_scalar = TDatabase::ParamDB->REACTOR_P11;
   limit_scalar = TDatabase::ParamDB->REACTOR_P12;

   UpdateStiffnessMat = FALSE; //check BilinearCoeffs in example file
   UpdateRhs = TRUE; //check BilinearCoeffs in example file
   ConvectionFirstTime=TRUE;

//    cout << "init " << Ddot(N_DOF, sol, sol)<< endl;
//======================================================================
// time disc loop
//======================================================================      
   /** time loop starts */
   
   while(TDatabase::TimeDB->CURRENTTIME< end_time)
    {
     m++;
     TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

     for(l=0;l<N_SubSteps;l++) // sub steps of fractional step theta
      {
       SetTimeDiscParameters(1);

#ifdef _MPI
      if(rank == out_rank)
#endif
      if(m==1)
       {
        OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
        OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
        OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
        OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
       }

      /** copy the sol to old sol */  
      memcpy(oldsol, sol, N_DOF*SizeOfDouble);  
      memcpy(V_oldsol, V_sol, N_V_DOF*SizeOfDouble);
      memcpy(W_oldsol, W_sol, N_W_DOF*SizeOfDouble);   
       
      tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;
/*    OutPut(endl << "CURRENT TIME: "<<TDatabase::TimeDB->CURRENTTIME << endl);
    exit(0);   */   
#ifdef _MPI
      if(rank == out_rank)
#endif
      OutPut(endl << "CURRENT TIME: "<<TDatabase::TimeDB->CURRENTTIME << endl);
      
      for(j=0;j<Max_It_scalar;j++)
      {
//===========================================================================================
// u system
//===========================================================================================		
      /** unless the stiffness matrix or rhs change in time, it is enough to assemble only once at the begning   */
      
      if(profiling){
	total_assembling += (end_assembling-start_assembling); 
#ifdef _MPI
	start_assembling = MPI_Wtime();
#else
	start_assembling = GetTime();
#endif
      }
      
      if(UpdateStiffnessMat || UpdateRhs ||  ConvectionFirstTime)
       {  
        if(UpdateRhs)
         { SystemMatrix->AssembleARhs(NULL, Sol_array, Rhs_array); }
        else
         { SystemMatrix->AssembleARhs(NULL, Sol_array, Rhs_array); }
         
        AssembleCancerDensity(mg_level, SystemMatrix->GetAMatrix(), Scalar_FeFunctions, Scalar_V_FeFunctions,Scalar_W_FeFunctions);
         
        /**  M:= M + (tau*TDatabase::TimeDB->THETA1)*A
         *   rhs: =(tau*THETA4)*rhs +(tau*THETA3)*oldrhs + [ M - (tau*THETA2)A]*oldsol **/
	SystemMatrix->AssembleSystMat(oldrhs, oldsol, rhs, sol
#ifdef _MPI
					    , Rhs_array
#endif
							);
        /** copy rhs to oldrhs before calling the solver, as rhs will change in multigrid solver */
        memcpy(oldrhs_itr, rhs, N_DOF*SizeOfDouble); 
       }
      
      if(profiling){
#ifdef _MPI
	end_assembling = MPI_Wtime();
#else
	end_assembling = GetTime();
#endif
      }
      
      residual_scalar = SystemMatrix->GetResidual(sol);
#ifdef _MPI
      if(rank == out_rank)
#endif
      {
        OutPut("Scalar nonlinear step " << setw(3) << j);
        OutPut(setw(14) << residual_scalar); // sqrt of residual_scalar is alread done in ScalarDefect
      }
      
#ifdef _MPI
	  if(rank == out_rank)
#endif
      if (j>0)
      { 
          if(fabs(oldresidual_scalar)>0)
          OutPut(setw(14) <<  residual_scalar/oldresidual_scalar ); 
      }
      
#ifdef _MPI
     if(rank == out_rank)
#endif      
      OutPut(endl);
//       cout<<"yes"<<endl;
      oldresidual_scalar = residual_scalar;
   
      if( ((residual_scalar<=limit_scalar)||(j==Max_It_scalar-1))  && (j>=Min_It_scalar) )
      {
	if(UpdateStiffnessMat || UpdateRhs)        
            SystemMatrix->RestoreMassMat(); 
        break;
      } 
         
      // solve the system matrix 
      if(profiling){
	total_solve += (end_solve-start_solve);
#ifdef _MPI
	start_solve = MPI_Wtime();
#else
	start_solve = GetTime();
#endif
      }

      SystemMatrix->Solve(sol);
      
      if(profiling){
#ifdef _MPI
	end_solve = MPI_Wtime();
#else
	end_solve = GetTime();
#endif
      }
      
      /** restore the mass matrix for the next time step unless the stiffness matrix 
       * or rhs change in time, it is not necessary to assemble the system matrix in every time step */
      if(UpdateStiffnessMat || UpdateRhs)
       {         
        SystemMatrix->RestoreMassMat();
       }
//===========================================================================================
// w system
//===========================================================================================		
      /** unless the stiffness matrix or rhs change in time, it is enough to assemble only once at the begning   */
      
      if(profiling){
	total_assembling += (end_assembling-start_assembling); 
#ifdef _MPI
	start_assembling = MPI_Wtime();
#else
	start_assembling = GetTime();
#endif
      }
      
      if(UpdateStiffnessMat || UpdateRhs ||  ConvectionFirstTime)
       {  
        if(UpdateRhs)
         { W_SystemMatrix->AssembleARhs(NULL, WSol_array, WRhs_array); }
        else
         { W_SystemMatrix->AssembleARhs(NULL, WSol_array, WRhs_array); }
        
        //###################################################doubt on this part ??????
        AssembleWArhs(mg_level, W_SystemMatrix->GetAMatrix(), Scalar_FeFunctions, WRhs_array);
	//###################################################doubt on this part ??????
         
        /**  M:= M + (tau*TDatabase::TimeDB->THETA1)*A
         *   rhs: =(tau*THETA4)*rhs +(tau*THETA3)*oldrhs + [ M - (tau*THETA2)A]*oldsol **/
	W_SystemMatrix->AssembleSystMat(W_oldrhs, W_oldsol, W_rhs, W_sol
#ifdef _MPI
					    , WRhs_array
#endif
							);
        /** copy rhs to oldrhs before calling the solver, as rhs will change in multigrid solver */
        memcpy(W_oldrhs_itr, W_rhs, N_W_DOF*SizeOfDouble); 
       }
      
      if(profiling){
#ifdef _MPI
	end_assembling = MPI_Wtime();
#else
	end_assembling = GetTime();
#endif
      }
         
      // solve the system matrix 
      if(profiling){
	total_solve += (end_solve-start_solve);
#ifdef _MPI
	start_solve = MPI_Wtime();
#else
	start_solve = GetTime();
#endif
      }

      W_SystemMatrix->Solve(W_sol);
      
      if(profiling){
#ifdef _MPI
	end_solve = MPI_Wtime();
#else
	end_solve = GetTime();
#endif
      }
      
      /** restore the mass matrix for the next time step unless the stiffness matrix 
       * or rhs change in time, it is not necessary to assemble the system matrix in every time step */
      if(UpdateStiffnessMat || UpdateRhs)
       {         
        W_SystemMatrix->RestoreMassMat();
       }
       
//===========================================================================================
// V system
//===========================================================================================		
      /** unless the stiffness matrix or rhs change in time, it is enough to assemble only once at the begning   */
      
      if(profiling){
	total_assembling += (end_assembling-start_assembling); 
#ifdef _MPI
	start_assembling = MPI_Wtime();
#else
	start_assembling = GetTime();
#endif
      }
      
      if(UpdateStiffnessMat || UpdateRhs ||  ConvectionFirstTime)
       {  
        if(UpdateRhs)
         { V_SystemMatrix->AssembleARhs(NULL, VSol_array, VRhs_array); }
        else
         { V_SystemMatrix->AssembleARhs(NULL, VSol_array, VRhs_array); }
        
        //###################################################doubt on this part ??????
        AssembleVMat(mg_level, V_SystemMatrix->GetAMatrix(), Scalar_W_FeFunctions, Scalar_FeFunctions, Scalar_V_FeFunctions);
	//###################################################doubt on this part ??????
         
        /**  M:= M + (tau*TDatabase::TimeDB->THETA1)*A
         *   rhs: =(tau*THETA4)*rhs +(tau*THETA3)*oldrhs + [ M - (tau*THETA2)A]*oldsol **/
	V_SystemMatrix->AssembleSystMat(V_oldrhs, V_oldsol, V_rhs, V_sol
#ifdef _MPI
					    , VRhs_array
#endif
							);
        /** copy rhs to oldrhs before calling the solver, as rhs will change in multigrid solver */
        memcpy(V_oldrhs_itr, V_rhs, N_V_DOF*SizeOfDouble); 
	
	ConvectionFirstTime = FALSE;
       }
      
      if(profiling){
#ifdef _MPI
	end_assembling = MPI_Wtime();
#else
	end_assembling = GetTime();
#endif
      }
         
      // solve the system matrix 
      if(profiling){
	total_solve += (end_solve-start_solve);
#ifdef _MPI
	start_solve = MPI_Wtime();
#else
	start_solve = GetTime();
#endif
      }

      V_SystemMatrix->Solve(V_sol);
      
      if(profiling){
#ifdef _MPI
	end_solve = MPI_Wtime();
#else
	end_solve = GetTime();
#endif
      }
      
      /** restore the mass matrix for the next time step unless the stiffness matrix 
       * or rhs change in time, it is not necessary to assemble the system matrix in every time step */
      if(UpdateStiffnessMat || UpdateRhs)
       {         
        V_SystemMatrix->RestoreMassMat();
       }
     } // for(j=o;j<Max_It_scalar;j++)
//=========================================================================================================
     if(UpdateStiffnessMat || UpdateRhs ||  ConvectionFirstTime)
       {
        memcpy(oldrhs, oldrhs_itr, N_DOF*SizeOfDouble); 
        memcpy(W_oldrhs, W_oldrhs_itr, N_W_DOF*SizeOfDouble); 
        memcpy(V_oldrhs, V_oldrhs_itr, N_V_DOF*SizeOfDouble);

       }
   } // for(l=0;l<N_SubSteps;l++) 
//======================================================================
// produce output
//======================================================================
#ifdef _MPI
    if(profiling){
      total_vtk += (end_vtk-start_vtk);
      start_vtk = MPI_Wtime();
    }
  if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
    if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
        img++;

    if(profiling)	end_vtk = MPI_Wtime();
#else
    if(profiling){
      total_vtk += (end_vtk-start_vtk);
      start_vtk = GetTime();
    }
   
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
    if(profiling)	end_vtk = GetTime();
#endif

//======================================================================
// measure errors to known solution u
//======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {      
      olderror  = errors[0];
      olderror1 = errors[1];
 
      Scalar_FeFunction->GetErrors(Exact, 4, AllDerivatives, 2, L2H1Errors,
                                   BilinearCoeffs, aux, 1, fesp, errors);

#ifdef _MPI
      MPI_Allreduce(errors, reduced_errors, 2, MPI_DOUBLE, MPI_SUM, Comm);      
      for(i=0;i<2;i++)
	errors[i] = sqrt(reduced_errors[i]);
      if(rank == out_rank){
	OutPut(endl);
	OutPut("time: " << TDatabase::TimeDB->CURRENTTIME<<endl);
	OutPut( "L2 (u): " << sqrt(reduced_errors[0]) << endl);
	OutPut( "H1-semi (u): " << sqrt(reduced_errors[1]) << endl);
      }
#else
      OutPut(endl);
      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME<<endl);
      OutPut( "L2 (u): " << errors[0] << endl);
      OutPut( "H1-semi (u): " << errors[1] << endl);
#endif 
       
      if(m>1)
      {      
        errors[3] += (errors[0]*errors[0] + olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
	errors[4] += (errors[1]*errors[1] +olderror1 * olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
	
#ifdef _MPI
	if(rank == out_rank)
#endif
	{
	  OutPut(TDatabase::TimeDB->CURRENTTIME <<  " u-L2(0,T;L2) " << sqrt(errors[3]) << " ");
	  OutPut( "u-L2(0,T;H1) " << sqrt(errors[4]) << endl);
        }
        
        if(Linfty<errors[0])
          Linfty=errors[0];

#ifdef _MPI
        if(rank == out_rank)
#endif
	OutPut( "u-Linfty " << Linfty << endl);    
     } //  if(TDatabase::ParamDB->MEASURE_ERRORS) 
  } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

//======================================================================
// measure errors to known solution w
//======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {      
      w_olderror  = w_errors[0];
      w_olderror1 = w_errors[1];
 
      Scalar_W_FeFunction->GetErrors(W_Exact, 4, AllDerivatives, 2, L2H1Errors,
                                   W_BilinearCoeffs, aux, 1, fesp, w_errors);

#ifdef _MPI
      MPI_Allreduce(w_errors, reduced_errors, 2, MPI_DOUBLE, MPI_SUM, Comm);      
      for(i=0;i<2;i++)
	w_errors[i] = sqrt(reduced_errors[i]);
      if(rank == out_rank){
	OutPut(endl);
	OutPut("time: " << TDatabase::TimeDB->CURRENTTIME<<endl);
	OutPut( "L2 (w): " << sqrt(reduced_errors[0]) << endl);
	OutPut( "H1-semi (w): " << sqrt(reduced_errors[1]) << endl);
      }
#else
      OutPut(endl);
      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME<<endl);
      OutPut( "L2 (w): " << w_errors[0] << endl);
      OutPut( "H1-semi (w): " << w_errors[1] << endl);
#endif 
       
      if(m>1)
      {      
        w_errors[3] += (w_errors[0]*w_errors[0] + w_olderror * w_olderror )*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
	w_errors[4] += (w_errors[1]*w_errors[1] +w_olderror1 * w_olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
	
#ifdef _MPI
	if(rank == out_rank)
#endif
	{
	  OutPut(TDatabase::TimeDB->CURRENTTIME <<  " w-L2(0,T;L2) " << sqrt(w_errors[3]) << " ");
	  OutPut( "w-L2(0,T;H1) " << sqrt(w_errors[4]) << endl);
        }
        
        if(w_Linfty<w_errors[0])
          w_Linfty=w_errors[0];

#ifdef _MPI
        if(rank == out_rank)
#endif
	OutPut( "w-Linfty " << w_Linfty << endl);  
     } //  if(TDatabase::ParamDB->MEASURE_ERRORS) 
  } // while(TDatabase::TimeDB->CURRENTTIME< end_time)
  
//======================================================================
// measure errors to known solution w
//======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {      
      v_olderror  = v_errors[0];
      v_olderror1 = v_errors[1];
 

      Scalar_V_FeFunction->GetErrors(V_Exact, 4, AllDerivatives, 2, L2H1Errors,
                                   V_BilinearCoeffs, aux, 1, fesp, v_errors);

#ifdef _MPI
      MPI_Allreduce(v_errors, reduced_errors, 2, MPI_DOUBLE, MPI_SUM, Comm);      
      for(i=0;i<2;i++)
	v_errors[i] = sqrt(reduced_errors[i]);
      if(rank == out_rank){
	OutPut(endl);
	OutPut("time: " << TDatabase::TimeDB->CURRENTTIME<<endl);
	OutPut( "L2 (v): " << sqrt(reduced_errors[0]) << endl);
	OutPut( "H1-semi (v): " << sqrt(reduced_errors[1]) << endl);
      }
#else
      OutPut(endl);
      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME<<endl);
      OutPut( "L2 (v): " << v_errors[0] << endl);
      OutPut( "H1-semi (v): " << v_errors[1] << endl);
#endif 
       
      if(m>1)
      {      
        v_errors[3] += (v_errors[0]*v_errors[0] + v_olderror * v_olderror )*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
	v_errors[4] += (v_errors[1]*v_errors[1] +v_olderror1 * v_olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
	
#ifdef _MPI
	if(rank == out_rank)
#endif
	{
	  OutPut(TDatabase::TimeDB->CURRENTTIME <<  " v-L2(0,T;L2) " << sqrt(v_errors[3]) << " ");
	  OutPut( "v-L2(0,T;H1) " << sqrt(v_errors[4]) << endl);
        }
        
        if(v_Linfty<v_errors[0])
          v_Linfty=v_errors[0];

#ifdef _MPI
        if(rank == out_rank)
#endif
	OutPut( "v-Linfty " << v_Linfty << endl);  
     } //  if(TDatabase::ParamDB->MEASURE_ERRORS) 
  } // while(TDatabase::TimeDB->CURRENTTIME< end_time)  
  
//======================================================================
// produce final output
//======================================================================
  #ifdef _MPI
    if(profiling){
      total_vtk += (end_vtk-start_vtk);
      start_vtk = MPI_Wtime();
    }
  if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
    if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
        img++;

    if(profiling){
      end_vtk = MPI_Wtime();
      total_vtk += (end_vtk-start_vtk);
    }
#else
    if(profiling){
      total_vtk += (end_vtk-start_vtk);
      start_vtk = GetTime();
    }
   
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
    if(profiling){
      end_vtk = GetTime();
      total_vtk += (end_vtk-start_vtk);
    }
#endif
    }
  if(profiling){
#ifdef _MPI
    stop_time = MPI_Wtime();
#else
    stop_time = GetTime();
#endif
  }
//======================================================================
// Time profiling Output
//======================================================================  
  if(profiling){
#ifdef _MPI
    int Total_cells, Total_dof;
    MPI_Reduce(&N_Cells, &Total_cells, 1, MPI_INT, MPI_SUM, out_rank, Comm);
    MPI_Reduce(&N_DOF, &Total_dof, 1, MPI_INT, MPI_SUM, out_rank, Comm);
    N_Cells = Total_cells;
    N_DOF   = Total_dof;
    if(rank == out_rank){
#endif
    OutPut(endl<<"#Levels :: "<<LEVELS<<"  #Uniform refinement :: "<<TDatabase::ParamDB->UNIFORM_STEPS <<"   Order :: "<<TDatabase::ParamDB->ANSATZ_ORDER<<endl);  
    OutPut("Total Cells :: "<<N_Cells<<"     Total_dof :: "<<N_DOF<<endl<<endl);  
    OutPut("----------------------------------------------------------------------------------------------------------------------"<<endl); 
    OutPut( "Total time taken for initializing System Matrix : " << (total_int) << "("<<100*(total_int)/(stop_time-start_time)<<"%)"<<endl);
    OutPut( "Total time taken for vtk writing : " << (total_vtk) << "("<<100*(total_vtk)/(stop_time-start_time)<<"%)"<<endl);
    OutPut( "Total time taken for assembling : " << (total_assembling) << "("<<100*(total_assembling)/(stop_time-start_time)<<"%)"<<endl);
    OutPut( "Total time taken for solving : " << (total_solve) << "("<<100*(total_solve)/(stop_time-start_time)<<"%)"<<endl);
    OutPut( "Total time taken for communication : " << timeC << "(" <<100*timeC/(stop_time-start_time) <<"%)"<< endl);
    OutPut( "Total time taken throughout : " << (stop_time-start_time) << endl);
    OutPut("----------------------------------------------------------------------------------------------------------------------"<<endl);
#ifdef _MPI
    }
#endif
  }
      
  CloseFiles();
#ifdef _MPI
  MPI_Finalize(); 
#endif
  return 0;
} // end main





