/** ==========================================================================
#    This file is part of the finite element software ParMooN.
# 
#    ParMooN (cmg.cds.iisc.ac.in/parmoon) is a free finite element software  
#    developed by the research groups of Prof. Sashikumaar Ganesan (IISc, Bangalore),
#    Prof. Volker John (WIAS Berlin) and Prof. Gunar Matthies (TU-Dresden):
#
#    ParMooN is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    ParMooN is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with ParMooN. If not, see <http://www.gnu.org/licenses/>.
#
#    If your company is selling a software using ParMooN, please consider 
#    the option to obtain a commercial license for a fee. Please send 
#    corresponding requests to sashi@iisc.ac.in
# =========================================================================*/ 
   
/** ************************************************************************ 
* @brief     source file for TSystemNSE3D
* @author    Sashikumaar Ganesan, 
* @date      27.01.15
* @History 
 ************************************************************************  */
#include <Database.h>
#include <SystemHyperElast3D.h>
#include <SquareStructure3D.h>
#include <DiscreteForm3D.h>
#include <Assemble3D.h>
#include <FEVectFunct3D.h>
#include <AuxParam3D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
// #include <NSE3D_ParamRout.h>
#include <MainUtilities.h>
#include <Upwind.h>
#include <MultiGridIte.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <Upwind3D.h>
#include <AssembleMat3D.h>
#include <LinAlg.h>
// THIVIN #include <HyperElastic_MGLevel1.h>
//#include <HyperElastic_MGLevel.h>
// #include <NSE_MultiGrid.h>
// #include <NSE_MGLevel1.h>
// #include <NSE_MGLevel2.h>
// #include <NSE_MGLevel3.h>
// #include <NSE_MGLevel4.h>

#ifdef _MPI
#include <ParFEMapper3D.h>
#include <ParFECommunicator3D.h>
#endif

#define AMG 0
#define GMG 1
#define DIRECT 2

#include <stdlib.h>
#include <string.h>
// #include <sstream>
// #include <MooNMD_Io.h>

void swap(int a, int  b){
    int c;
    c = a; a = b; b = c;
}


void Galerkin3D(double Mult, double *coeff, double *param, double hK, double **OrigValues, int *N_BaseFuncts, double ***LocMatrices, double **LocRhs);
void HyperParamsVelo(double *in, double *out);
double Piola_Kir(double *param, double test100, double test010, double test001, double test000, double ansatz100, double ansatz010, double ansatz001, double ansatz000, int row, int col);

/**  parameters for Auxparam */
 int Hyper_FEFctIndex[9] = { 0, 1, 2, 0, 1, 2, 0, 1, 2};
 int Hyper_BeginParam[1] = { 0 };  
 MultiIndex3D Hyper_FEMultiIndex[9] = { D100, D100, D100, D010, D010, D010, D001, D001, D001 };  
 ParamFct *Hyper_ParamFct[1] = { HyperParamsVelo };


/** parameters for Discrete form */
  MultiIndex3D Derivatives[4] = { D100, D010, D001, D000};
  int SpaceNumbers[4] = { 0, 0, 0, 0};
  int RowSpace[9]    = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int ColumnSpace[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int RhsSpace[3] = { 0, 0, 0 };

TSystemHyperElast3D::TSystemHyperElast3D(int N_levels, TFESpace3D **disp_fespace, TFEVectFunct3D **displacement, double **sol, double **rhs, int disctype, int solver)
{
  int i, zerostart;
  int profiling = TDatabase::ParamDB->timeprofiling;
 
  /** need it for solver */
  sqmatrices = (TSquareMatrix **)SQMATRICES;
      
  //set number of multigrid levels
  N_Levels = N_levels;
   
//   cout << " N_levels " << N_levels <<endl;
//   exit(0);
//   
  //set the discretization type
  Disctype = disctype;
  
#ifdef _MPI
  Comm = TDatabase::ParamDB->Comm;
#endif  
  
  //set the solver type
  SOLVER = solver;
  
  Displacement = displacement;
 
  SolArray = sol;
  RhsArray = rhs;
  
  U_Space = disp_fespace;
  
  N_U = disp_fespace[N_levels-1]->GetN_DegreesOfFreedom();
  N_TotalDOF = 3*N_U;
       
  N_Active =  disp_fespace[N_levels-1]->GetActiveBound();
  N_DirichletDof = N_U - N_Active;  
 
  sqstructureA = new TSquareStructure3D *[N_levels];
   
  SqmatrixA11 = new TSquareMatrix3D*[N_levels];
  SqmatrixA12 = new TSquareMatrix3D*[N_levels]; 
  SqmatrixA13 = new TSquareMatrix3D*[N_levels];  
  SqmatrixA21 = new TSquareMatrix3D*[N_levels];
  SqmatrixA22 = new TSquareMatrix3D*[N_levels]; 
  SqmatrixA23 = new TSquareMatrix3D*[N_levels];    
  SqmatrixA31 = new TSquareMatrix3D*[N_levels];
  SqmatrixA32 = new TSquareMatrix3D*[N_levels]; 
  SqmatrixA33 = new TSquareMatrix3D*[N_levels];    
  
  AMatRhsAssemble = new TAssembleMat3D*[N_levels];  
  AMatAssembleNonLinear = new TAssembleMat3D*[N_levels];  
  
  if(TDatabase::ParamDB->INTERFACE_FLOW)
   {  
    SqmatrixF11 = new TSquareMatrix3D*[N_levels]; 
    SqmatrixF22 = new TSquareMatrix3D*[N_levels]; 
    SqmatrixF33 = new TSquareMatrix3D*[N_levels]; 
    
    N_FreeSurfFaces = new int[N_levels];
    FreeSurfCellNumbers = new int*[N_levels];
    FreeSurfJointNumbers = new int*[N_levels];
   }
  
  
  if(SOLVER==AMG_SOLVE || SOLVER==DIRECT)
   {
    Start_Level=N_Levels-1;
   }
  else 
   {
    Start_Level=0;
   }
   
  // build matrices   
  for(i=Start_Level;i<N_Levels;i++)
   {
    if(SOLVER==GMG)
     OutPut("MULTIGRID LEVEL : " << i<<endl;)
    
     // first build matrix structure
     sqstructureA[i] = new TSquareStructure3D(U_Space[i]);
     sqstructureA[i]->Sort();  // sort column numbers: numbers are in increasing order

     SqmatrixA11[i] = new TSquareMatrix3D(sqstructureA[i]);
     SqmatrixA12[i] = new TSquareMatrix3D(sqstructureA[i]);
     SqmatrixA13[i] = new TSquareMatrix3D(sqstructureA[i]);
     SqmatrixA21[i] = new TSquareMatrix3D(sqstructureA[i]);
     SqmatrixA22[i] = new TSquareMatrix3D(sqstructureA[i]);
     SqmatrixA23[i] = new TSquareMatrix3D(sqstructureA[i]);
     SqmatrixA31[i] = new TSquareMatrix3D(sqstructureA[i]);
     SqmatrixA32[i] = new TSquareMatrix3D(sqstructureA[i]);
     SqmatrixA33[i] = new TSquareMatrix3D(sqstructureA[i]);

//         MatVect = MatVect_NSE4;
          Defect = Defect_HE;

        if(TDatabase::ParamDB->INTERFACE_FLOW)
         {
          SqmatrixF11[i] = new TSquareMatrix3D(sqstructureA[i]);
          SqmatrixF22[i] = new TSquareMatrix3D(sqstructureA[i]);
          SqmatrixF33[i] = new TSquareMatrix3D(sqstructureA[i]);
         }
    } //  for(i=Start_Level;i<N_Levels;i++)
    Defect = Defect_HE;
#ifdef _MPI
  i = N_Levels-1; 
 

  double t1,t2,tdiff;
  if(profiling)  t1 = MPI_Wtime();

  ParMapper_U = new TParFEMapper3D*[N_levels]; 
  ParComm_U = new TParFECommunicator3D*[N_levels];
  
  for(i=Start_Level;i<N_levels;i++)
  {
    ParMapper_U[i] = new TParFEMapper3D(3, U_Space[i], sqstructureA[i]->GetRowPtr(),sqstructureA[i]->GetKCol());

    ParComm_U[i] = new TParFECommunicator3D(ParMapper_U[i]);
  }


  if(profiling)
  {
     t2 = MPI_Wtime();
     tdiff = t2-t1;
     int out_rank=TDatabase::ParamDB->Par_P0;
     int rank;
     MPI_Comm_rank(Comm, &rank);
     MPI_Reduce(&tdiff, &t1, 1, MPI_DOUBLE, MPI_MIN, out_rank, Comm);
     if(rank == out_rank)
     {
      printf("Time taken for ParFEMapper3D for all levels is %e\n", t1);
     }
  }
#endif

   //initialize multigrid solver
//    if(SOLVER==GMG)
//    {    
//     Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
//     Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
//     MG = new TNSE_MultiGrid(1, 2, Parameters);
//     
//      switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
//       {
//        case 11:
//           zerostart = 0;
//        break;
//        case 16:
//            zerostart = 1;
//         break;
//       default:
//          zerostart = 1;
//        }    
//     
//     // build preconditioner
//     switch (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE)
//      {
//       case 5:
//        prec = new TMultiGridIte(MatVect, Defect, NULL, 0, N_TotalDOF, MG, zerostart);
//        Itmethod_sol = new double[N_TotalDOF];
//        Itmethod_rhs = new double[N_TotalDOF];
//       break;
//       default:
//          OutPut("Unknown preconditioner !!!" << endl);
//          exit(4711);
//      }
//      
//     // build solver
//     switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
//      {
//         case 11:
//             Itmethod = new TFixedPointIte(MatVect, Defect, prec, 0, N_TotalDOF, 0);
//          break;
//          case 16:
//             Itmethod = new TFgmresIte(MatVect, Defect, prec, 0, N_TotalDOF, 0
//             #ifdef _MPI 
// ,ParComm_U[N_levels-1], ParComm_P[N_levels-1]
// #endif
// );
//          break;
//          default:
//             OutPut("Unknown solver !!!" << endl);
//             exit(4711);
//       }         
//      
//         
//      } //  if(SOLVER==GMG)  
       
  
    fefct_aux = new  TFEFunction3D*[N_levels*3]; 
    Hyperaux = new TAuxParam3D*[N_levels];
    Hyperaux_error = new TAuxParam3D*[N_levels];    

   for(i=Start_Level;i<N_Levels;i++)
   {
    Hyperaux_error[i] = NULL;
    Hyperaux[i] =  NULL;         
   }
   
} // TSystemHyperElast3D

// TSystemHyperElast3D::~TSystemHyperElast3D()
// {
//     delete [] Hyperaux; 
//        
//     if(Hyperaux_error)
//       delete [] Hyperaux_error;
// 
// }


void TSystemHyperElast3D::Init(CoeffFct3D *lincoeffs, BoundCondFunct3D *BoundCond, BoundValueFunct3D *U1BoundValue, 
                           BoundValueFunct3D *U2BoundValue, BoundValueFunct3D *U3BoundValue)
{ 
  int i, N_SquareMatrices, N_Rhs, N_FESpaces;
  int N_U_Current;
  int N_Terms = 4;
  int N_Matrices = 9;
   
  char GalerkinString[] = "Galerkin";
  char rhs[] = "rhs";
  char all[] = "all";
  
//   int mg_type = TDatabase::ParamDB->SC_MG_TYPE_SADDLE;
//     
//   double alpha[2];  
  
  TDiscreteForm3D *DiscreteFormGalerkin;
 
    
  /**  save the boundary condition */
  BoundaryConditions[0] = BoundCond;
  BoundaryConditions[1] = BoundCond;  
  BoundaryConditions[2] = BoundCond;  
  
  /**  save the boundary values   */
  BoundaryValues[0] = U1BoundValue;
  BoundaryValues[1] = U2BoundValue;
  BoundaryValues[2] = U3BoundValue;
  
  /** save the nse bilinear coefficient   */
  LinCoeffs[0] = lincoeffs;
  
  N_Rhs = 3;
  DiscreteFormGalerkin = new TDiscreteForm3D(GalerkinString, all,  N_Terms, Derivatives, SpaceNumbers,
                                             N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace, Galerkin3D, LinCoeffs[0], NULL); 
  
    /** find discrete form */
    switch(Disctype)
       {
          case GALERKIN:
            DiscreteFormARhs = DiscreteFormGalerkin;
          break;
 
          default:
            Error("Unknown DISCTYPE" << endl);
            exit(-1);
        } 
 
 
#ifdef _OMPONLY
   if(SOLVER == DIRECT && TDatabase::ParamDB->DSType == 1)
   {
     cout<<"NOT YET IMPLEMENTED !!!"<<endl;
     exit(0);
   }
#endif  
              
   // initilize the assemble    
   for(i=Start_Level;i<N_Levels;i++)
    {     
     SQMATRICES[0] = SqmatrixA11[i];
     SQMATRICES[1] = SqmatrixA12[i];
     SQMATRICES[2] = SqmatrixA13[i];	  
     SQMATRICES[3] = SqmatrixA21[i];
     SQMATRICES[4] = SqmatrixA22[i];
     SQMATRICES[5] = SqmatrixA23[i]; 
     SQMATRICES[6] = SqmatrixA31[i];
     SQMATRICES[7] = SqmatrixA32[i];
     SQMATRICES[8] = SqmatrixA33[i];  
 
     N_SquareMatrices = 9;

      N_Rhs = 3;
      N_FESpaces = 1;   
      N_U_Current = U_Space[i]->GetN_DegreesOfFreedom();   
      
      this->InitHyperAuxParm(i);
           
      fesp[0] =  U_Space[i];
 
      fefct[0] = Displacement[i]->GetComponent(0);
      fefct[1] = Displacement[i]->GetComponent(1);
      fefct[2] = Displacement[i]->GetComponent(2);      
     
      fesprhs[0] =  U_Space[i];
      fesprhs[1] =  U_Space[i];
      fesprhs[2] =  U_Space[i];

      RHSs[0] = RhsArray[i];
      RHSs[1] = RhsArray[i] + N_U_Current;
      RHSs[2] = RhsArray[i] + 2*N_U_Current;
//    cout << " Hyper_FEMultiIndex   3D: " << Hyper_FEMultiIndex[0] << endl;
  
     // array of assemble objects
     AMatRhsAssemble[i] = new TAssembleMat3D(N_FESpaces, fesp, N_SquareMatrices, SQMATRICES, 0, NULL,
                                             N_Rhs, RHSs, fesprhs, DiscreteFormARhs, BoundaryConditions, BoundaryValues, Hyperaux[i]);
     AMatRhsAssemble[i]->Init();    
        
//      if(TDatabase::ParamDB->INTERFACE_FLOW)
//       {
//        this->FindFreeSurfJoints(i, 0);
//       }
        
#ifdef _MPI
   if(i == N_Levels-1) {
    if(SOLVER == DIRECT)
     {
      DS = new TParDirectSolver(ParComm_U[N_Levels-1], ParComm_P[N_Levels-1], SQMATRICES, NULL);
     }
   }
#endif
    
     } // for(i=Start_Level;i<N_Levels;i++)      
//               
  cout << " TSystemHyperElast3D::Init done ! " << endl; 
              
} // TSystemHyperElast3D::Init


void TSystemHyperElast3D::Assemble()
{
  int i, N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
  int N_U_Current,  N_Active_Current, N_DirichletDof;
    
  double alpha[2];

//   
   for(i=Start_Level;i<N_Levels;i++)
    {     
      N_U_Current = U_Space[i]->GetN_DegreesOfFreedom();
      N_Active_Current  = U_Space[i]->GetActiveBound();     
      N_DirichletDof = N_U_Current - N_Active_Current;
      //cout<<"N_U_Current Active Dirichlet: "<<N_U_Current<<" "<<N_Active_Current<<" "<<N_DirichletDof<<endl;
      // initialize matrices
      AMatRhsAssemble[i]->Reset();

      /** assemble */
      AMatRhsAssemble[i]->Assemble3D();
      
       /** free surface/interface integration */
 //       if(TDatabase::ParamDB->INTERFACE_FLOW)
 //        {      
 //         FreeSurfInt(U_Space[i]->GetCollection(), N_FreeSurfFaces[i], FreeSurfCellNumbers[i], FreeSurfJointNumbers[i],
     
 	 
 /*		      if(N_FreeJoints)
    {
     FreeSurfCellNumbers[level] = new int[N_FreeJoints];
     FreeSurfJointNumbers[level] = new int[N_FreeJoints];*/  
 
 //   void FreeSurfInt(TCollection *Coll, int N_BoundFaces,
 //                  int *CellNumbers, int *JointNumbers,
 //                  TFEFunction3D *potential, double dt,
 //                  TSquareMatrix3D *Aii,
 //                  double *rhs1, double *rhs2, double *rhs3)
   
 //        }
 
           
 //          AMatRhsAssemble[i]->AssembleNavierSlip(); 
           
           // prepare everything for the assembling of slip with friction bc
           // on all levels
/*           N_FESpaces = 1;
           N_SquareMatrices = 9;
           N_RectMatrices = 0;
           N_Rhs = 3; 
 
           fesp[0] =  U_Space[i];
 	  
           fesprhs[0] =  U_Space[i];
           fesprhs[1] =  U_Space[i];
           fesprhs[2] =  U_Space[i];
       
           SQMATRICES[0] = SqmatrixA11[i];
           SQMATRICES[1] = SqmatrixA22[i];
           SQMATRICES[2] = SqmatrixA33[i];
           SQMATRICES[3] = SqmatrixA12[i];
           SQMATRICES[4] = SqmatrixA13[i];
           SQMATRICES[5] = SqmatrixA21[i];
           SQMATRICES[6] = SqmatrixA23[i];
           SQMATRICES[7] = SqmatrixA31[i];
           SQMATRICES[8] = SqmatrixA32[i];
 
           RHSs[0] = RhsArray[i];
           RHSs[1] = RhsArray[i] + N_U_Current;
           RHSs[2] = RhsArray[i] + 2*N_U_Current;
	   cout<<"Size of RHS "<<sizeof(SQMATRICES[0])<<endl;
       
//           Assemble3DSlipBC(N_FESpaces, fesp,
//                            N_SquareMatrices, SQMATRICES,
//                            N_RectMatrices, NULL,
//                            N_Rhs, RHSs, fesprhs,
//                            NULL,
//                            BoundaryConditions,
//                            BoundaryValues,
//                            NSEaux[i]);
//         } //  if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRIC  */
//     
//      set rhs for Dirichlet nodes
/*       memcpy(SolArray[i]+N_Active_Current, RhsArray[i]+N_Active_Current, N_DirichletDof*SizeOfDouble);
       memcpy(SolArray[i]+N_U_Current+N_Active_Current, RhsArray[i]+N_U_Current+N_Active_Current, N_DirichletDof*SizeOfDouble); 
       memcpy(SolArray[i]+2*N_U_Current+N_Active_Current, RhsArray[i]+2*N_U_Current+N_Active_Current, N_DirichletDof*SizeOfDouble);   */  
        
     } // for(i=Start_Level;i<N_Levels;i++)
//     cout << "Test Assemble " << endl; 
//     exit(0);
    
} // TSystemHyperElast3D::Assemble(T

void TSystemHyperElast3D::InitHyperAuxParm(int i)
{
  int Hyper_N_FESpace = 1;
  int Hyper_N_FEFunction = 3;
  int Hyper_N_ParamFct = 1;
  int Hyper_N_FEValues = 9;
  int Hyper_N_ParamValues = 54;                                      
                                            
  fesp_aux[0] =  U_Space[i];
     
  fefct_aux[i*3 ] = Displacement[i]->GetComponent(0);
  fefct_aux[(i*3)+1] = Displacement[i]->GetComponent(1);
  fefct_aux[(i*3)+2] = Displacement[i]->GetComponent(2);

  Hyperaux[i] =  new TAuxParam3D(Hyper_N_FESpace, Hyper_N_FEFunction, Hyper_N_ParamFct, Hyper_N_FEValues,
                                fesp_aux, fefct_aux+(i*3), Hyper_ParamFct, Hyper_FEFctIndex, Hyper_FEMultiIndex,
                                Hyper_N_ParamValues, Hyper_BeginParam);
}



void Galerkin3D(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA13, **MatrixA21, **MatrixA22, **MatrixA23; 
  double **MatrixA31, **MatrixA32, **MatrixA33;
  double *Rhs1, *Rhs2, *Rhs3, val, val1;
  double *MatrixRowA11, *MatrixRowA12, *MatrixRowA13, *MatrixRowA21, *MatrixRowA22, *MatrixRowA23;
  double *MatrixRowA31, *MatrixRowA32, *MatrixRowA33;
  double ansatz100, ansatz010, ansatz001, ansatz000;
  double test000, test100, test010, test001;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U;
  double nu, f1, f2, f3;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA13 = LocMatrices[2];
  MatrixA21 = LocMatrices[3];
  MatrixA22 = LocMatrices[4];
  MatrixA23 = LocMatrices[5];    
  MatrixA31 = LocMatrices[6];  
  MatrixA32 = LocMatrices[7];  
  MatrixA33 = LocMatrices[8];


  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // phi_x
  Orig1 = OrigValues[1]; // phi_y
  Orig2 = OrigValues[2]; // phi_z
  Orig3 = OrigValues[3]; // phi


  nu = coeff[0]; // nu
  f1 = coeff[1]; // f1
  f2 = coeff[2]; // f2
  f3 = coeff[3]; // f3
 
/*  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  u3 = param[2]; // u3old
  u1x = param[3]; // u1old
  u2x = param[4]; // u2old
  u3x = param[5]; // u3old 
  u1y = param[6]; // u1old
  u2y = param[7]; // u2old
  u3y = param[8]; // u3old  
  u1z = param[9]; // u1old
  u2z = param[10]; // u2old
  u3z = param[11]; // u3old */ 
 
  double Sder[6];
  double g[9];
  double h[9];
  double F[9];
  double v1[3];
  double H_symm[6];
  int m, k, l, d;
  double S_der[9];
  
  for(i=0;i<N_U;i++)
  {
      
    MatrixRowA11 = MatrixA11[i];
    MatrixRowA12 = MatrixA12[i];
    MatrixRowA13 = MatrixA13[i];
    MatrixRowA21 = MatrixA21[i];    
    MatrixRowA22 = MatrixA22[i];    
    MatrixRowA23 = MatrixA23[i];    
    MatrixRowA31 = MatrixA31[i];    
    MatrixRowA32 = MatrixA32[i];
    MatrixRowA33 = MatrixA33[i];    
    
    test100 = Orig0[i];
    test010 = Orig1[i];
    test001 = Orig2[i];
    test000 = Orig3[i];
    
    val1 = Mult*test000;
    
    Rhs1[i] += val1*f1;
    Rhs2[i] += val1*f2;
    Rhs3[i] += val1*f3;
    //cout<<"RHS i= "<<i<<" "<<Rhs1[i]<<" "<<Rhs2[i]<<" "<<Rhs3[i]<<endl;
    
    for(j=0;j<N_U;j++)
    { 
      
      ansatz100 = Orig0[j];
      ansatz010 = Orig1[j];
      ansatz001 = Orig2[j];
      ansatz000 = Orig3[j];

      for(m =0; m<9; m++)
         F[m] = param[45+m]; 
      //===============Blocks A11, A21, and A31 ================//
      //Construction of derivative of (F^t)F matrix ; H;
      v1[0] = param[45]; v1[1] = param[48]; v1[2] = param[51];
      
      for(k =0; k < 3; k++){
       d = 3*k;
       h[d] = ansatz100*v1[k];
       h[d+1] = ansatz010*v1[k];
       h[d+2] = ansatz001*v1[k];
       }
      for(k =0; k < 3; k++){
       h[k] += ansatz100*v1[k];
       h[k+3] += ansatz010*v1[k];
       h[k+6] += ansatz001*v1[k];
       }
       d = 0;
       for(k =0; k< 3; k++)
       {
            for(l =0; l< 3; l++){
               if( k>=l){
                   H_symm[d] = h[3*l +k];
                   d++;
               }
           }
       } 
       
       
      //=====================================================//
      d = 9;
      for( k =0; k < 6; k++){
        Sder[k] = 0;
        for(l =0 ; l< 6; l++){ 
        Sder[k] += (H_symm[l]*param[d]);
        d++;
         }
       }
       
//        cout<<param[15]<<" "<<param[16]<<" "<<param[17]<<" "<<param[18]<<" "<<param[19]<<" "<<param[10]<<endl;
//        cout<<H_symm[0]<<" "<<H_symm[1]<<" "<<H_symm[2]<<" "<<H_symm[3]<<" "<<H_symm[4]<<" "<<H_symm[5]<<endl;
//        cout<<Sder[0]<<" "<<Sder[1]<<" "<<Sder[2]<<" "<<Sder[3]<<" "<<Sder[4]<<" "<<Sder[5]<<endl;
       
       S_der[0] = Sder[0]; S_der[1] = Sder[1]; S_der[2] = Sder[2];
       S_der[3] = Sder[1]; S_der[4] = Sder[3]; S_der[5] = Sder[4];
       S_der[6] = Sder[2]; S_der[7] = Sder[4]; S_der[8] = Sder[5];
      memset(g, 0, 9*sizeof(double));
      //MatrixMult(F, S_der,g, 'n', 'n') ;
      for(int x = 0 ; x < 3; x++){
          for(int y = 0; y < 3; y++){
              for(int w = 0; w< 3; w++){
                  g[3*y + x] += F[3*w + x] * S_der[3*y + w];
              }
          }
      }
      //cout<<"Matrix F and S_der"<<endl;
      //for(int as = 0; as < 9; as++) cout<<F[as]<<" "<<S_der[as]<<" "<<g[as]<<"\t";
      g[0] += (ansatz100*param[0] + ansatz010*param[1] + ansatz001*param[2]);
      g[3] += (ansatz100*param[3] + ansatz010*param[4] + ansatz001*param[5]);
      g[6] += (ansatz100*param[6] + ansatz010*param[7] + ansatz001*param[8]);
      
      MatrixRowA11[j] += Mult * (test100*g[0] + test010*g[3] + test001*g[6]);
      MatrixRowA21[j] += Mult * (test100*g[1] + test010*g[4] + test001*g[7]);
      MatrixRowA31[j] += Mult * (test100*g[2] + test010*g[5] + test001*g[8]);
      //cout<<MatrixRowA11[j]<<" "<<MatrixRowA21[j]<<" "<<MatrixRowA31[j]<<endl;
//       if(i == 2){
//           cout<<test100<<" "<<test010<<" "<<test001<<endl;
//           for(int as = 0; as < 9; as++) cout<<F[as]<<" "<<S_der[as]<<" "<<g[as] <<"\t";
//       }
      
//       if( i== N_U-1){
//           cout<<"Full matrix at the end "<<N_U<<endl;
//           for(int x = 0 ; x < N_U; x++) cout<<MatrixRowA11[x]<<endl;
//       }
      //===============Blocks A12, A22, and A32 ================//      
      //Construction of derivative of (F^t)F matrix ; H;
      v1[0] = param[46]; v1[1] = param[49]; v1[2] = param[52];
      for(k =0; k < 3; k++){
       d = 3*k;
       h[d] = ansatz100*v1[k];
       h[d+1] = ansatz010*v1[k];
       h[d+2] = ansatz001*v1[k];
       }

      for( k =0; k < 3; k++){
       h[k] += ansatz100*v1[k];
       h[k+3] += ansatz010*v1[k];
       h[k+6] += ansatz001*v1[k];
       }
      
       d = 0;
       for(k =0; k< 3; k++){
           for(l =0; l< 3; l++){
               if( k>=l){
                   H_symm[d] = h[3*l +k];
                   d++;
               }
           }
       }

//        cout<<ansatz100<<" "<<ansatz010<<" "<<ansatz001<<endl;
//        for(k =0; k< 9; k++) cout<<h[k]<<"\t";    
//        cout<<endl;
//        for(k =0; k< 6; k++) cout<<H_symm[k]<<"\t";
//        cout<<endl;
      //=====================================================//
      d = 9;
      for( k =0; k < 6; k++){
        Sder[k] = 0;
        for(l =0 ; l< 6; l++){ 
        Sder[k] += H_symm[l]*param[d];
        d++;
         }
       }
      memset(g, 0, 9*sizeof(double));

       S_der[0] = Sder[0]; S_der[1] = Sder[1]; S_der[2] = Sder[2];
       S_der[3] = Sder[1]; S_der[4] = Sder[3]; S_der[5] = Sder[4];
       S_der[6] = Sder[2]; S_der[7] = Sder[4]; S_der[8] = Sder[5];
       MatrixMult(F, S_der,g, 'n', 'n') ;
       g[1] += (ansatz100*param[0] + ansatz010*param[1] + ansatz001*param[2]);
       g[4] += (ansatz100*param[3] + ansatz010*param[4] + ansatz001*param[5]);
       g[7] += (ansatz100*param[6] + ansatz010*param[7] + ansatz001*param[8]);
       MatrixRowA12[j] += Mult * (test100*g[0] + test010*g[3] + test001*g[6]);
       MatrixRowA22[j] += Mult * (test100*g[1] + test010*g[4] + test001*g[7]);
       MatrixRowA32[j] += Mult * (test100*g[2] + test010*g[5] + test001*g[8]);
       
       
//        cout<<param[15]<<" "<<param[16]<<" "<<param[17]<<" "<<param[18]<<" "<<param[19]<<" "<<param[10]<<endl;
//        cout<<H_symm[0]<<" "<<H_symm[1]<<" "<<H_symm[2]<<" "<<H_symm[3]<<" "<<H_symm[4]<<" "<<H_symm[5]<<endl;
//        cout<<Sder[0]<<" "<<Sder[1]<<" "<<Sder[2]<<" "<<Sder[3]<<" "<<Sder[4]<<" "<<Sder[5]<<endl;       

      //===============Blocks A13, A23, and A33 ================//

      v1[0] = param[47]; v1[1] = param[50]; v1[2] = param[53];
      for( k =0; k < 3; k++){
       d = 3*k;
       h[d] = ansatz100*v1[k];
       h[d+1] = ansatz010*v1[k];
       h[d+2] = ansatz001*v1[k];
       }
      for( k =0; k < 3; k++){
       h[k] += ansatz100*v1[k];
       h[k+3] += ansatz010*v1[k];
       h[k+6] += ansatz001*v1[k];
       }
       d = 0;
       for(k =0; k< 3; k++){
           for(l =0; l< 3; l++){
               if( k>=l){
                   H_symm[d] = h[3*l +k];
                   d++;
               }
           }
       }

      //=====================================================//
      d = 9;
      for( k =0; k < 6; k++){
        Sder[k] = 0;
        for(l =0 ; l< 6; l++){ 
        Sder[k] += H_symm[l]*param[d];
        d++;
         }
       }
        memset(g, 0, 9*sizeof(double));
                             
       S_der[0] = Sder[0]; S_der[1] = Sder[1]; S_der[2] = Sder[2];
       S_der[3] = Sder[1]; S_der[4] = Sder[3]; S_der[5] = Sder[4];
       S_der[6] = Sder[2]; S_der[7] = Sder[4]; S_der[8] = Sder[5];
       MatrixMult(F, S_der,g, 'n', 'n') ;
       g[2] += (ansatz100*param[0] + ansatz010*param[1] + ansatz001*param[2]);
       g[5] += (ansatz100*param[3] + ansatz010*param[4] + ansatz001*param[5]);
       g[8] += (ansatz100*param[6] + ansatz010*param[7] + ansatz001*param[8]);
       MatrixRowA13[j] += Mult * (test100*g[0] + test010*g[3] + test001*g[6]);
       MatrixRowA23[j] += Mult * (test100*g[1] + test010*g[4] + test001*g[7]);
       MatrixRowA33[j] += Mult * (test100*g[2] + test010*g[5] + test001*g[8]);
      //cout<<ansatz100<<"\t"<<ansatz010<<"\t"<<ansatz001<<"\t"<<test100<<"\t"<<test010<<"\t"<<test001<<endl;
      //cout<<"MR "<<i<<" "<<j<<" "<<MatrixRowA11[j]<<"\t";//<<MatrixRowA12[j]<<"\t"<<MatrixRowA13[j]<<"\n"
      //                           <<MatrixRowA21[j]<<"\t"<<MatrixRowA22[j]<<"\t"<<MatrixRowA23[j]<<"\n"
      //                           <<MatrixRowA31[j]<<"\t"<<MatrixRowA32[j]<<"\t"<<MatrixRowA33[j]<<endl;
    } // endfor j
    //cout<<endl;


  } // endfor i

//      for (i = 0; i < 54; i++)
//         cout << "param " << i << ": " << param[i] <<endl;
   
// cout << " Galerkin3D " <<endl;
// exit(0);
 
 
}


void HyperParamsVelo(double *in, double *out)
{  
  /** based on the elastic model, S_ij, 0< i,j < 9 vales will be returned */

  int i, j, k, l, J, K, L, d;
  const double c1= 100;
  const double c2= 6000 ;
  const double D1=pow(10, -3);
 
  double I1, I3, I3pow, k1, k2, I1_bar, *OUT;
  double  F[9], C[9];
//   double cij, cik, cjk, ckl, cil, cjl; 
  int flag_ik = 0, flag_jk = 0, flag_il = 0, flag_jl = 0 ;
  int del_ij, del_kl;
  
  //cout<<" Input parameter "<<in[2]<<endl;
  // grad U
  for(int i =0; i<9; i++){
   // in[3+i] = 0.001*i + 0.001;
    F[i] = in[3+i];
   }
  
  // grad U + I
  F[0] += 1.0; 
  F[4] += 1.0;
  F[8] += 1.0;
  
  // C = F^T F
  MatrixMult(F, F, C, 't', 'n');
 
  // I1 = tr(C)
  I1 = C[0] +C[4] + C[8];
  
  // I3=det(C), Blas will retun LU decomposed value in C
  I3 = MatrixDeterminant(C, 3);
  I3pow = pow(I3, -1.0/3.0);


  //cout<<"trace "<<I1<<" "<<I3<<" "<<I3pow<<endl;
  // compute C^{-1}
  MatrixInverse(C, 3, 1);
  
  I1_bar = I1 * I3pow;
  
  for(i=0; i<9; i++)
   out[i] = 0.;
       
  k1 = 2 * (c1 + 2*c2*I1_bar - 6.*c2);

  k2 = (4.0 * (I3-1.) * I3)/D1;
  //cout<<"k1 = "<<k1 <<" k2 ="<<k2<<endl;
  for(i=0; i<9; i++)
   {
     if(i == 0 || i == 4 || i == 8)  
      out[i] = k1*(I3pow - (I1_bar*C[i])/3.0) + (k2* C[i]);
     else 
      out[i] = k1*(- (I1_bar*C[i])/3.0) + (k2* C[i]);
    //cout<<"S["<<i<<"] = "<<out[i]<<endl;
   }
 OUT = out+9;
 d =0;
   k1 = (c1 + 2*c2*I1_bar - 6.*c2);

   for(i = 0; i< 3; i++)
  { 
       
     for( j = 0; j< 3; j++)
      {
          
          del_ij = 0;
          if(i >= j){
              
              if(i == j) del_ij =1;
              k2 = 2 * c2 * (I3pow * del_ij - ( I1_bar*C[3*j+i] )/3.0 );
              //cout<<k1<<endl;
              for ( k =0; k < 3; k++)
              {
              
                  if( k > i){
                      swap(i, k); 
                      flag_ik = 1;
                  }
                  if( k > j){
                      swap(j, k); 
                      flag_jk = 1;
                  }        
                 for( l =0; l < 3; l++)
                  {
                      //cout<<"Check";
                      del_kl = 0;
                      if(k == l) del_kl =1;
                      if(k >= l){
                          
                          if( l > j){
                              swap(j, l); 
                              flag_jl = 1;
                          }
                          
                          if( l > i){
                              swap(i, l); 
                              flag_il = 1;
                          }
                          double a = k2*(I3pow * del_kl - (I1_bar * C[3*l + k]) /3.0 );
                          double b = k1 * ( (I1_bar * del_ij * C[3*l + k]) +  ( (-I1_bar * C[3*l +k])/3.0 + I3pow * del_kl ) * C[3*j +i]  +
                                    I1_bar * ( C[3*k +i] * C[3*l +j] + C[3*k +j] * C[3*i +l] ) /2.0 ) / 3.0;
                          double c = (2./D1) * ( I3*(I3-1)*C[3*l + k] *C[3 *j +i] + I3*I3*C[3*l +k] *C[3 *j +i] + I3*(I3-1)*(C[3*k +i] * C[3*l +j] + C[3*k +j] * C[3*i +l])/2.0);
                          
                          OUT[d] =  4 * ( a - b + c  );
                          d++;
                          //cout<<d<<" "<<a<<" "<<b<<" "<<c<<" "<<out[d-1]<<endl;
                      }
                      if(flag_jl == 1){
                          swap(j, l); 
                          flag_jl = 0;
                      }
                      if(flag_il == 1){
                          swap(i, l); 
                          flag_il = 0;
                      }
                  }
                  if(flag_ik == 1){
                          swap(i, k); 
                          flag_ik = 0;
                  }
                  if(flag_jk == 1){
                          swap(j, k); 
                          flag_jk = 0;
                  }
              }
          }
      }
  }
 

  
 

    for(i =0; i<9; i++)
      out[45+i] = F[i];
    
//   out[9] = in[0]; // x - coordinate  
//   out[10] = in[1]; // y - coordinate  
//   out[11] = in[2]; // z - coordinate 
/*  for (i = 0; i < 54; i++)
    cout << "out " << i << ": " << out[i] <<endl;
 */  
//  cout << " stress cal " <<endl;
//  exit(0);
 
}

void TSystemHyperElast3D::GetResidual(double *sol, double *rhs, double *res, double &impuls_residual, double &residual)
{

     
     SQMATRICES[0] = SqmatrixA11[N_Levels-1];
     SQMATRICES[1] = SqmatrixA12[N_Levels-1];
     SQMATRICES[2] = SqmatrixA13[N_Levels-1];	  
     SQMATRICES[3] = SqmatrixA21[N_Levels-1];
     SQMATRICES[4] = SqmatrixA22[N_Levels-1];
     SQMATRICES[5] = SqmatrixA23[N_Levels-1]; 
     SQMATRICES[6] = SqmatrixA31[N_Levels-1];
     SQMATRICES[7] = SqmatrixA32[N_Levels-1];
     SQMATRICES[8] = SqmatrixA33[N_Levels-1];  

#ifdef _MPI      
    ParComm_U[N_Levels-1]->CommUpdate(sol);
    ParComm_P[N_Levels-1]->CommUpdate(sol+3*N_U);
#endif
    
    //THIVIN Defect(sqmatrices, matrices, sol, rhs, res); 

#ifdef _MPI
   double residual_scalar = 0.0;
   double sum =0.;
   int i,j,rank;
   MPI_Comm_rank(Comm, &rank);
   int *master = ParComm_U[N_Levels-1]->GetMaster();
   
   for(i=0;i<N_U;i++)
   {
     if(master[i]!=rank)    continue;
      
      residual_scalar += res[i      ]*res[i      ];
      residual_scalar += res[i+  N_U]*res[i+  N_U];
      residual_scalar += res[i+2*N_U]*res[i+2*N_U];

    }

    MPI_Allreduce(&residual_scalar, &sum, 1, MPI_DOUBLE, MPI_SUM, Comm);
    impuls_residual = (sum);

    master = ParComm_P[N_Levels-1]->GetMaster();
    for(i=0;i<N_P;i++)
    {
      if(master[i]!=rank)    continue;
      
      residual_scalar += res[i+3*N_U]*res[i+3*N_U];

    }
    
    sum = 0;
    MPI_Allreduce(&residual_scalar, &sum, 1, MPI_DOUBLE, MPI_SUM, Comm);
    residual = (sum);

#else
    impuls_residual  =  Ddot(3*N_U, res, res);
    residual         =  Ddot(N_TotalDOF, res, res); 

#endif
  
} // TSystemHyperElast3D::GetResidual



void TSystemHyperElast3D::Solve(double *sol, double *rhs)
 {
   int N_LinIter=0;
   double summ = 0;
   double residual, residual_scalar = 0.0;
   double sum =0.;
   int i,j,rank;
   int *master;

   
   double *EntriesA11, *EntriesA12, *EntriesA13, *EntriesA21;
  double *EntriesA22, *EntriesA23, *EntriesA31, *EntriesA32, *EntriesA33;
   EntriesA11 = SqmatrixA11[N_Levels-1]->GetEntries();
   EntriesA12 = SqmatrixA12[N_Levels-1]->GetEntries();
   EntriesA13 = SqmatrixA13[N_Levels-1]->GetEntries();
   EntriesA21 = SqmatrixA21[N_Levels-1]->GetEntries();
   EntriesA22 = SqmatrixA22[N_Levels-1]->GetEntries();
   EntriesA23 = SqmatrixA23[N_Levels-1]->GetEntries();
   EntriesA31 = SqmatrixA31[N_Levels-1]->GetEntries();
   EntriesA32 = SqmatrixA32[N_Levels-1]->GetEntries();
   EntriesA33 = SqmatrixA33[N_Levels-1]->GetEntries();

  
  std::ofstream file;
  file.open ("/home/stuti/ParMooN_All/ParMooN_Raw/Output/mat_before_sol.csv", std::ios::app);
  for(i=0;i<27;i++)
     file << EntriesA11[i] <<",";
  file<<"\n";

  for(j=0;j<27;j++)
     file << EntriesA12[j] <<",";
  file<<"\n";

  for(j=0;j<27;j++)
     file << EntriesA13[j] <<",";
  file<<"\n";

  for(j=0;j<27;j++)
     file << EntriesA21[j] <<",";
  file<<"\n";

  for(j=0;j<27;j++)
     file << EntriesA22[j] <<",";
  file<<"\n";

  for(j=0;j<27;j++)
     file << EntriesA23[j] <<",";
  file<<"\n";

  for(j=0;j<27;j++)
     file << EntriesA31[j] <<",";
  file<<"\n";
 
  for(j=0;j<27;j++)
     file << EntriesA32[j] <<",";
  file<<"\n";
           
  for(j=0;j<27;j++)
     file << EntriesA33[j] <<",";
  file<<"\n";
    
    
    file.close();
     switch(SOLVER)
      {
       /*case AMG_SOLVE:
         cout << "AMG_SOLVE not yet implemented " <<endl;
       break;
 
       case GMG:
 
           if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
            {
             memcpy(Itmethod_sol, sol, N_TotalDOF*SizeOfDouble);
             memcpy(Itmethod_rhs, rhs, N_TotalDOF*SizeOfDouble);
            }
           // solve the linear system
           N_LinIter += Itmethod->Iterate(sqmatrices, matrices, Itmethod_sol, Itmethod_rhs);
 	  
           if(TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
            {
             memcpy(sol, Itmethod_sol, N_TotalDOF*SizeOfDouble);
             memcpy(rhs, Itmethod_rhs, N_TotalDOF*SizeOfDouble);
            }
           MG->RestrictToAllGrids();
 
       break;*/
 
              
       case DIRECT:
     
              DirectSolver(SqmatrixA11[N_Levels-1], SqmatrixA12[N_Levels-1], SqmatrixA13[N_Levels-1], 
                           SqmatrixA21[N_Levels-1], SqmatrixA22[N_Levels-1], SqmatrixA23[N_Levels-1],  
                           SqmatrixA31[N_Levels-1], SqmatrixA32[N_Levels-1], SqmatrixA33[N_Levels-1],  
                           rhs, sol, 0);

       break;
          

       default:
             OutPut("Unknown Solver" << endl);
             exit(4711);;
      }   

 } 
