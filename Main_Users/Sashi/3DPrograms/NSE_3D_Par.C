// =======================================================================
//
// Purpose:     NSE3D main program with parallel direct solver (root will not take part in computation)
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 17.09.2010
// 
// =======================================================================
#ifdef _MPI
# include "mpi.h"
#endif

#ifdef _OMPONLY
#include <omp.h>
#endif

#include <Domain.h>
#include <Database.h>
#include <DirectSolver.h>
#include <DiscreteForm3D.h>
#include <FEDatabase3D.h>
#include <FESpace3D.h>
#include <FEVectFunct3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <AuxParam3D.h>
#include <Solver.h>
#include <Assemble3D.h>
#include <Output3D.h>
#include <LinAlg.h>
#include <NSE3D_ParamRout.h>
#include <BoundFace.h>
#include <MainUtilities.h>
#include <Upwind3D.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
// #include <malloc.h>
#include <math.h>

#ifdef _MPI
#include <MeshPartition.h>
#include <ParFECommunicator3D.h>

#include <NSE_ParSolver.h>
#include <NSE_ParSolver2.h>
// #include <MumpsSolver.h>
#include <ParVectorNSE3D.h>
#include <ParmsSolver.h>
#endif

#include <tetgen.h>
#include <GridCell.h>
#include <MacroCell.h>
#include <BdPlane.h>
#include <BdSphere.h>
#include <IsoBoundFace.h>
#include <InterfaceJoint3D.h>
#include <IsoInterfaceJoint3D.h>

double bound = 0;

#define AMG 0
#define GMG 1
#define DIRECT 2

#define MUMPS 101
#define PARMS 103
// =======================================================================
// include current example
// =======================================================================
// #include "../Examples/NSE_3D/DrivenCavity3D.h"
// #include "../Examples/NSE_3D/AnsatzLinConst.h"
//#include "Examples/NSE_3D/AnsatzQuadLin.h"
//#include "Examples/NSE_3D/BSExample_alpha.h"
#include "../Examples/NSE_3D/BSExample.h"
//#include "Examples/NSE3D/StatBSExample.h"
// #include "Examples/NSE_3D/Bench3DQuaderStat.h"
// #include "Examples/NSE_3D/Bench3DQuaderStatNeum.h"
// #include "Examples/NSE_3D/Bench3DCylinderStatNeum.h"
//#include "Examples/NSE_3D/Bench3DQuaderStatNeum.h"
//#include "Examples/NSE_3D/ChannelQuadLong.h"
//#include "Examples/NSE_3D/LES_Bench3DQuaderStatNeum.h"
//#include "../Examples/TNSE_3D/LES_Bench3DQuaderNeum.h"
//#include "../Examples/TNSE_3D/ChannelStepSlip3D.h"
//#include "../Examples/TNSE_3D/ChannelStepSlipDiri3D.h"
//#include "Examples/NSE_3D/BrennstoffzelleParabol.h"
//#include "Examples/NSE_3D/BrennstoffzelleParabolGravi.h"
//#include "Examples/NSE_3D/BrennStutzenNSE.h"
//#include "Examples/NSE_3D/NoSlipExample.h"
//#include "../Examples/TNSE_3D/ChannelSlipDiri3D.h"
//#include "../Examples/TNSE_3D/ChannelStepSlipTopFreeSlipLatDiri.h"
//#include "Examples/NSE_3D/NoFlow.h"
//#include "Examples/TNSE_3D/DrivenCavity3D_SFB_Periodic.h"
//#include "Examples/TNSE_3D/ChannelStepDiri3D.h"
//#include "../Examples/TNSE_3D/WallMountedCube.h"
//#include "../Examples/TNSE_3D/UnitCubeSlip3D.h"
//#include "../Examples/TNSE_3D/windchannel.00_bilin_inter.h"
//#include "../Examples/TNSE_3D/Calotte.h"
// #include "../Examples/NSE_3D/CircularChannel.h"
// #include "../Examples/TNSE_3D/Urea.h"
// #include "../Examples/NSE_3D/UShape3D.h"

int main(int argc, char* argv[])
{
#ifdef _MPI
  const int root = 0;
  int rank, size, len;
  double t_par1, t_par2, total_time;
  char  name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);

  MPI_Comm Comm = MPI_COMM_WORLD;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  MPI_Get_processor_name(name, &len);
#endif 

  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
#ifdef _MPI
  TDatabase::ParamDB->Comm = Comm;
  total_time = MPI_Wtime();
#else
 double total_time = GetTime();
#endif 

  TFEDatabase3D *FEDatabase = new TFEDatabase3D();
  TCollection *coll, *mortarcoll = NULL;
  TBaseCell *cell;
  TFESpace3D *velocity_space, *pressure_space;
  TFESpace3D **USpaces, **PSpaces;
  TOutput3D *Output;
  TAuxParam3D *aux;
  TFEVectFunct3D *u;
  TFEFunction3D *u1, *u2, *u3, *p, *fefct[7];
  TFESpace3D *fesp[4], *ferhs[3];

  TDiscreteForm3D *DiscreteFormGalerkin;
  TDiscreteForm3D *DiscreteFormSDFEM;
  TDiscreteForm3D *DiscreteFormUpwind;
  TDiscreteForm3D *DiscreteFormSmagorinsky;
  TDiscreteForm3D *DiscreteFormVMSProjection;

  TDiscreteForm3D *DiscreteFormNLGalerkin;
  TDiscreteForm3D *DiscreteFormNLSDFEM;
  TDiscreteForm3D *DiscreteFormNLUpwind;
  TDiscreteForm3D *DiscreteFormNLSmagorinsky;
  TDiscreteForm3D *DiscreteFormNLSDFEM_DivDiv;
  TDiscreteForm3D *DiscreteFormNLVMSProjection;

  TDiscreteForm3D *DiscreteFormPressSep;
  TDiscreteForm3D *DiscreteFormAuxProbPressSep;

  TDiscreteForm3D *DiscreteFormNSRFBRhs;

  TDiscreteForm3D *DiscreteForm, *DiscreteForm0;


  BoundCondFunct3D *BoundaryConditions[3];
  BoundValueFunct3D *BoundValues[3];

  TSquareStructure3D *sqstructureA;
  TStructure3D *structureB, *structureBT;
  TSquareMatrix3D *MatrixA, *SQMATRICES[9];
  TSquareMatrix3D *MatrixA11, *MatrixA12, *MatrixA13;
  TSquareMatrix3D *MatrixA21, *MatrixA22, *MatrixA23;
  TSquareMatrix3D *MatrixA31, *MatrixA32, *MatrixA33;
  TMatrix3D *MatrixB1, *MatrixB2, *MatrixB3, *MATRICES[9];
  TMatrix3D *MatrixB1T, *MatrixB2T, *MatrixB3T;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TMatrix **matrices = (TMatrix **)MATRICES;

  MatVecProc *MatVect;
  DefectProc *Defect;

  #ifdef _MPI
  int  out_rank, MaxCpV;
  int  MaxSubDomainPerDof;

  double l2, H1;

  bool  ActiveProcess;

  TNSE_ParSolver2 *Par_Solver;
  
  TParVectorNSE3D  *ParSolVect, *ParRhsNSEVect, *ParDefectVect;
  TParFECommunicator3D **ParComm = new TParFECommunicator3D*[2];

//   TMumpsSolver *MUMPS_NSESolver;
  TParmsSolver *Parms_NSESolver;
  ParDefectProc *ParDefect;
  #endif

//   bool Initialize_ScalarSolver = TRUE;

  char *PRM, *GEO;
  char *PsBaseName, *VtkBaseName;

  // strings
  char Readin[] = "readin.dat";
  char NameString[] = "name";
  char UString[] = "u";
  char PString[] = "p";
  char PsiString[] = "psi";
  char DivergenceString[] = "divergence";
  char SubID[] = "_";

  std::ostringstream os;
  os << " ";

  int i,j,k,l,m,n, N_, Len, low, N_Active, CurrentDiscType, N_BData=1;
  int N_Rows, N_Columns, N_U, N_P, N_Unknowns, N_V, N_Vort;
  int ret, pressure_space_code, velocity_space_code, N_Cells;
  int N_Rhs, N_FESpaces, N_SquareMatrices, N_RectMatrices;
  int *RowPtr, N_LinIter, N_SubSteps, time_discs, img=1;
  int Max_It, very_first_time=0, methods, last_sq, N_LinIterCurr;
  int N_NonActive, GlobalN_U, GlobalN_P;

  double RE_NR, t1, t2, t3, t4, solver_time=0., gamma;
  double *sol, *rhs, *oldsol, *sol_timestep_m1, *RHSs[4];
  double *defect, *defect_own, *startsol, *frac_step_sol, *oldrhs;
  double oldtau, end_time, hmin, hmax, tau;
  double tau2, tau1, limit, solver_time_curr, residual, impuls_residual, oldresidual;
  double errors[7], total_time1, total_time2;
  double l_infty_l_2 = 0, l_infty_l_2_time=-4711.0;
  double olderror = 0, l_2_l_2Du=0, l_2_l_2u=0 , olderror_l_2_l_2u=0;
  double l_2_h_1u=0, olderror_l_2_h_1u=0;


 double U1max=-1e8;
  //======================================================================
  // read parameter file
  //======================================================================

  if (argc >= 2)
    ret = Domain->ReadParam(argv[1]);
  else
    ret = Domain->ReadParam(Readin);

#ifdef  __BENCH3_CYLINDER__
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 1;
#endif

  if (ret == -1)
  {
    exit(-1);
  }

  OpenFiles();
  OutFile.setf(std::ios::scientific);

#ifndef _MPI
  TDatabase::ParamDB->SOLVER_TYPE = DIRECT;
  total_time1 = GetTime();
#else
  total_time1 = MPI_Wtime();
//   TDatabase::ParamDB->SOLVER_TYPE = PARMS;
//   Initialize_ScalarSolver = FALSE;
  out_rank=TDatabase::ParamDB->Par_P0;

  if(rank==out_rank)
#endif
   Database->WriteParamDB(argv[0]);

   Database->CheckParameterConsistencyNSE();
#ifndef __UREA__
#ifdef _MPI
  if(rank==TDatabase::ParamDB->Par_P0)
#endif
#endif
   ExampleFile();

  //======================================================================
  // copy read parameters into local variables
  //======================================================================

  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;

  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;

  //======================================================================
  // initialize discrete forms
  //======================================================================
  InitializeDiscreteForms(
    DiscreteFormGalerkin, DiscreteFormSDFEM,
    DiscreteFormUpwind, DiscreteFormSmagorinsky,
    DiscreteFormVMSProjection,
    DiscreteFormNLGalerkin, DiscreteFormNLSDFEM,
    DiscreteFormNLUpwind, DiscreteFormNLSmagorinsky, 
    DiscreteFormNLVMSProjection,
    DiscreteFormNLSDFEM_DivDiv,
    DiscreteFormPressSep,
    DiscreteFormAuxProbPressSep,
    DiscreteFormNSRFBRhs,
    LinCoeffs, TDatabase::ParamDB->NSTYPE);


  //======================================================================
  // read boundary parameterization and initialize coarse grid
  //======================================================================
  Domain->Init(PRM, GEO);

  
//    /** Using tetgen with smesh mesh */
//    TetrameshGen(Domain);  
  
//   TDatabase::ParamDB->UNIFORM_STEPS +=TDatabase::ParamDB->LEVELS;
  
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();
//   printf("rank %d MainPrg N_Neibs  %d   \n", rank, TDatabase::ParamDB->Par_P5);

  //======================================================================
  // Partition grid using Metis
  //======================================================================
#ifdef _MPI
  Domain->GenerateEdgeInfo();
 
  
  t_par1 = MPI_Wtime();
  Partition_Mesh3D(Comm, Domain, MaxCpV);
  t_par2 = MPI_Wtime();
 
  if(rank==0)
    printf("Time taken for Domain Decomposition is %e\n", (t_par2-t_par1));

  MaxSubDomainPerDof = MIN(MaxCpV, size);
#endif

 
  // *****************************************************************************
  // read boundary conditions and their values
  // *****************************************************************************
  BoundaryConditions[0] = BoundCondition;
  BoundaryConditions[1] = BoundCondition;
  BoundaryConditions[2] = BoundCondition;

  BoundValues[0] = U1BoundValue;
  BoundValues[1] = U2BoundValue;
  BoundValues[2] = U3BoundValue;

  
//======================================================================
// Generating FE spaces and allocating memory for their matrices
//======================================================================
#ifdef _MPI
   coll=Domain->GetOwnCollection(It_Finest, 0, rank);
#else
   coll=Domain->GetCollection(It_Finest, 0);
#endif

  N_Cells = coll->GetN_Cells();  
//   printf("Rank %d: N_Cells  : %d \n",rank, N_Cells);    
 

   Output = new TOutput3D(2, 2, 1, 1,Domain, coll);

   GetVelocityAndPressureSpace3D(coll,BoundCondition, velocity_space, pressure_space,
                                 &pressure_space_code, TDatabase::ParamDB->VELOCITY_SPACE,
                                 TDatabase::ParamDB->PRESSURE_SPACE);
/*    #ifdef _MPI     
    if(rank==0)
      printf("Rank %d: N_Cells  : %d \n",rank, N_Cells);    
      MPI_Finalize(); 
#endif  
     exit(0); */  
    velocity_space_code =   TDatabase::ParamDB->VELOCITY_SPACE;
    TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;

    N_U = velocity_space->GetN_DegreesOfFreedom();
    N_P = pressure_space->GetN_DegreesOfFreedom();
    N_Active = velocity_space->GetActiveBound();
    
# ifdef _MPI    
     i = TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE;     
     MPI_Allreduce(&i, &TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE,
                    1, MPI_INT, MPI_MIN, Comm);                           
     
    velocity_space->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
    pressure_space->SetMaxSubDomainPerDof(MaxSubDomainPerDof);


    t_par1 = MPI_Wtime();
    ParComm[0] = new TParFECommunicator3D(Comm, velocity_space);
    t_par2 = MPI_Wtime();
    if(rank==out_rank)
    printf("Time taken for Velo FEComm constructure: %e\n", (t_par2-t_par1));  
 
    t_par1 = MPI_Wtime();
    ParComm[1] = new TParFECommunicator3D(Comm, pressure_space);
    t_par2 = MPI_Wtime();
    if(rank==out_rank)
    printf("Time taken for Pressure FEComm constructure: %e\n", (t_par2-t_par1));
 
    GlobalN_U = ParComm[0]->GetN_GlobalDegreesOfFreedom();
    GlobalN_P = ParComm[1]->GetN_GlobalDegreesOfFreedom();
#else
    GlobalN_U = N_U;
    GlobalN_P = N_P;
#endif

    // build matrices
    structureB = new TStructure3D(pressure_space, velocity_space);
    structureBT = new TStructure3D(velocity_space, pressure_space);
    sqstructureA = new TSquareStructure3D(velocity_space);
    sqstructureA->Sort();// necessary for ParDirectSolver

      // allocate matrices
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          MatrixB1 = new TMatrix3D(structureB);
          MatrixB2 = new TMatrix3D(structureB);
          MatrixB3 = new TMatrix3D(structureB);

          MatrixA = new TSquareMatrix3D(sqstructureA);

          MatVect = MatVect_NSE1;
          Defect = Defect_NSE1;
          break;

        case 2:
          MatrixB1 = new TMatrix3D(structureB);
          MatrixB2 = new TMatrix3D(structureB);
          MatrixB3 = new TMatrix3D(structureB);
          MatrixB1T = new TMatrix3D(structureBT);
          MatrixB2T = new TMatrix3D(structureBT);
          MatrixB3T = new TMatrix3D(structureBT);

          MatrixA = new TSquareMatrix3D(sqstructureA);

          MatVect = MatVect_NSE2;
          Defect = Defect_NSE2;

# ifdef _MPI
          Par_Solver = new TNSE_ParSolver2(ParComm[0], sqstructureA,
                                ParComm[1], structureBT, structureB);

          ParDefect = ParDefect_NSE2;  
#else
	  Defect = Defect_NSE2;  
# endif
          break;

        case 3:
          MatrixB1 = new TMatrix3D(structureB);
          MatrixB2 = new TMatrix3D(structureB);
          MatrixB3 = new TMatrix3D(structureB);

          MatrixA11 = new TSquareMatrix3D(sqstructureA);
          MatrixA12 = new TSquareMatrix3D(sqstructureA);
          MatrixA13 = new TSquareMatrix3D(sqstructureA);
          MatrixA21 = new TSquareMatrix3D(sqstructureA);
          MatrixA22 = new TSquareMatrix3D(sqstructureA);
          MatrixA23 = new TSquareMatrix3D(sqstructureA);
          MatrixA31 = new TSquareMatrix3D(sqstructureA);
          MatrixA32 = new TSquareMatrix3D(sqstructureA);
          MatrixA33 = new TSquareMatrix3D(sqstructureA);

          MatVect = MatVect_NSE3;
          Defect = Defect_NSE3;
          break;

        case 4:
          MatrixB1 = new TMatrix3D(structureB);
          MatrixB2 = new TMatrix3D(structureB);
          MatrixB3 = new TMatrix3D(structureB);
          MatrixB1T = new TMatrix3D(structureBT);
          MatrixB2T = new TMatrix3D(structureBT);
          MatrixB3T = new TMatrix3D(structureBT);

          MatrixA11 = new TSquareMatrix3D(sqstructureA);
          MatrixA12 = new TSquareMatrix3D(sqstructureA);
          MatrixA13 = new TSquareMatrix3D(sqstructureA);
          MatrixA21 = new TSquareMatrix3D(sqstructureA);
          MatrixA22 = new TSquareMatrix3D(sqstructureA);
          MatrixA23 = new TSquareMatrix3D(sqstructureA);
          MatrixA31 = new TSquareMatrix3D(sqstructureA);
          MatrixA32 = new TSquareMatrix3D(sqstructureA);
          MatrixA33 = new TSquareMatrix3D(sqstructureA);

          MatVect = MatVect_NSE4;
          Defect = Defect_NSE4;

# ifdef _MPI
          Par_Solver = new TNSE_ParSolver2(ParComm[0], sqstructureA,
                                ParComm[1], structureBT, structureB);
# endif

          break;
      }

# ifdef _MPI
   if(rank==out_rank)
# endif
    {
      OutPut("dof velocity : "<< setw(10) << 3*GlobalN_U << endl);
      OutPut("dof pressure : "<< setw(10) << GlobalN_P << endl);
      OutPut("dof all      : "<< setw(10) << 3*GlobalN_U+GlobalN_P << endl);
   }

   N_Unknowns = 3*N_U + N_P;
   sol = new double[N_Unknowns];
   rhs = new double[N_Unknowns];
   defect  = new double [N_Unknowns];

   memset(sol, 0, N_Unknowns*SizeOfDouble);
   memset(rhs, 0, N_Unknowns*SizeOfDouble);     //working rhs
   memset(defect, 0, N_Unknowns*SizeOfDouble);
 
#ifdef _MPI
   //  velocity and pressure
   ParSolVect =  new TParVectorNSE3D(Comm, sol, N_U, N_P, 3, ParComm[0], ParComm[1]);
   ParRhsNSEVect =  new TParVectorNSE3D(Comm, rhs, N_U, N_P, 3, ParComm[0], ParComm[1]);
   ParDefectVect =  new TParVectorNSE3D(Comm, defect, N_U, N_P, 3, ParComm[0], ParComm[1]);
# endif

   u = new TFEVectFunct3D(velocity_space, UString, UString, sol, N_U, 3);
   u1 = u->GetComponent(0);
   u2 = u->GetComponent(1);
   u3 = u->GetComponent(2);

   p = new TFEFunction3D(pressure_space, PString, PString, sol+3*N_U, N_P);
   
   Output->AddFEVectFunct(u);
   Output->AddFEFunction(p);

   oldsol = new double[N_Unknowns];
   memset(oldsol, 0, N_Unknowns*SizeOfDouble);

    // Assemble the matrices
    // find discrete form
    switch(TDatabase::ParamDB->DISCTYPE)
       {
          case GALERKIN:
            DiscreteForm = DiscreteFormGalerkin;
          break;

          default:
            Error("Unknown DISCTYPE" << endl);
            return -1;
        }


      if (TDatabase::ParamDB->STOKES_PROBLEM)
        DiscreteForm = DiscreteFormUpwind;

      // initialize matrices
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          SQMATRICES[0] = MatrixA;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB3;

          SQMATRICES[0]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();

          N_SquareMatrices = 1;
          N_RectMatrices = 3;

          N_Rhs = 3;
          N_FESpaces = 2;

          break;

        case 2:
          SQMATRICES[0] = MatrixA;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB3;
          MATRICES[3] = MatrixB1T;
          MATRICES[4] = MatrixB2T;
          MATRICES[5] = MatrixB3T;

          SQMATRICES[0]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();
          MATRICES[4]->Reset();
          MATRICES[5]->Reset();

          N_SquareMatrices = 1;
          N_RectMatrices = 6;

          N_Rhs = 3;
          N_FESpaces = 2;
          break;

        case 3:
          SQMATRICES[0] = MatrixA11;
          SQMATRICES[1] = MatrixA12;
          SQMATRICES[2] = MatrixA13;
          SQMATRICES[3] = MatrixA21;
          SQMATRICES[4] = MatrixA22;
          SQMATRICES[5] = MatrixA23;
          SQMATRICES[6] = MatrixA31;
          SQMATRICES[7] = MatrixA32;
          SQMATRICES[8] = MatrixA33;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB3;

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          SQMATRICES[4]->Reset();
          SQMATRICES[5]->Reset();
          SQMATRICES[6]->Reset();
          SQMATRICES[7]->Reset();
          SQMATRICES[8]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();

          N_SquareMatrices = 9;
          N_RectMatrices = 3;

          N_Rhs = 3;
          N_FESpaces = 2;
          break;

        case 4:
          SQMATRICES[0] = MatrixA11;
          SQMATRICES[1] = MatrixA12;
          SQMATRICES[2] = MatrixA13;
          SQMATRICES[3] = MatrixA21;
          SQMATRICES[4] = MatrixA22;
          SQMATRICES[5] = MatrixA23;
          SQMATRICES[6] = MatrixA31;
          SQMATRICES[7] = MatrixA32;
          SQMATRICES[8] = MatrixA33;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB3;
          MATRICES[3] = MatrixB1T;
          MATRICES[4] = MatrixB2T;
          MATRICES[5] = MatrixB3T;

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          SQMATRICES[4]->Reset();
          SQMATRICES[5]->Reset();
          SQMATRICES[6]->Reset();
          SQMATRICES[7]->Reset();
          SQMATRICES[8]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();
          MATRICES[4]->Reset();
          MATRICES[5]->Reset();

          N_SquareMatrices = 9;
          N_RectMatrices = 6;

          N_Rhs = 3;
          N_FESpaces = 2;
          break;
      }

      N_NonActive = N_U - N_Active;

      RHSs[0] = rhs;
      RHSs[1] = rhs + N_U;
      RHSs[2] = rhs + 2*N_U;
      RHSs[3] = rhs + 3*N_U;
      memset(rhs, 0, (3*N_U+N_P)*SizeOfDouble);

      fesp[0] = velocity_space;
      fesp[1] = pressure_space;

      fefct[0] = u1;
      fefct[1] = u2;
      fefct[2] = u3;
      ferhs[0] = velocity_space;
      ferhs[1] = velocity_space;
      ferhs[2] = velocity_space;

       aux =  new TAuxParam3D(NSN_FESpacesVelo, NSN_FctVelo, NSN_ParamFctVelo,
          NSN_FEValuesVelo,
          fesp, fefct,
          NSFctVelo,
          NSFEFctIndexVelo, NSFEMultiIndexVelo,
          NSN_ParamsVelo, NSBeginParamVelo);

  
      // assemble
      Assemble3D(N_FESpaces, fesp,
        N_SquareMatrices, SQMATRICES,
        N_RectMatrices, MATRICES,
        N_Rhs, RHSs, ferhs,
        DiscreteForm,
        BoundaryConditions,
        BoundValues,
        aux);
        
#ifdef _MPI
     i = TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE;     
     MPI_Allreduce(&i, &TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE, 1, MPI_INT, MPI_MIN, Comm);  
#endif

      if ((DiscreteForm == DiscreteFormUpwind)
        &&(!TDatabase::ParamDB->STOKES_PROBLEM))
      {
        switch(TDatabase::ParamDB->NSTYPE)
        {
          case 1:
          case 2:
            // do upwinding with one matrix
            UpwindForNavierStokes3D(SQMATRICES[0], u1, u2, u3);
            cout << "UPWINDING DONE : " << endl;
            break;

          case 3:
          case 4:
            // do upwinding with two matrices
            UpwindForNavierStokes3D(SQMATRICES[0], u1, u2, u3);
            UpwindForNavierStokes3D(SQMATRICES[4], u1, u2, u3);
            UpwindForNavierStokes3D(SQMATRICES[8], u1, u2, u3);
            cout << "UPWINDING DONE :   "  << endl;
            break;
        }// endswitch
      } // endif

    // set rhs for Dirichlet nodes
    memcpy(sol+N_Active, rhs+N_Active, N_NonActive*SizeOfDouble);
    memcpy(sol+N_U+N_Active, rhs+N_U+N_Active, N_NonActive*SizeOfDouble);
    memcpy(sol+2*N_U+N_Active, rhs+2*N_U+N_Active, N_NonActive*SizeOfDouble);

 
    
//    u1->Interpolate(ExactU1);     
//    u2->Interpolate(ExactU2);     
//    u3->Interpolate(ExactU3);   
//    p->Interpolate(ExactP);      

#ifdef _MPI
    t_par1 = MPI_Wtime();
    if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(Comm, img, SubID);
      img++;
    t_par2 = MPI_Wtime();

    if(rank==out_rank)
     printf("Time taken for writing the parvtk file %e\n", (t_par2-t_par1));
#else

    if(TDatabase::ParamDB->WRITE_VTK)
     {
       os.seekp(std::ios::beg);
       if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
       else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
       else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
       else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
       else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());

       img++;
      }
#endif  
  
    // compute defect
    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
      IntoL20FEFunction3D(sol+3*N_U, N_P, pressure_space);

     memset(defect, 0, N_Unknowns*SizeOfDouble);
    
#ifdef _MPI
     ParDefect(sqmatrices, matrices, ParSolVect, ParRhsNSEVect, ParDefectVect);             
     ParDefectVect->ParDdot(BYOWN, residual, impuls_residual);
#else
    Defect(sqmatrices, matrices, sol, rhs, defect);
    residual =  Ddot(N_Unknowns,defect,defect);
    impuls_residual = Ddot(3*N_U,defect,defect);
#endif

#ifdef _MPI
    if(rank==out_rank)
#endif
    {
    OutPut("Nonlinear iteration step   0");
    OutPut(setw(14) << impuls_residual);
    OutPut(setw(14) << residual-impuls_residual);
    OutPut(setw(14) << sqrt(residual) << endl);
    }


     // solve the system
#ifndef _MPI
    t1 = GetTime();
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
         t1 = GetTime();

         switch(TDatabase::ParamDB->NSTYPE)
          {
           case 2:
            DirectSolver(MatrixA, MatrixB1T, MatrixB2T, MatrixB3T, 
                         MatrixB1, MatrixB2, MatrixB3, rhs, sol);
          
           break;
           case 4:    
            DirectSolver(MatrixA11, MatrixA12, MatrixA13,
			 MatrixA21, MatrixA22, MatrixA23,
			 MatrixA31, MatrixA32, MatrixA33,
			 MatrixB1T, MatrixB2T, MatrixB3T, 
                         MatrixB1, MatrixB2, MatrixB3, rhs, sol, 3);      
           break;	     
           default:
            OutPut("Parallel solver for NSTYPE " << TDatabase::ParamDB->NSTYPE
              << " not implemented !!!" << endl);
           exit(4711);
          }  
         t2 = GetTime();
       break;

#ifdef _OMP
       case 100:
             // pardiso
             cout << "solver type not implemented yet !!!!!!!!!!!!!" << endl;
             exit(0);

       break;
#endif
       default: 
        cout << "wrong  solver type !!!!!!!!!!!!!" << endl;
        exit(0);
       break;
      }
    t2 = GetTime();


#else // Parallel solver
    t1 = MPI_Wtime();
    switch(TDatabase::ParamDB->NSTYPE)
     {
      case 2:

        Par_Solver->Solve(MatrixA, MatrixB1T, MatrixB2T, MatrixB3T,
                          MatrixB1, MatrixB2, MatrixB3, ParRhsNSEVect, ParSolVect);

      break;
      default:
       OutPut("Parallel solver for NSTYPE " << TDatabase::ParamDB->NSTYPE
              << " not implemented !!!" << endl);
       exit(4711);
      }
    t2 = MPI_Wtime();
#endif

   solver_time += (t2-t1);

#ifdef _MPI
     if(rank==TDatabase::ParamDB->Par_P0)
#endif
     {
      OutPut("Time taken by SOLVER : " << solver_time << "s" << endl);
      OutPut("MEMORY: " << setw(10) << GetMemory()/(1048576.0));
      OutPut(" MB" << endl);
     }

     
     
     
     

//     // set rhs for Dirichlet nodes
//     memcpy(sol+N_Active, rhs+N_Active, N_NonActive*SizeOfDouble);
//     memcpy(sol+N_U+N_Active, rhs+N_U+N_Active, N_NonActive*SizeOfDouble);
//     memcpy(sol+2*N_U+N_Active, rhs+2*N_U+N_Active, N_NonActive*SizeOfDouble);
// 
//   

  // ************************************************************* //
  // end of first nonlinear step
  // the nonlinear iteration
  // ************************************************************* //
  for(j=1;j<=Max_It;j++)
   {
    // assemble the non-linear matrix
    switch(TDatabase::ParamDB->DISCTYPE)
     {
      case GALERKIN:
        DiscreteForm = DiscreteFormNLGalerkin;
      break;

      case UPWIND:
        DiscreteForm = DiscreteFormNLUpwind;
      break;

      default:
        Error("Unknown DISCTYPE" << endl);
       return -1;
     }                      // endswitch

    switch (TDatabase::ParamDB->NSTYPE)
     {
       case 1:
        SQMATRICES[0] = MatrixA;
        SQMATRICES[0]->Reset();

        N_SquareMatrices = 1;
        N_RectMatrices = 0;

        N_Rhs = 0;
        N_FESpaces = 1;
       break;

       case 2:
        SQMATRICES[0] = MatrixA;
        SQMATRICES[0]->Reset();

        N_SquareMatrices = 1;
        N_RectMatrices = 0;
        N_Rhs = 0;
        N_FESpaces = 1;
       break;

       case 3:
        SQMATRICES[0] = MatrixA11;
        SQMATRICES[1] = MatrixA22;
        SQMATRICES[2] = MatrixA33;
        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        SQMATRICES[2]->Reset();

        N_SquareMatrices = 3;
        N_RectMatrices = 0;
        last_sq = 2;
        N_Rhs = 0;
        N_FESpaces = 1;
       break;

       case 4:
        SQMATRICES[0] = MatrixA11;
        SQMATRICES[1] = MatrixA22;
        SQMATRICES[2] = MatrixA33;
        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        SQMATRICES[2]->Reset();

        N_SquareMatrices = 3;
        N_RectMatrices = 0;
        last_sq = 2;
        N_Rhs = 0;
        N_FESpaces = 1;
       break;

     } //  switch (TDatabase::ParamDB->NSTYPE)

    fesp[0] = velocity_space;
    fesp[1] = pressure_space;

    fefct[0] = u1;
    fefct[1] = u2;
    fefct[2] = u3;

    aux =  new TAuxParam3D(NSN_FESpacesVelo, NSN_FctVelo, NSN_ParamFctVelo,
            NSN_FEValuesVelo,
            fesp, fefct,
            NSFctVelo,
            NSFEFctIndexVelo, NSFEMultiIndexVelo,
            NSN_ParamsVelo, NSBeginParamVelo);

     // assembling
     Assemble3D(N_FESpaces, fesp,
          N_SquareMatrices, SQMATRICES,
          N_RectMatrices, MATRICES,
          N_Rhs, RHSs, ferhs,
          DiscreteForm,
          BoundaryConditions,
          BoundValues,
          aux);
     
#ifdef _MPI // due to BC Cond it has to be set after every assemble
     i = TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE;     
     MPI_Allreduce(&i, &TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE, 1, MPI_INT, MPI_MIN, Comm);  
#endif
     
     delete aux;

      if(DiscreteForm == DiscreteFormNLUpwind)
      {
        switch(TDatabase::ParamDB->NSTYPE)
        {
          case 1:
          case 2:
            // do upwinding with one matrix
            UpwindForNavierStokes3D(SQMATRICES[0], u1, u2, u3);
            cout << "UPWINDING DONE : " << endl;
            break;

          case 3:
          case 4:
            // do upwinding with two matrices
            UpwindForNavierStokes3D(SQMATRICES[0], u1, u2, u3);
            UpwindForNavierStokes3D(SQMATRICES[1], u1, u2, u3);
            UpwindForNavierStokes3D(SQMATRICES[2], u1, u2, u3);
            cout << "UPWINDING DONE :   "  << endl;
            break;
        }// endswitch
      } // endif

    // end of assembling
    memcpy(sol+N_Active, rhs+N_Active, N_NonActive*SizeOfDouble);
    memcpy(sol+N_U+N_Active, rhs+N_U+N_Active, N_NonActive*SizeOfDouble);
    memcpy(sol+2*N_U+N_Active, rhs+2*N_U+N_Active, N_NonActive*SizeOfDouble);


    // compute defect
    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
      IntoL20FEFunction3D(sol+3*N_U, N_P, pressure_space);

     memset(defect, 0, N_Unknowns*SizeOfDouble);
    
#ifdef _MPI
     ParDefect(sqmatrices, matrices, ParSolVect, ParRhsNSEVect, ParDefectVect);             
     ParDefectVect->ParDdot(BYOWN, residual, impuls_residual);
#else
    Defect(sqmatrices, matrices, sol, rhs, defect);
    residual =  Ddot(N_Unknowns,defect,defect);
    impuls_residual = Ddot(3*N_U,defect,defect);
#endif

#ifdef _MPI
    if(rank==out_rank)
#endif
//     {
//      OutPut("nonlinear iteration step " << setw(3) << j);
//      OutPut(setw(14) << impuls_residual);
//      OutPut(setw(14) << residual-impuls_residual);
//      OutPut(setw(14) << sqrt(residual) << endl);
//     }

 
/*#ifdef _MPI
  MPI_Finalize();
#endif
  exit(0);  */   
       
    
    
   if ((sqrt(residual)<=limit)||(j==Max_It))
    {
#ifndef _MPI
     total_time2 = GetTime();
#else
     total_time2 = MPI_Wtime();
     if(rank==out_rank)
#endif
      {
       OutPut("ITE : " << setw(3) << j);
       OutPut(" T/SOLVER : " << solver_time << "s");
       OutPut(" T/TOTAL : " << total_time2-total_time1 << "s");
       OutPut("("<<setprecision(2)<<solver_time/(total_time2-total_time1)<<")");
       OutPut(setprecision(6)<<" RES : " <<  sqrt(residual));
       OutPut(endl);
      }

      break;
     }

  
     
     // solve the system with a Direct seq/parallel solver
#ifndef _MPI
    t1 = GetTime();
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
         t1 = GetTime();
             DirectSolver(MatrixA, MatrixB1T, MatrixB2T, MatrixB3T, 
                          MatrixB1, MatrixB2, MatrixB3, rhs, sol);
         t2 = GetTime();
       break;

#ifdef _OMP
       case 100:
             // pardiso
             cout << "solver type not implemented yet !!!!!!!!!!!!!" << endl;
             exit(0);

       break;
#endif
       default: 
        cout << "wrong  solver type !!!!!!!!!!!!!" << endl;
        exit(0);
       break;
      }
    t2 = GetTime();
#else // Parallel solver
    t1 = MPI_Wtime();
    switch(TDatabase::ParamDB->NSTYPE)
     {
      case 2:
        Par_Solver->Solve(MatrixA, ParRhsNSEVect, ParSolVect);
      break;
      default:
       OutPut("Parallel solver for NSTYPE " << TDatabase::ParamDB->NSTYPE
              << " not implemented !!!" << endl);
       exit(4711);
      }
    t2 = MPI_Wtime();
#endif
 
     solver_time += (t2-t1);
     
#ifdef _MPI
     ParSolVect->ParDdot(BYOWN, residual, impuls_residual);     
#else
    Defect(sqmatrices, matrices, sol, rhs, defect);
    residual =  Ddot(N_Unknowns,sol,sol);
    impuls_residual = Ddot(3*N_U,sol,sol);
#endif

#ifdef _MPI
    if(rank==out_rank)
#endif
    {
     OutPut("nonlinear iteration step " << setw(3) << j);
     OutPut(setw(14) << impuls_residual);
     OutPut(setw(14) << residual-impuls_residual);
     OutPut(setw(14) << sqrt(residual) << endl);
    }

#ifdef _MPI
  MPI_Finalize();
#endif
//   exit(0);   
     
     
     
   } //for(j=1;j<=Max_It;j++) 

   
   
  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    IntoL20FEFunction3D(sol+3*N_U, N_P, pressure_space);

#ifdef _MPI
    t_par1 = MPI_Wtime();
    if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(Comm, img, SubID);
      img++;
    t_par2 = MPI_Wtime();

    if(rank==out_rank)
     printf("Time taken for writing the parvtk file %e\n", (t_par2-t_par1));
#else

    if(TDatabase::ParamDB->WRITE_VTK)
     {
       os.seekp(std::ios::beg);
       if(img<10) os << VtkBaseName<<".0000"<<img<<".vtk" << ends;
       else if(img<100) os << VtkBaseName<<".000"<<img<<".vtk" << ends;
       else if(img<1000) os << VtkBaseName<<".00"<<img<<".vtk" << ends;
       else if(img<10000) os << VtkBaseName<<".0"<<img<<".vtk" << ends;
       else  os << VtkBaseName<<"."<<img<<".vtk" << ends;
       Output->WriteVtk(os.str().c_str());

       img++;
      }
#endif

    // write solution for later use
//     u->WriteSol(0);
    
// #ifdef _MPI
//  if(rank==TDatabase::ParamDB->Par_P0)
//  printf("Rank  %d main  \n",rank );
//  MPI_Finalize();
// #else
//  cout << "test main " << endl;
// #endif
// 
//  exit(0);

#ifdef _MPI
  MPI_Finalize();
#endif
  CloseFiles();

  return 0;
}








