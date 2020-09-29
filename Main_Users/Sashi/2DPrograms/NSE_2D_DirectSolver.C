// =======================================================================
//
// Purpose:     NSE2D main program with parallel solver (root will not take part in computation)
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 12.10.2009
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
#include <DiscreteForm2D.h>

#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <AuxParam2D.h>
#include <QuadAffin.h>
#include <Solver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <Matrix2D.h>
#include <NSE2D_ParamRout.h>
#include <TNSE2D_ParamRout.h>
#include <VMS.h>
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include <MainUtilities.h>
#include <Upwind.h>
#include <DirectSolver.h>
#include <BoundEdge.h>
#include <QuadBilinear.h>

#ifdef _MPI
 #include <MeshPartition.h>
 #include <ParFECommunicator2D.h>
 #include <MumpsSolver.h>
#include <ParVectorNSE.h>
#endif


double bound = 0;
// =======================================================================
// include current example
// =======================================================================
// #include "Examples/NSE_2D/Linear.h"
// #include "Examples/NSE_2D/Const.h"
// #include "Examples/NSE_2D/FSHabil.h"
//#include "Examples/NSE_2D/FSHabil_slip.h"
// #include "Examples/NSE_2D/DC2.h"
// #include "Examples/NSE_2D/DrivenCavity.h"
#include "../Examples/NSE_2D/Benchmark.h"
// #include "Examples/NSE_2D/Benchmark_Neum.h"
// #include "Examples/NSE_2D/Frequence.h"
// #include "Examples/NSE_2D/Poiseuille.h"
// #include "Examples/NSE_2D/Poiseuille_Neum.h"
// #include "Examples/NSE_2D/Poiseuille2.h"
// #include "Examples/NSE_2D/Einfach.h"
// #include "Examples/NSE_2D/SinCos.h"
// #include "Examples/NSE_2D/BraessSarazin.h"
// #include "Examples/NSE_2D/BraessSarazinGeneralized.h"
// #include "Examples/NSE_2D/ZeroSolution.h"
//#include "Examples/NSE_2D/FerroMagnet.h"
//#include "../Examples/TNSE_2D/ChannelStep.h"
//#include "../Examples/NSE_2D/PressureGradient.h"
//#include "../Examples/NSE_2D/constant_velo_free_slip.h"
//#include "Examples/TNSE_2D/TrapezCavity.h"
//#include "Examples/TNSE_2D/StepCavity.h"
//#include "../Examples/TNSE_2D/ChannelStep.h"
//#include "Examples/NSE_2D/BenchmarkQuadNeum.h"
//#include "Examples/NSE_2D/atife_diri.h"
//#include "Examples/NSE_2D/atife_slip.h"
//#include "Examples/NSE_2D/Channel30_diri.h"
//#include "Examples/NSE_2D/Channel30_slip.h"
//#include "Examples/NSE_2D/atife_diri_sincos.h"
//#include "Examples/NSE_2D/atife_slip_sincos.h"
//#include "Examples/NSE_2D/Brennstoffzelle.h"
//#include "BOSCH/data/bosch_0433175329_00.h"
//#include "Examples/NSE_2D/NonSymmUnitSq.h"
//#include "Examples/NSE_2D/ZeroVeloLinPress.h"
//#include "Examples/NSE_2D/GleitRB_HalfCircle.h"
//#include "Examples/NSE_2D/Example_sashikumar.h"
//#include "Examples/NSE_2D/Example_sashikumar.01.h"
//#include "Examples/NSE_2D/Example_sashikumar.02.h"
//#include "Examples/NSE_2D/Example_sashikumar.03.h"
//#include "Examples/NSE_2D/NoFlow.h"
//#include "Examples/TNSE_2D/UnitSquareHoleDrag.h"
//#include "Examples/TNSE_2D/SFBCavity.h"
//#include "Examples/TNSE_2D/SFBCavity_2.h"
//#include "Examples/TNSE_2D/SFBCavity_3.h"
//#include "Examples/NSE_2D/Gartling2.h"
//#include "Examples/NSE_2D/CTVW05.h"
// #include "Examples/NSE_2D/SinCos.h"
//#include "Examples/NSE_2D/Calotte_test.h"
//#include "Examples/NSE_2D/LinkeNoFlow.h"
//  #include "Examples/NSE_2D/SimPaTurS.h"
// #include "Examples/NSE_2D/BraessSarazinNSE.h"
// #include "Examples/NSE_2D/CirculateNSE.h"
// #include "Examples/NSE_2D/Benchmark_FSI.h"

// ======================================================================
// utilities for main program
// ======================================================================

int main(int argc, char* argv[])
{
  #ifdef _MPI
  const int root = 0;
  int rank, size, len;
  double t_par1, t_par2;
  char  name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Get_processor_name(name, &len);
  #endif 

  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *coll, *mortarcoll = NULL;
  TBaseCell *cell;
  TFESpace2D *velocity_space, *pressure_space;
  TFESpace2D **USpaces, **PSpaces;
  TOutput2D *Output;
  TAuxParam2D *aux;
  TFEVectFunct2D *u;
  TFEFunction2D *u1, *u2, *p, *fefct[7];
  TFESpace2D *fesp[4], *ferhs[3];

  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormSDFEM;
  TDiscreteForm2D *DiscreteFormUpwind;
  TDiscreteForm2D *DiscreteFormSmagorinsky;
  TDiscreteForm2D *DiscreteFormVMSProjection;

  TDiscreteForm2D *DiscreteFormNLGalerkin;
  TDiscreteForm2D *DiscreteFormNLSDFEM;
  TDiscreteForm2D *DiscreteFormNLUpwind;
  TDiscreteForm2D *DiscreteFormNLSmagorinsky;
  TDiscreteForm2D *DiscreteFormNLVMSProjection;

  TDiscreteForm2D *DiscreteFormPressSep;
  TDiscreteForm2D *DiscreteFormAuxProbPressSep;

  TDiscreteForm2D *DiscreteFormNSRFBRhs;

  TDiscreteForm2D *DiscreteForm, *DiscreteForm0;

  BoundCondFunct2D *BoundaryConditions[2];
  BoundValueFunct2D *BoundValues[2];

  #ifdef _MPI
  int  out_rank, MaxCpV;
  int  MaxSubDomainPerDof;

  double l2, H1;

  TMumpsSolver *MUMPS_NSESolver;
  TParVectorNSE  *ParSolVect, *ParRhsNSEVect, *ParDefectVect;
  TParFECommunicator2D **ParComm = new TParFECommunicator2D*[2];
  ParDefectProc *ParDefect;
  #endif

  TStructure2D *structureB, *structureBT;
  TSquareStructure2D *sqstructureA;
  TMatrix2D *MatrixB1, *MatrixB2, *MatrixB1T , *MatrixB2T,  *MATRICES[10];
  TSquareMatrix2D *SqmatrixA, *SqmatrixM, *SQMATRICES[8];
  TSquareMatrix2D *SqmatrixA11, *SqmatrixA12, *SqmatrixA21, *SqmatrixA22;
  TSquareMatrix2D *SqmatrixM11, *SqmatrixM12, *SqmatrixM21, *SqmatrixM22;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TMatrix **matrices = (TMatrix **)MATRICES;
  MatVecProc *MatVect;
  DefectProc *Defect;
  CoeffFct2D *Coefficients[1];

  bool Initialize_ScalarSolver = TRUE;

  char *PRM, *GEO;
  char *PsBaseName, *VtkBaseName;
  // strings
  char ReadinDat[] = "readin.dat";
  char NameString[] = "name";
  char UString[] = "u";
  char PString[] = "p";
  char PsiString[] = "psi";
  char DivergenceString[] = "divergence";
  char *SubID = NULL;
  
  std::ostringstream os;
  os << " ";

  int i,j,k,l,m,n, N_, Len, low, N_Active, CurrentDiscType, N_BData=1;
  int N_Rows, N_Columns, N_U, N_P, N_Unknowns, N_V, N_Vort;
  int ret, pressure_space_code, velocity_space_code;
  int N_Rhs, N_FESpaces, N_SquareMatrices, N_RectMatrices;
  int *RowPtr, N_LinIter, N_SubSteps, time_discs, img=1;
  int Max_It, very_first_time=0, methods, last_sq, N_LinIterCurr;
  int N_NonActive;

  double RE_NR, total_time, t1, t2, t3, t4, solver_time=0., gamma;
  double *sol, *rhs, *oldsol, *sol_timestep_m1, *RHSs[4];
  double *defect, *defect_own, *startsol, *frac_step_sol, *oldrhs;
  double oldtau, end_time, hmin, hmax, tau;
  double tau2, tau1, limit, solver_time_curr, residual, impuls_residual, oldresidual;
  double errors[7], total_time1, total_time2;
  double l_infty_l_2 = 0, l_infty_l_2_time=-4711.0;
  double olderror = 0, l_2_l_2Du=0, l_2_l_2u=0 , olderror_l_2_l_2u=0;
  double l_2_h_1u=0, olderror_l_2_h_1u=0;
  //======================================================================
  // read parameter file
  //======================================================================
  total_time = GetTime();
  if(argc>=2)
    ret=Domain->ReadParam(argv[1]);
  else
    ret=Domain->ReadParam(ReadinDat);

  RE_NR=TDatabase::ParamDB->RE_NR;

  if(ret==-1)
  {
    exit(-1);
  }

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->CheckParameterConsistencyNSE();
  #ifdef _MPI
  Initialize_ScalarSolver = FALSE;
  out_rank=TDatabase::ParamDB->Par_P0;
  TDatabase::ParamDB->Comm = MPI_COMM_WORLD;

  if(rank==out_rank)
  #endif
  {
   Database->WriteParamDB(argv[0]);
   Database->WriteTimeDB();
  }

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
    DiscreteFormPressSep,
    DiscreteFormAuxProbPressSep,
    DiscreteFormNSRFBRhs,
    LinCoeffs, TDatabase::ParamDB->NSTYPE);

   Coefficients[0] = LinCoeffs;
/*    printf("Test  Decomposition \n");
   MPI_Finalize();
   exit(0);*/ 
  //======================================================================
  // read boundary parameterization and initialize coarse grid
  //======================================================================
  Domain->Init(PRM, GEO);

  /** MEdit Mesh*/
//   MEditMeshCreate(Domain);
    
  
//   Domain->ConvertQuadToTri(TDatabase::ParamDB->CONVERT_QUAD_TO_TRI);
  # ifdef __SIMPATURS__
  for(i=0;i<TDatabase::ParamDB->P9;i++)
    Domain->RefineallxDirection();  
  # endif  

  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();

  # ifdef _MPI
  if(rank==out_rank)
  # endif
  {
   if(TDatabase::ParamDB->WRITE_PS)
    {
      // write grid into an Postscript file
     os.seekp(std::ios::beg);
     os << "Domain" << ".ps" << ends;
     Domain->PS(os.str().c_str(),It_Finest,0);
   }
  }
exit(0);

  //======================================================================
  // Partition grid using Metis
  //======================================================================
  #ifdef _MPI
  bool  ActiveProcess = TDatabase::ParamDB->ActiveProcess;

  t_par1 = MPI_Wtime();
  Partition_Mesh(MPI_COMM_WORLD, Domain, MaxCpV);
  t_par2 = MPI_Wtime();

  if(rank==out_rank)
   printf("Time taken for Domain Decomposition is %e\n", (t_par2-t_par1));
  MaxSubDomainPerDof = MIN(MaxCpV, size);
  #endif //   _MPI


  # ifdef _MPI
  if(rank!=out_rank)
  # endif
  {
   if(TDatabase::ParamDB->WRITE_PS)
    {
      // write grid into an Postscript file
     os.seekp(std::ios::beg);
     os << "Domain" << ".ps" << ends;
     Domain->PS(os.str().c_str(),It_Finest,0);
   }
  }



  // *****************************************************************************
  // read boundary conditions and their values
  // *****************************************************************************
  BoundaryConditions[0] = BoundCondition;
  BoundaryConditions[1] = BoundCondition;

  BoundValues[0] = U1BoundValue;
  BoundValues[1] = U2BoundValue;

  t3 = GetTime();
  total_time = t3 - total_time;

//======================================================================
// Generating FE spaces and allocating memory for their matrices
//======================================================================
    #ifdef _MPI
    total_time1 = MPI_Wtime();
    #else
    total_time1 = GetTime();
    #endif

    coll=Domain->GetCollection(It_Finest, 0);
    Output = new TOutput2D(2, 2, 1, 1,Domain);
    cout << endl << endl;
    OutPut("N_Cells : "<<   coll->GetN_Cells() << endl);

    GetVelocityAndPressureSpace(coll,BoundCondition, mortarcoll, velocity_space,
                                pressure_space, &pressure_space_code,
                                TDatabase::ParamDB->VELOCITY_SPACE,
                                TDatabase::ParamDB->PRESSURE_SPACE);

    velocity_space_code =   TDatabase::ParamDB->VELOCITY_SPACE;
    TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;

    N_U = velocity_space->GetN_DegreesOfFreedom();
    N_P = pressure_space->GetN_DegreesOfFreedom();
    N_Active = velocity_space->GetActiveBound();

    // build matrices
    structureB = new TStructure2D(pressure_space, velocity_space);
    structureBT = new TStructure2D(velocity_space, pressure_space);
    sqstructureA = new TSquareStructure2D(velocity_space);

    sqstructureA->Sort();// necessary for ParDirectSolver

    # ifdef _MPI
    velocity_space->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
    pressure_space->SetMaxSubDomainPerDof(MaxSubDomainPerDof);

    t_par1 = MPI_Wtime();
    ParComm[0] = new TParFECommunicator2D(MPI_COMM_WORLD, velocity_space);
    ParComm[1] = new TParFECommunicator2D(MPI_COMM_WORLD, pressure_space);
    t_par2 = MPI_Wtime();
    if(rank==out_rank)
    printf("Time taken for Velo and Pressure dof mapping %e\n", (t_par2-t_par1));

    if(ActiveProcess)
    # endif
    {
    switch(TDatabase::ParamDB->NSTYPE)
     {
      case 1:
        MatrixB1 = new TMatrix2D(structureB);
        MatrixB2 = new TMatrix2D(structureB);

        SqmatrixA = new TSquareMatrix2D(sqstructureA);
        MatVect = MatVect_NSE1;
//         Defect = Defect_NSE1;
      break;

      case 2:
        MatrixB1 = new TMatrix2D(structureB);
        MatrixB2 = new TMatrix2D(structureB);
        MatrixB1T = new TMatrix2D(structureBT);
        MatrixB2T = new TMatrix2D(structureBT);

        SqmatrixA = new TSquareMatrix2D(sqstructureA);
        MatVect = MatVect_NSE2;

        #ifdef _MPI
        ParDefect = ParDefect_NSE2;
        #else
        Defect = Defect_NSE2;
        #endif

      break;

      case 3:
        MatrixB1 = new TMatrix2D(structureB);
        MatrixB2 = new TMatrix2D(structureB);

        SqmatrixA11 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA12 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA21 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA22 = new TSquareMatrix2D(sqstructureA);

        MatVect = MatVect_NSE3;
//         Defect = Defect_NSE3;
      break;

      case 4:
        MatrixB1 = new TMatrix2D(structureB);
        MatrixB2 = new TMatrix2D(structureB);
        MatrixB1T = new TMatrix2D(structureBT);
        MatrixB2T = new TMatrix2D(structureBT);

        SqmatrixA11 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA12 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA21 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA22 = new TSquareMatrix2D(sqstructureA);

        MatVect = MatVect_NSE4;
//         Defect = Defect_NSE4;
      break;
     }
    }

   N_Unknowns = 2*N_U + N_P;

   # ifdef _MPI
   if(rank==0)
   # endif
    {
    OutPut("dof velocity : "<< setw(10) << 2* N_U << endl);
    OutPut("dof pressure : "<< setw(10) << N_P << endl);
    OutPut("dof all      : "<<  setw(10) << N_Unknowns  << endl);
   }

   sol = new double[N_Unknowns];
   rhs = new double[N_Unknowns];
   defect  = new double [N_Unknowns];

   memset(sol, 0, N_Unknowns*SizeOfDouble);
   memset(rhs, 0, N_Unknowns*SizeOfDouble);     //working rhs
   memset(defect, 0, N_Unknowns*SizeOfDouble);  

   #ifdef _MPI
   //  velocity and pressure
   ParSolVect =  new TParVectorNSE(MPI_COMM_WORLD, sol, N_U, N_P, 2, ParComm[0], ParComm[1]);
   ParRhsNSEVect =  new TParVectorNSE(MPI_COMM_WORLD, rhs, N_U, N_P, 2, ParComm[0], ParComm[1]);
   ParDefectVect =  new TParVectorNSE(MPI_COMM_WORLD, defect, N_U, N_P, 2, ParComm[0], ParComm[1]);
   # endif


   u = new TFEVectFunct2D(velocity_space, UString,  UString,  sol, N_U, 2);
   u1 = u->GetComponent(0);
   u2 = u->GetComponent(1);
   p = new TFEFunction2D(pressure_space, PString,  PString,  sol+2*N_U, N_P);

   Output->AddFEVectFunct(u);
   Output->AddFEFunction(p);


//   u1->Interpolate(ExactU1);
//   u2->Interpolate(ExactU2);
//   p->Interpolate(ExactP);


   #ifdef _MPI
   if(ActiveProcess)
   #endif
   {
    oldsol = new double[N_Unknowns];
    memset(oldsol, 0, N_Unknowns*SizeOfDouble);

    // Assemble the matrices
    // find discrete form
    switch(TDatabase::ParamDB->DISCTYPE)
       {
          case GALERKIN:
            DiscreteForm = DiscreteFormGalerkin;
          break;

//           case SDFEM:
//             DiscreteForm = DiscreteFormSDFEM;
// 	    break;
// 
//           case UPWIND:
//             DiscreteForm = DiscreteFormUpwind;
//             break;
// 
//           case SMAGORINSKY:
//             DiscreteForm = DiscreteFormSmagorinsky;
//             break;

//           case VMS_PROJECTION:
// 	      DiscreteForm = DiscreteFormVMSProjection;
// 	      if (TDatabase::ParamDB->NSTYPE != 1)
// 	    {
// 		OutPut("VMS only for NSTYPE 1 implemented !!!"<<endl);
// 		exit(4711);
// 	    }
//             break;

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
          SQMATRICES[0] = SqmatrixA;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;

          SQMATRICES[0]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();

          N_SquareMatrices = 1;
          N_RectMatrices = 2;

          N_Rhs = 2;
          N_FESpaces = 2;
          break;

        case 2:
          SQMATRICES[0] = SqmatrixA;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;

          SQMATRICES[0]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();

          N_SquareMatrices = 1;
          N_RectMatrices = 4;

          N_Rhs = 2;
          N_FESpaces = 2;
          break;

        case 3:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[1] = SqmatrixA12;
          SQMATRICES[2] = SqmatrixA21;
          SQMATRICES[3] = SqmatrixA22;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();

          N_SquareMatrices = 4;
          N_RectMatrices = 2;

          N_Rhs = 2;
          N_FESpaces = 2;
          break;

        case 4:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[1] = SqmatrixA12;
          SQMATRICES[2] = SqmatrixA21;
          SQMATRICES[3] = SqmatrixA22;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();

          N_SquareMatrices = 4;
          N_RectMatrices = 4;

          N_Rhs = 2;
          N_FESpaces = 2;
          break;

      } //  switch(TDatabase::ParamDB->NSTYPE)


      N_NonActive = N_U - N_Active;

      RHSs[0] = rhs;
      RHSs[1] = rhs + N_U;
      RHSs[2] = rhs + 2*N_U;
      memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);

      fesp[0] = velocity_space;
      fesp[1] = pressure_space;

      fefct[0] = u1;
      fefct[1] = u2;
      ferhs[0] = velocity_space;
      ferhs[1] = velocity_space;
      ferhs[2] = pressure_space;

     // get auxiliary values
     aux =  new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo, NSN_ParamFctVelo,
                            NSN_FEValuesVelo,
                            fesp, fefct,
                            NSFctVelo,
                            NSFEFctIndexVelo, NSFEMultiIndexVelo,
                            NSN_ParamsVelo, NSBeginParamVelo);

      // assemble
      Assemble2D(N_FESpaces, fesp,
        N_SquareMatrices, SQMATRICES,
        N_RectMatrices, MATRICES,
        N_Rhs, RHSs, ferhs,
        DiscreteForm,
        BoundaryConditions,
        BoundValues,
        aux);


      if ((DiscreteForm == DiscreteFormUpwind)
        &&(!TDatabase::ParamDB->STOKES_PROBLEM))
      {
        switch(TDatabase::ParamDB->NSTYPE)
        {
          case 1:
          case 2:
            // do upwinding with one matrix
            UpwindForNavierStokes(Coefficients[0], SQMATRICES[0], u1, u2);
            cout << "UPWINDING DONE : level " << endl;
            break;

          case 3:
          case 4:
            // do upwinding with two matrices
            cout << "UPWINDING DONE : level " << endl;
            UpwindForNavierStokes(Coefficients[0], SQMATRICES[0], u1, u2);
            UpwindForNavierStokes(Coefficients[0],SQMATRICES[3], u1, u2);
	    break;
        }                        // endswitch
      }                          // endif


      if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
      {
        if (TDatabase::ParamDB->NSTYPE <4)
        {
          OutPut("For slip with friction bc NSTYPE 4 is ");
          OutPut("necessary !!!!! " << endl);
          exit(4711);
        }
	OutPut("sl"<<endl);
        // prepare everything for the assembling of slip with friction bc
        // on all levels
        N_FESpaces = 1;
        N_SquareMatrices = 4;
        N_RectMatrices = 2;
        N_Rhs = 2;
        DiscreteForm0 = NULL;

        SQMATRICES[0] = SqmatrixA11;
        SQMATRICES[1] = SqmatrixA12;
        SQMATRICES[2] = SqmatrixA21;
        SQMATRICES[3] = SqmatrixA22;

        MATRICES[0] = MatrixB1T;
        MATRICES[1] = MatrixB2T;

        fesp[0] = velocity_space;
        ferhs[0] = velocity_space;
        ferhs[1] = velocity_space;

        RHSs[0] = rhs;
        RHSs[1] = rhs + N_U;

        Assemble2DSlipBC(N_FESpaces, fesp,
          N_SquareMatrices, SQMATRICES,
          N_RectMatrices, MATRICES,
          N_Rhs, RHSs, ferhs,
          DiscreteForm0,
          BoundaryConditions,
          BoundValues,
          aux,
          u1, u2);

        // set MATRICES for solver
        SQMATRICES[0] = SqmatrixA11;
        SQMATRICES[1] = SqmatrixA12;
        SQMATRICES[2] = SqmatrixA21;
        SQMATRICES[3] = SqmatrixA22;
        MATRICES[0] = MatrixB1;
        MATRICES[1] = MatrixB2;
        MATRICES[2] = MatrixB1T;
        MATRICES[3] = MatrixB2T;
       }
      delete aux;

     // set Dirichlet nodes
     memcpy(sol+N_Active, rhs+N_Active, N_NonActive*SizeOfDouble);
     memcpy(sol+N_U+N_Active, rhs+N_U+N_Active, N_NonActive*SizeOfDouble);
    } //if(ActiveProcess)

#ifdef _MPI
  // check all dof mapping informations are available
    ParComm[0]->WaitForMakeDofMapping();
    //ParComm[0]->SetOwnDof();

    ParComm[1]->WaitForMakeDofMapping();
    //ParComm[1]->SetOwnDof();
#endif


#ifdef _MPI
     t_par1 = MPI_Wtime();

     if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
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

//   exit(0);

    // compute defect
    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       IntoL20FEFunction(sol+2*N_U, N_P, pressure_space,
                         velocity_space_code, pressure_space_code
                         #ifdef _MPI
                         , MPI_COMM_WORLD
                         #endif
                        );

#ifdef _MPI
    if(ActiveProcess)
#endif
     {
      defect = new double[N_Unknowns];
      memset(defect,0,N_Unknowns*SizeOfDouble);
     }  //    if(ActiveProcess)

#ifdef _MPI
    if(ActiveProcess)
     ParDefect(sqmatrices,matrices, ParSolVect, rhs, defect);
#else
     Defect(sqmatrices,matrices, sol, rhs, defect);
#endif

    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
      IntoL20Vector2D(defect+2*N_U, N_P, pressure_space_code
#ifdef _MPI
                      , ParComm[1]
#endif
                     );
//         exit(0);
#ifdef _MPI
    ParDdotNSE2D(defect, N_U, N_P, ParComm[0], ParComm[1], residual, impuls_residual);
#else
    residual =  Ddot(N_Unknowns,defect,defect);
    impuls_residual = Ddot(2*N_U,defect,defect);
#endif

    #ifdef _MPI
    if(rank==out_rank)
    #endif
     {
      OutPut("nonlinear iteration step   0");
      OutPut(setw(14) << impuls_residual);
      OutPut(setw(14) << residual-impuls_residual);
      OutPut(setw(14) << sqrt(residual) << endl);
     }
//    MPI_Finalize();
//    exit(0); 
   // initilalize the solver
     // solve the system with a parallel solver    
     switch(int(TDatabase::ParamDB->SOLVER_TYPE))
      {
       #ifndef _MPI
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
             DirectSolver(SqmatrixA, MatrixB1T, MatrixB2T, MatrixB1, MatrixB2, rhs, sol);    
       break;

      #ifdef _OMP
       case 100:
             // pardiso
             cout << "solver type not implemented yet !!!!!!!!!!!!!" << endl;
             exit(0);

       break;
      #endif
      #endif
      #ifdef _MPI
       case 101:
         t1 = MPI_Wtime();
             // MUMPS Parallel solver
         MUMPS_NSESolver = new TMumpsSolver(MPI_COMM_WORLD, sqstructureA, structureB, structureBT, ParComm[0], ParComm[1]);

         if(rank==out_rank)
           OutPut(" Time taken for MUMPS Analysis : " << MPI_Wtime()- t1 << "s" << endl);

         t1 = MPI_Wtime();
         switch(TDatabase::ParamDB->NSTYPE)
          {
            case 2:
               MUMPS_NSESolver->FactorizeAndSolve(SqmatrixA, MatrixB1T, MatrixB2T,  MatrixB1, MatrixB2,
                                                  ParRhsNSEVect,  ParComm[0], ParComm[1], ParSolVect);
               break;
           default:
              OutPut("Direct solver for NSTYPE " << TDatabase::ParamDB->NSTYPE
                      << " not implemented !!!" << endl);
              exit(4711);
          }
         t2 = MPI_Wtime();
         solver_time += (t2-t1);
       break;

       case 102:
   
           printf("Parallel SuperLU solver not yet implemented \n" );
           MPI_Finalize();
           exit(0);

       break;
      #endif
       default: 
        cout << "wrong parallel solver type !!!!!!!!!!!!!" << endl;
        exit(0);
       break;
      }

    // ************************************************************* //
    // end of first nonlinear step
    // ************************************************************* //
    #ifdef _MPI
    if(rank==out_rank)
    #endif
     {
      OutPut(" Time taken by SOLVER : " << solver_time << "s" << endl);
      OutPut("MEMORY: " << setw(10) << GetMemory()/(1048576.0));
      OutPut(" MB" << endl);
     }

    // ************************************************************* //
    // the nonlinear iteration
    // ************************************************************* //
    for(j=1;j<=Max_It;j++)
     {
      #ifdef _MPI
      if(ActiveProcess)
       #endif 
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


      switch(TDatabase::ParamDB->NSTYPE)
       {
        case 1:
          SQMATRICES[0] = SqmatrixA;
          SQMATRICES[0]->Reset();

          N_SquareMatrices = 1;
          N_RectMatrices = 0;

          N_Rhs = 0;
          N_FESpaces = 1;
        break;

        case 2:

          SQMATRICES[0] = SqmatrixA;
          SQMATRICES[0]->Reset();

          N_SquareMatrices = 1;
          N_RectMatrices = 0;
          N_Rhs = 0;
          N_FESpaces = 1;

        break;

        case 3:
          if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
           {
            SQMATRICES[0] = SqmatrixA11;
            SQMATRICES[1] = SqmatrixA22;
            SQMATRICES[0]->Reset();
            SQMATRICES[1]->Reset();

            N_SquareMatrices = 2;
            N_RectMatrices = 0;
            last_sq = 1;

            N_Rhs = 0;
            N_FESpaces = 1;
           }
          else
           {
                      // Newton method
            cout<< "Newton method not tested " <<endl;
            exit(0);
           }


         break;

         case 4:
          if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
           {
            N_SquareMatrices = 2;
            SQMATRICES[0] = SqmatrixA11;
            SQMATRICES[1] = SqmatrixA22;
            SQMATRICES[0]->Reset();
            SQMATRICES[1]->Reset();

            N_RectMatrices = 0;
            N_Rhs = 0;
            N_FESpaces = 1;
            last_sq = 1;
           }
          else
           {                      // Newton method
            cout<< "Newton method not tested " <<endl;
            exit(0);
           }
          break;
        } // switch(TDatabase::ParamDB->NSTYPE)


      aux =  new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo, NSN_ParamFctVelo,
                             NSN_FEValuesVelo,
                             fesp, fefct,
                             NSFctVelo,
                             NSFEFctIndexVelo, NSFEMultiIndexVelo,
                             NSN_ParamsVelo, NSBeginParamVelo);


      Assemble2D(N_FESpaces, fesp,
                 N_SquareMatrices, SQMATRICES,
                 N_RectMatrices, MATRICES,
                 N_Rhs, RHSs, ferhs,
                 DiscreteForm,
                 BoundaryConditions,
                 BoundValues,
                 aux);

    if ((DiscreteForm == DiscreteFormNLUpwind)
      &&(!TDatabase::ParamDB->STOKES_PROBLEM))
     {
      switch(TDatabase::ParamDB->NSTYPE)
       {
        case 1:
        case 2:
          // do upwinding with one matrix
          UpwindForNavierStokes(Coefficients[0], SQMATRICES[0], u1, u2);
          cout << "UPWINDING DONE :  "  << endl;
          break;

            case 3:
            case 4:
              // do upwinding with two matrices
              UpwindForNavierStokes(Coefficients[0], SQMATRICES[0], u1, u2);
              UpwindForNavierStokes(Coefficients[0], SQMATRICES[last_sq], u1, u2);
              cout << "UPWINDING DONE(2) :  " << endl;
              break;
          }                      // endswitch
        }                        // endif

      if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
      {
        // prepare everything for the assembling of slip with friction bc
        // on all levels
        N_FESpaces = 1;
        N_SquareMatrices = 4;
        N_RectMatrices = 0;
        N_Rhs = 0;
        DiscreteForm0 = NULL;

        SQMATRICES[0] = SqmatrixA11;
        SQMATRICES[1] = SqmatrixA12;
        SQMATRICES[2] = SqmatrixA21;
        SQMATRICES[3] = SqmatrixA22;

        fesp[0] = velocity_space;
        ferhs[0] = velocity_space;
        ferhs[1] = velocity_space;

        RHSs[0] = rhs;
        RHSs[1] = rhs + N_U;

        Assemble2DSlipBC(N_FESpaces, fesp,
          N_SquareMatrices, SQMATRICES,
          N_RectMatrices, MATRICES,
          N_Rhs, RHSs, ferhs,
          DiscreteForm0,
          BoundaryConditions,
          BoundValues,
          aux,
          u1, u2);
       }
      delete aux;

      // end of assembling

      // set Dirichlet nodes
      memcpy(sol+N_Active, rhs+N_Active, N_NonActive*SizeOfDouble);
      memcpy(sol+N_U+N_Active, rhs+N_U+N_Active, N_NonActive*SizeOfDouble);
     } //if(ActiveProcess)


      #ifdef _MPI
      if(ActiveProcess)
      #endif
       {
        memset(defect,0,N_Unknowns*SizeOfDouble);
       }  //    if(ActiveProcess)

      #ifdef _MPI
      if(ActiveProcess)
       ParDefect(sqmatrices,matrices, ParSolVect, rhs, defect);
      #else
       Defect(sqmatrices,matrices, sol, rhs, defect);
      #endif

      if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
        IntoL20Vector2D(defect+2*N_U, N_P, pressure_space_code
                        #ifdef _MPI
                        , ParComm[1]
                        #endif
                       );
      #ifdef _MPI
      ParDdotNSE2D(defect, N_U, N_P, ParComm[0], ParComm[1], residual, impuls_residual);
      #else
      residual =  Ddot(N_Unknowns,defect,defect);
      impuls_residual = Ddot(2*N_U,defect,defect);
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

      if ((sqrt(residual)<=limit)||(j==Max_It))
      {
        #ifdef _MPI
        total_time2 = MPI_Wtime();
        #else
        total_time2 = GetTime();
        #endif

        #ifdef _MPI
        if(rank==out_rank)
        #endif
         {
          OutPut("ITE : " << setw(3) << j);
//          OutPut(" (" << setw(3) << N_LinIter << " LINITE)");
          OutPut(" T/SOLVER : " << solver_time << "s");
          OutPut(" T/TOTAL : " << total_time2-total_time1 << "s");
          OutPut("("<<setprecision(2)<<solver_time/(total_time2-total_time1)<<")");
          OutPut(setprecision(6)<<" RES : " <<  sqrt(residual));
          OutPut(endl);
          }

          break;
       }//if ((sqrt(residual)<=limit)||(j==Max_It))

    // solve the system with a parallel solver    
     switch(int(TDatabase::ParamDB->SOLVER_TYPE))
      {
       #ifndef _MPI
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
             DirectSolver(SqmatrixA, MatrixB1T, MatrixB2T, MatrixB1, MatrixB2, rhs, sol);
       break;

      #ifdef _OMP
       case 100:
             // pardiso
             cout << "solver type not implemented yet !!!!!!!!!!!!!" << endl;
             exit(0);

       break;
      #endif
      #endif
      #ifdef _MPI
       case 101:
         t1 = MPI_Wtime();
         switch(TDatabase::ParamDB->NSTYPE)
          {
            case 2:
//                MUMPS_NSESolver->Solve_NSENL(SqmatrixA, ParRhsNSEVect,  ParComm[0], ParSolVect);
               MUMPS_NSESolver->FactorizeAndSolve_NSENL(SqmatrixA, ParRhsNSEVect,  ParComm[0], ParSolVect);
//                MUMPS_NSESolver->FactorizeAndSolve(SqmatrixA, MatrixB1T, MatrixB2T,  MatrixB1, MatrixB2,
//                                                   ParRhsNSEVect,  ParComm[0], ParComm[1], ParSolVect);
               break;
           default:
              OutPut("Direct solver for NSTYPE " << TDatabase::ParamDB->NSTYPE
                      << " not implemented !!!" << endl);
              exit(4711);
          }
         t2 = MPI_Wtime();
         solver_time += (t2-t1);
       break;

       case 102:
           printf("Parallel SuperLU solver not yet implemented \n" );
           MPI_Finalize();
           exit(0);

       break;
      #endif
       default: 
        cout << "wrong parallel solver type !!!!!!!!!!!!!" << endl;
        exit(0);
       break;
      }
     } //  for(j=1;j<=Max_It;j++)

    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       IntoL20FEFunction(sol+2*N_U, N_P, pressure_space,
                         velocity_space_code, pressure_space_code
                         #ifdef _MPI
                         , MPI_COMM_WORLD
                         #endif
                        );


//      #ifdef _MPI
//     if(rank==0)
//      #endif
     u->WriteSol(0);

     #ifdef _MPI
     if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
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

   #ifdef _MPI
   if(rank==out_rank)
    printf("Main root used bytes %f\n", GetMemory());
   #endif

  CloseFiles();
  #ifdef _MPI
  MPI_Finalize();
  #endif
  return 0;
}
