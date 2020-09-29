#include <MooNMD.h>

#include <stdlib.h>
#include <stdio.h>

// #include "../Examples/TNSE_3D/TwoPhase_NavierStokes.h"
// #include "../Examples/TNSE_3D/Benchmark3D.h"
#include "../Examples/TNSE_3D/AnsatzLinConst.h"

#define EXITMAIN(x) MPI_Barrier(comm); MPI_Finalize();  return (x)

#define ROOTOUT(x) if ( rank == 0 ) { OutPut(x); } 

int g_rank = 0;

typedef struct _MUTABLE_DATA
{
  TCollection *Coll;
  TCollection *Coll_P0;
  TCollection *Coll_P1;
  
  int *GlobalCellNo_P0;
  int *GlobalCellNo_P1;
  
  int N_SurfaceJoints_P0;
  int N_SurfaceJoints_P1;
  
  int *CellNumbers_P0;
  int *CellNumbers_P1;
  int *JointNumbers_P0;
  int *JointNumbers_P1;
  
  TFESpace3D *velocity_space;
  TFESpace3D *pressure_space;
  TFESpace3D *grid_space;
  TFESpace3D *grid_space_P0;
  TFESpace3D *grid_space_P1;
  int pressure_space_code;
  
  TParFECommunicator3D *velocity_comm;
  TParFECommunicator3D *pressure_comm;
  TParFECommunicator3D *grid_comm_p0, *grid_comm_p1;
  
  int N_U, N_P, N_Grid, N_Grid_P0, N_Grid_P1;
  int N_Unknowns, ActiveBound;
  
  double *sol, *rhs, *res;
  double *grid_velo, *grid_pos;
  double *grid_rhs_p0, *grid_rhs_p1;
  double *grid_sol_p0, *grid_sol_p1;
  double *aux_array;
  
  TParVector3D *Sol_ux, *Sol_uy, *Sol_uz, *Sol_p;
  TParVector3D *Rhs_ux, *Rhs_uy, *Rhs_uz, *Rhs_p;
  TParVector3D *Res_ux, *Res_uy, *Res_uz, *Res_p;
  TParVector3D *GridRhs_P0_x, *GridRhs_P0_y, *GridRhs_P0_z;
  TParVector3D *GridSol_P0_x, *GridSol_P0_y, *GridSol_P0_z;
  TParVector3D *GridRhs_P1_x, *GridRhs_P1_y, *GridRhs_P1_z;
  TParVector3D *GridSol_P1_x, *GridSol_P1_y, *GridSol_P1_z;
  
  TFEVectFunct3D *u;
  TFEFunction3D *ux, *uy, *uz;
  TFEVectFunct3D *g;
  TFEFunction3D *gx, *gy, *gz;
  TFEFunction3D *p;
  TFEVectFunct3D *g_p0, *g_p1;
  TFEVectFunct3D *grid;
  
  TFEFunction3D *FEFCT[6];
  
  TFEFunction3D *phase_marker;
  
  TSquareStructure3D *sqstructureA, *sqstructureGrid_P0, *sqstructureGrid_P1;
  TStructure3D *structureB, *structureBT;
  
  TSquareMatrix3D *A11, *A12, *A13;
  TSquareMatrix3D *A21, *A22, *A23;
  TSquareMatrix3D *A31, *A32, *A33;
  TSquareMatrix3D *M11, *M22, *M33;
  
  TSquareMatrix3D *G11_P0, *G12_P0, *G13_P0;
  TSquareMatrix3D *G21_P0, *G22_P0, *G23_P0;
  TSquareMatrix3D *G31_P0, *G32_P0, *G33_P0;
  
  TSquareMatrix3D *G11_P1, *G12_P1, *G13_P1;
  TSquareMatrix3D *G21_P1, *G22_P1, *G23_P1;
  TSquareMatrix3D *G31_P1, *G32_P1, *G33_P1;
  
  TMatrix3D *B1, *B2, *B3;
  TMatrix3D *B1T, *B2T, *B3T;
  
  TAuxParam3D *aux;
  
  TOutput3D *Output;
  
  TMumpsSolver *Solver;
  TMumpsSolver *Solver_Grid_P0;
  TMumpsSolver *Solver_Grid_P1;
  
  void Create(TDomain *Domain, MPI_Comm comm);
  
} MUTABLE_DATA;

void WriteGrid(const char *basename, int number, TDomain *Domain, TCollection *Phase);
void WriteGrid(const char *basename, TDomain *Domain, TCollection *Phase);
void WriteSolution(const char *basename, int number, TOutput3D *Output);

void PrintMemory();
void PrintMemory(int bmem);

void Remesh(MPI_Comm comm, TDomain *Domain, MUTABLE_DATA *Data);
void GetFacets(MUTABLE_DATA *Data, int &N_Points, double *&Points, int *&Facets);
double CheckForRemesh(TCollection *Coll);
void FreeDomain(TTetGenMeshLoader::TDummyDomain *Domain);

int main (int argc, char **argv)
{
  // init mpi
  int size, rank, status;
  MPI_Comm comm = MPI_COMM_WORLD;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  g_rank = rank;
  
//   #ifdef __DEBUG
//   MPI_Barrier(comm);
//   if ( rank == 0 )
//   {
//     int i = 0;
//     char hostname[256];
//     gethostname(hostname, sizeof(hostname));
//     OutPut("PID "<<getpid()<<" on "<<hostname<<" ready for attach"<<endl);
//     fflush(stdout);
//     while (0 == i)
//         sleep(5);
//   }  
//   MPI_Barrier(comm);
//   #endif
  
  TDatabase Database;
  TFEDatabase3D FEDatabase;
  TDomain Domain;
  
  MUTABLE_DATA Data;
  
  TDiscreteForm3D *DiscreteFormGalerkin, *DiscreteFormNLGalerkin;
  TDiscreteForm3D *DiscreteFormRhs, *DiscreteFormGrid;
  
  int N_Remesh = 0, MaxCpV;
  
  char Name[] = "name";
  char Desc[] = "desc";
  char Inter[] = "Isosurface";
  
  double Residual_norm, mesh_quality, mesh_quality_loc;
  double qualitybound;
  
  if ( argc > 1 )
    Domain.ReadParam(argv[1]);
  else
  {
    cerr << "no readin file given" << endl;
    EXITMAIN(0);
  }
  
  OpenFiles(comm);
  
  char *SMESH = TDatabase::ParamDB->SMESHFILE;
  char *VTK = TDatabase::ParamDB->VTKBASENAME;
  int NSTYPE = TDatabase::ParamDB->NSTYPE;
  
  qualitybound = TDatabase::ParamDB->P7*TDatabase::ParamDB->TETGEN_QUALITY;
  
  if ( NSTYPE != 4 )
  {
    cerr << "Use NSTYPE = 4 !" << endl;
    EXITMAIN(0);
  }
  
  //
  BoundCondFunct3D      *BOUNDCOND[4], *GRIDBOUNDCOND[4];
  BoundValueFunct3D     *BOUNDVALUE[4], *GRIDBOUNDVALUE[4];
  
  BOUNDCOND[0] = BoundCondition;
  BOUNDCOND[1] = BoundCondition;
  BOUNDCOND[2] = BoundCondition;
  
  BOUNDVALUE[0] = U1BoundValue;
  BOUNDVALUE[1] = U2BoundValue;
  BOUNDVALUE[2] = U3BoundValue;
  
  GRIDBOUNDCOND[0] = GridBoundCond;
  GRIDBOUNDCOND[1] = GridBoundCond;
  GRIDBOUNDCOND[2] = GridBoundCond;
  
  GRIDBOUNDVALUE[0] = GridXBoundValues;
  GRIDBOUNDVALUE[1] = GridYBoundValues;
  GRIDBOUNDVALUE[2] = GridZBoundValues;
  
  ExampleFile();
  
  // discreteforms
  InitializeDiscreteForms(DiscreteFormGalerkin, DiscreteFormNLGalerkin,
			  DiscreteFormRhs, LinCoeffs, NSTYPE);
			  
  DiscreteFormGrid = new TDiscreteForm3D (Name, Desc, ESGridN_Terms, ESGridDerivatives,
					  ESGridSpaceNumbers, ESGridN_Matrices,
					  ESGridN_Rhs, ESGridRowSpace, ESGridColumnSpace,
					  ESGridRhsSpace, ESGridAssemble, GridCoeffs, NULL);
					  
  // init grid
  Domain.Tetgen(SMESH);
  
  for (int i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;++i)
    Domain.RegRefineAll();
  
  Partition_Mesh(comm, &Domain, MaxCpV);
// EXITMAIN(0);
  Data.Create(&Domain, comm);
// EXITMAIN(0);
  TAuxParam3D aux_dummy (0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
  
  TSquareMatrix3D *SQMATRICES[12];
  TMatrix3D *MATRICES[6];
  TFESpace3D *FESP[4], *FESPRHS[4];
  int N_FESP, N_SQMAT, N_MAT, N_RHS;
  double *RHS[3];
  double *dptr1, *dptr2;
  
  /// time loop (back euler);
  int img = 0;
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  
  if ( Data.velocity_comm->IsMember() )
  {
    Data.ux->Interpolate(InitialU1);
    Data.uy->Interpolate(InitialU2);
    Data.uz->Interpolate(InitialU3);
    Data.p->Interpolate(InitialP);
  }
  
  if ( rank == 0 ) // root
  {
    // print some informations
    OutPut(endl);
    OutPut("RE_NR: " << TDatabase::ParamDB->RE_NR << endl);
    OutPut("WE_NR: " << TDatabase::ParamDB->WB_NR << endl);
    
    OutPut(endl);
    OutPut("STARTTIME: " << TDatabase::TimeDB->CURRENTTIME << endl);
    OutPut("ENDTIME: " << TDatabase::TimeDB->ENDTIME << endl);
    OutPut("TIMESTEPLENGTH: " << tau << endl);
//     double Volumen = GetVolume(Data.Coll);
//     double Volumen_P0 = GetVolume(Data.Coll_P0);
//     double Volumen_P1 = GetVolume(Data.Coll_P1);
//     OutPut("Volume: " << Volumen << endl);
//     OutPut("Volume phase 0: " << Volumen_P0 << endl);
//     OutPut("Volume phase 1: " << Volumen_P1 << endl);
    OutPut("Density ratio: " << TDatabase::ParamDB->P6 << endl);
    OutPut("Viscosity ratio: " << TDatabase::ParamDB->P5 << endl << endl);
  }
  
  WriteSolution(VTK, img++, Data.Output);
  
  PrintMemory();
// EXITMAIN(0);
//   MPI_Barrier(comm);
  /// time loop begin
  for(int timestep=0;TDatabase::TimeDB->CURRENTTIME <= TDatabase::TimeDB->ENDTIME;++timestep)
  {
    MPI_Barrier(comm);
    
    TDatabase::TimeDB->CURRENTTIME += tau;
    
    ROOTOUT (endl << "======================================================="<<endl);
    ROOTOUT (        "CURRENTTIME: " << TDatabase::TimeDB->CURRENTTIME << endl);
    ROOTOUT (        "======================================================="<<endl);
    
    if ( Data.velocity_comm->IsMember() )
      memset(Data.rhs, 0, Data.N_Unknowns*sizeof(Data.rhs[0]));
;
    /// assemble systems (Navier-Stokes)
    if ( Data.velocity_comm->IsWorker() )
    {
      FESP[0] = Data.velocity_space;
      FESP[1] = Data.pressure_space;
      FESP[2] = Data.grid_space;
      N_FESP = 3;
      
      FESPRHS[0] = Data.velocity_space;
      FESPRHS[1] = Data.velocity_space;
      FESPRHS[2] = Data.velocity_space;
      
      RHS[0] = Data.rhs        ;
      RHS[1] = Data.rhs +   Data.N_U;
      RHS[2] = Data.rhs + 2*Data.N_U;
      N_RHS = 3;
      
      SQMATRICES[0] = Data.A11;
      SQMATRICES[1] = Data.A12;
      SQMATRICES[2] = Data.A13;
      SQMATRICES[3] = Data.A21;
      SQMATRICES[4] = Data.A22;
      SQMATRICES[5] = Data.A23;
      SQMATRICES[6] = Data.A31;
      SQMATRICES[7] = Data.A32;
      SQMATRICES[8] = Data.A33;
      SQMATRICES[9] = Data.M11;
      SQMATRICES[10] = Data.M22;
      SQMATRICES[11] = Data.M33;
      N_SQMAT = 12;
      
      MATRICES[0] = Data.B1;
      MATRICES[1] = Data.B2;
      MATRICES[2] = Data.B3;
      MATRICES[3] = Data.B1T;
      MATRICES[4] = Data.B2T;
      MATRICES[5] = Data.B3T;
      N_MAT = 6;
      
      MatReset(N_SQMAT, SQMATRICES, N_MAT, MATRICES);
      
      // assemble
      Assemble3D (N_FESP, FESP,
		  N_SQMAT, SQMATRICES,
		  N_MAT, MATRICES,
		  N_RHS, RHS, FESPRHS,
		  DiscreteFormGalerkin,
		  BOUNDCOND, BOUNDVALUE, Data.aux);
		       
      OutPut("Assemble done" << endl);
      
      // scale BT's
      Dscal(Data.B1T->GetN_Entries(), tau, Data.B1T->GetEntries());
      Dscal(Data.B2T->GetN_Entries(), tau, Data.B2T->GetEntries());
      Dscal(Data.B3T->GetN_Entries(), tau, Data.B3T->GetEntries());
      
      // add M*u_old to rhs and scale rhs
      SQMATRICES[0] = Data.M11;
      SQMATRICES[1] = Data.M22;
      SQMATRICES[2] = Data.M33;
      for (int i=0;i<3;++i)
      {
	dptr1 = Data.sol + i*Data.N_U;
	dptr2 = Data.rhs + i*Data.N_U;
	
	MatVectActive(SQMATRICES[i], dptr1, Data.aux_array);
	Dscal(Data.ActiveBound, tau, dptr2);
	Daxpy(Data.ActiveBound, 1.0, Data.aux_array, dptr2);
      }
      
      // add surface integrals to M
      FreeSurfInt_TwoPhase(Data.N_SurfaceJoints_P1, Data.CellNumbers_P1,
			  Data.JointNumbers_P1, Data.GlobalCellNo_P1, tau, SQMATRICES,
			  RHS[0], RHS[1], RHS[2]);
			  
      // add M block to A blocks and scale A blocks
      MatAdd2(Data.A11, Data.M11, tau);
      MatAdd2(Data.A22, Data.M22, tau);
      MatAdd2(Data.A33, Data.M33, tau);
      
      SQMATRICES[0] = Data.A12;
      SQMATRICES[1] = Data.A13;
      SQMATRICES[2] = Data.A21;
      SQMATRICES[3] = Data.A23;
      SQMATRICES[4] = Data.A31;
      SQMATRICES[5] = Data.A32;
      
      for (int i=0;i<6;++i)
	Dscal(SQMATRICES[i]->GetN_Entries(), tau, SQMATRICES[i]->GetEntries());
     
    } // end assemble systems (Navier-Stokes)
    
    /// assemble grid equation
    // grid phase 0;
    if ( Data.grid_comm_p0->IsWorker() )
    {
      FESP[0] = Data.grid_space_P0;
      N_FESP = 1;
		  
      SQMATRICES[0] = Data.G11_P0;
      SQMATRICES[1] = Data.G12_P0;
      SQMATRICES[2] = Data.G13_P0;
      SQMATRICES[3] = Data.G21_P0;
      SQMATRICES[4] = Data.G22_P0;
      SQMATRICES[5] = Data.G23_P0;
      SQMATRICES[6] = Data.G31_P0;
      SQMATRICES[7] = Data.G32_P0;
      SQMATRICES[8] = Data.G33_P0;
      N_SQMAT = 9;
      N_MAT = 0;
      
      MatReset(N_SQMAT, SQMATRICES, N_MAT, NULL);
      
      Assemble3D (N_FESP, FESP,
		  N_SQMAT, SQMATRICES,
		  N_MAT, NULL,
		  0, NULL, NULL,
		  DiscreteFormGrid,
		  GRIDBOUNDCOND, GRIDBOUNDVALUE, &aux_dummy);
    }	  
    // grid phase 1
    if ( Data.grid_comm_p1->IsWorker() )
    {
      FESP[0] = Data.grid_space_P1;
      N_FESP = 1;
		  
      SQMATRICES[0] = Data.G11_P1;
      SQMATRICES[1] = Data.G12_P1;
      SQMATRICES[2] = Data.G13_P1;
      SQMATRICES[3] = Data.G21_P1;
      SQMATRICES[4] = Data.G22_P1;
      SQMATRICES[5] = Data.G23_P1;
      SQMATRICES[6] = Data.G31_P1;
      SQMATRICES[7] = Data.G32_P1;
      SQMATRICES[8] = Data.G33_P1;
      N_SQMAT = 9;
      N_MAT = 0;
      
      MatReset(N_SQMAT, SQMATRICES, N_MAT, NULL);
      
      Assemble3D (N_FESP, FESP,
		  N_SQMAT, SQMATRICES,
		  N_MAT, NULL,
		  0, NULL, NULL,
		  DiscreteFormGrid,
		  GRIDBOUNDCOND, GRIDBOUNDVALUE, &aux_dummy);
    }

    MPI_Barrier(comm);
    
    // collect rhs on root and solve linear system
    if ( Data.velocity_comm->IsMember() )
    {		
      Data.Rhs_ux->AssembleAtRoot();
      Data.Rhs_uy->AssembleAtRoot();
      Data.Rhs_uz->AssembleAtRoot();
      Data.Rhs_p->AssembleAtRoot();
      
      /// Solve sytem
      Data.Solver->SetMatrix(Data.A11, Data.A12, Data.A13, Data.A21, Data.A22, Data.A23,
			      Data.A31, Data.A32, Data.A33, Data.B1T, Data.B2T, Data.B3T,
			      Data.B1, Data.B2, Data.B3);
	  
      Data.Solver->Factorize();
      Data.Solver->Solve(Data.sol, Data.rhs);
      
      // distribute solution to workers
      Data.Sol_ux->ScatterFromRoot();
      Data.Sol_uy->ScatterFromRoot();
      Data.Sol_uz->ScatterFromRoot();
      Data.Sol_p->ScatterFromRoot();
    }
   
    ///===================================================================
    /// START nonlinear loop =============================================
    ///===================================================================
    for (int iter=0;iter<TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;++iter)
    {
      /// solve grid equation
      // phase 0
      if ( Data.grid_comm_p0->IsMember() )
      {
	memset(Data.grid_rhs_p0, 0, 3*Data.N_Grid_P0*sizeof(Data.grid_rhs_p0[0]));
	
	if ( 0 <= (int) TDatabase::ParamDB->P1 && Data.grid_comm_p0->IsMember() )
	{
	  SetGridRhs(Data.u, Data.grid_space_P0, Data.N_SurfaceJoints_P0, Data.CellNumbers_P0,
		  Data.JointNumbers_P0, Data.GlobalCellNo_P0, tau, Data.grid_rhs_p0, 
		  Data.grid_rhs_p0+Data.N_Grid_P0, Data.grid_rhs_p0+2*Data.N_Grid_P0);
	}
      }	
	// collect rhs on root and solve
      if ( 0 < (int) TDatabase::ParamDB->P1 && Data.grid_comm_p0->IsMember() )
      {
	Data.GridRhs_P0_x->AssembleAtRoot();
	Data.GridRhs_P0_y->AssembleAtRoot();
	Data.GridRhs_P0_z->AssembleAtRoot();
		  
	SQMATRICES[0] = Data.G11_P0;
	SQMATRICES[1] = Data.G12_P0;
	SQMATRICES[2] = Data.G13_P0;
	SQMATRICES[3] = Data.G21_P0;
	SQMATRICES[4] = Data.G22_P0;
	SQMATRICES[5] = Data.G23_P0;
	SQMATRICES[6] = Data.G31_P0;
	SQMATRICES[7] = Data.G32_P0;
	SQMATRICES[8] = Data.G33_P0;
	  
	Data.Solver_Grid_P0->SetMatrix(SQMATRICES);
	Data.Solver_Grid_P0->Factorize();
	Data.Solver_Grid_P0->Solve(Data.grid_sol_p0, Data.grid_rhs_p0);
	
	// distribute solution to worker
	Data.GridSol_P0_x->ScatterFromRoot();
	Data.GridSol_P0_y->ScatterFromRoot();
	Data.GridSol_P0_z->ScatterFromRoot();
      }
      else if ( Data.grid_comm_p0->IsMember() )
      {
// 	Data.GridRhs_P0_x->AssembleAtRoot();
// 	Data.GridRhs_P0_y->AssembleAtRoot();
// 	Data.GridRhs_P0_z->AssembleAtRoot();
// 	
// 	Data.GridSol_P0_x->ScatterFromRoot();
// 	Data.GridSol_P0_y->ScatterFromRoot();
// 	Data.GridSol_P0_z->ScatterFromRoot();
	
	memcpy(Data.grid_sol_p0, Data.grid_rhs_p0, 3*Data.N_Grid_P0*sizeof(Data.grid_sol_p0[0]));
      }
            
      // phase 1
      if ( Data.grid_comm_p1->IsMember() )
      {
	memset(Data.grid_rhs_p1, 0, 3*Data.N_Grid_P1*sizeof(Data.grid_rhs_p1[0]));
      
	if ( 0 <= (int) TDatabase::ParamDB->P1 && Data.grid_comm_p1->IsMember() )
	{		  
	  SetGridRhs(Data.u, Data.grid_space_P1, Data.N_SurfaceJoints_P1, Data.CellNumbers_P1,
		  Data.JointNumbers_P1, Data.GlobalCellNo_P1, tau, Data.grid_rhs_p1,
		  Data.grid_rhs_p1+Data.N_Grid_P1, Data.grid_rhs_p1+2*Data.N_Grid_P1);
	}
      }
      
      if ( 0 < (int) TDatabase::ParamDB->P1 && Data.grid_comm_p1->IsMember() )
      {
	Data.GridRhs_P1_x->AssembleAtRoot();
	Data.GridRhs_P1_y->AssembleAtRoot();
	Data.GridRhs_P1_z->AssembleAtRoot();
	
	SQMATRICES[0] = Data.G11_P1;
	SQMATRICES[1] = Data.G12_P1;
	SQMATRICES[2] = Data.G13_P1;
	SQMATRICES[3] = Data.G21_P1;
	SQMATRICES[4] = Data.G22_P1;
	SQMATRICES[5] = Data.G23_P1;
	SQMATRICES[6] = Data.G31_P1;
	SQMATRICES[7] = Data.G32_P1;
	SQMATRICES[8] = Data.G33_P1;
	
	Data.Solver_Grid_P1->SetMatrix(SQMATRICES);
	Data.Solver_Grid_P1->Factorize();
	Data.Solver_Grid_P1->Solve(Data.grid_sol_p1, Data.grid_rhs_p1);
	
	// distribute solution to worker
	Data.GridSol_P1_x->ScatterFromRoot();
	Data.GridSol_P1_y->ScatterFromRoot();
	Data.GridSol_P1_z->ScatterFromRoot();
      }
      else if ( Data.grid_comm_p1->IsMember() )
      {
// 	Data.GridRhs_P1_x->AssembleAtRoot();
// 	Data.GridRhs_P1_y->AssembleAtRoot();
// 	Data.GridRhs_P1_z->AssembleAtRoot();
// 	
// 	Data.GridSol_P1_x->ScatterFromRoot();
// 	Data.GridSol_P1_y->ScatterFromRoot();
// 	Data.GridSol_P1_z->ScatterFromRoot();
	
	memcpy(Data.grid_sol_p1, Data.grid_rhs_p1, 3*Data.N_Grid_P1*sizeof(Data.grid_sol_p1[0]));
      }
      
      if ( Data.grid_comm_p0->IsMember() )
	MapGridVelo(Data.grid_space_P0, Data.g, Data.GlobalCellNo_P0, Data.grid_sol_p0);
      if ( Data.grid_comm_p1->IsMember() )
	MapGridVelo(Data.grid_space_P1, Data.g, Data.GlobalCellNo_P1, Data.grid_sol_p1);	
      
      if ( Data.velocity_comm->IsMember() )
	Dscal(3*Data.N_Grid, 1/tau, Data.grid_velo);
      
      MPI_Barrier(comm);

      /// assemble nonlinear part
      if ( Data.velocity_comm->IsWorker() ) // worker
      {
	FESP[0] = Data.velocity_space;
	FESP[1] = Data.grid_space;
	N_FESP = 2;
	
	SQMATRICES[0] = Data.A11;
	SQMATRICES[1] = Data.A22;
	SQMATRICES[2] = Data.A33;
	N_SQMAT = 3;
	N_MAT = 0,
	
	MatReset(N_SQMAT, SQMATRICES, N_MAT, NULL);
	
	Assemble3D (N_FESP, FESP,
		    N_SQMAT, SQMATRICES,
		    0, NULL,
		    0, NULL, NULL,
		    DiscreteFormNLGalerkin,
		    BOUNDCOND, BOUNDVALUE, Data.aux);
		    
	// add M blocks to A blocks 
	MatAdd2(Data.A11, Data.M11, tau);
	MatAdd2(Data.A22, Data.M22, tau);
	MatAdd2(Data.A33, Data.M33, tau);
	
	// calculate defect
	CoupledDefect(Data.A11, Data.A12, Data.A13, Data.A21, Data.A22, Data.A23,
		      Data.A31, Data.A32, Data.A33, Data.B1, Data.B2, Data.B3,
		      Data.B1T, Data.B2T, Data.B3T, Data.sol, Data.rhs, Data.res);
      } // end assemble nonlinear part
      
//       MPI_Barrier(comm);
      // collect residual on root
      if ( Data.velocity_comm->IsMember() )
      {
	Data.Res_ux->AssembleAtRoot();
	Data.Res_uy->AssembleAtRoot();
	Data.Res_uz->AssembleAtRoot();
	Data.Res_p->AssembleAtRoot();
      }
      
//       MPI_Barrier(comm);
      if ( Data.velocity_comm->IsRoot() ) // root
      {
	if ( TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE )
	  IntoL20Vector3D(Data.res+3*Data.N_U, Data.N_P, TDatabase::ParamDB->PRESSURE_SPACE);
	
	Residual_norm = Dnorm(Data.N_Unknowns, Data.res);
      }
      
//       MPI_Barrier(comm);
      if ( Data.velocity_comm->IsMember() )
	MPI_Bcast (&Residual_norm, 1, MPI_DOUBLE, 0, Data.velocity_comm->GetMPI_Comm());
      
//       MPI_Barrier(comm);
      ROOTOUT ("Nonlinear iterationstep: " << iter);
      ROOTOUT ("  ; RES: " << Residual_norm << endl);
      
//       MPI_Barrier(comm);
      if (  Residual_norm < TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE )
	break;
      
//       MPI_Barrier(comm);
      if ( Data.velocity_comm->IsMember() )
      {
	Data.Solver->SetMatrix (Data.A11, Data.A12, Data.A13, Data.A21, Data.A22, Data.A23,
				Data.A31, Data.A32, Data.A33, Data.B1T, Data.B2T, Data.B3T,
				Data.B1, Data.B2, Data.B3);
				
	MPI_Barrier(comm);
	Data.Solver->Factorize();
	MPI_Barrier(comm);
	Data.Solver->Solve(Data.sol, Data.rhs);
	
	MPI_Barrier(comm);
	// distribute solution to workers
	Data.Sol_ux->ScatterFromRoot();
	Data.Sol_uy->ScatterFromRoot();
	Data.Sol_uz->ScatterFromRoot();
	Data.Sol_p->ScatterFromRoot();
      }
    }
    ///===================================================================
    /// END nonlinear loop ===============================================
    ///===================================================================
// EXITMAIN(0);
    if ( Data.velocity_comm->IsRoot() ) // root
    {
      OutPut(endl);
      if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
      {
	OutPut("PROJECTING PRESSURE" << endl);
	IntoL20FEFunction3D(Data.p->GetValues(), Data.p->GetLength(), Data.p->GetFESpace3D(),
			    TDatabase::ParamDB->VELOCITY_SPACE, TDatabase::ParamDB->PRESSURE_SPACE);
      }
    }
 
//     MPI_Barrier(comm);
//     if ( velocity_comm->IsMember() )
      Data.Sol_p->ScatterFromRoot();
    
    // Move grid
    Daxpy(3*Data.N_Grid, tau, Data.grid_velo, Data.grid_pos);
    Data.grid->DataToGrid();
    if ( (int) TDatabase::ParamDB->USE_ISOPARAMETRIC && TDatabase::ParamDB->P1 >= 0)
    {
      MoveIsoVertices(Data.u, Data.N_SurfaceJoints_P1, Data.CellNumbers_P1,
		      Data.JointNumbers_P1, Data.GlobalCellNo_P1, tau);
    }
    
    WriteSolution(VTK, img++, Data.Output);
//     if ( rank == 0 ) // root
//     {
//       DumpIsosurface(Inter, img, Data.velocity_space, Data.sol, Data.Coll_P1, Data.N_SurfaceJoints_P1,
// 		 Data.CellNumbers_P1, Data.JointNumbers_P1, Data.GlobalCellNo_P1);
//     }
     
    PrintMemory();
    
    mesh_quality_loc = 0;
    if ( rank != 0 ) // worker
    {
      mesh_quality_loc = CheckForRemesh(Data.Coll);
    }
    
    MPI_Barrier(comm);
    MPI_Reduce(&mesh_quality_loc, &mesh_quality, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    ROOTOUT("mesh quality: " << mesh_quality << endl);
    
    status = 0;
    if ( rank == 0 && mesh_quality > qualitybound )
      status = 1;      
    
    MPI_Bcast(&status, 1, MPI_INT, 0, comm);
    
    if ( status )
    {
      MPI_Barrier(comm);
      Remesh(comm, &Domain, &Data);
    }
  } // end timeloop
  
  CloseFiles();
  MPI_Finalize();
  
  return 0;
}

void WriteGrid(const char *basename, int number, TDomain *Domain, TCollection *Phase)
{
  std::ostringstream os;
  
  os << basename << "." << g_rank << ".";
  os.width(5);
  os.fill('0');
  os << number << ".vtk" << ends;
  
  TOutput3D Output(0,0,0,0, Domain, Phase);
  Output.WriteVtk(os.str().c_str());
}

void WriteGrid(const char *basename, TDomain *Domain, TCollection *Phase)
{
  std::ostringstream os;
  
  os << basename << ".";
  os.width(5);
  os.fill('0');
  os << g_rank << ".vtk" << ends;
  
  TOutput3D Output(0,0,0,0, Domain, Phase);
  Output.WriteVtk(os.str().c_str());
}

void MUTABLE_DATA::Create(TDomain *Domain, MPI_Comm comm)
{
  int rank, size;
//   MPI_Comm comm = *Comm;
  
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  
  char Name[] = "name";
  char Desc[] = "desc";
  char UName[] = "u";
  char PName[] = "p";
  char GName[] = "g";
  
  SplitGrid_TwoPhase(Domain, Coll, Coll_P0, Coll_P1, phase_marker,
		     GlobalCellNo_P0, GlobalCellNo_P1,
		     N_SurfaceJoints_P0, CellNumbers_P0, JointNumbers_P0,
		     N_SurfaceJoints_P1, CellNumbers_P1, JointNumbers_P1);
		    
  Output = new TOutput3D (3, 1, 2, 0, Domain);
  Output->AddFEFunction(phase_marker);
  Output->WriteVtk("phase.vtk");
  
  /// velocity and pressure space
  GetVelocityAndPressureSpace3D (Coll, BoundCondition, velocity_space, pressure_space,
				 &pressure_space_code,
				 TDatabase::ParamDB->VELOCITY_SPACE,
				 TDatabase::ParamDB->PRESSURE_SPACE);
 
  /// grid spaces
  grid_space = new TFESpace3D (Coll, Name, Desc, GridBoundCond, 1);
  grid_space_P0 = new TFESpace3D (Coll_P0, Name, Desc, GridBoundCond, 1);
  grid_space_P1 = new TFESpace3D (Coll_P1, Name, Desc, GridBoundCond, 1);
  
  if ( TDatabase::ParamDB->USE_ISOPARAMETRIC )
  {
    SetInterfaceJoints(N_SurfaceJoints_P0, CellNumbers_P0, JointNumbers_P0,
		     N_SurfaceJoints_P1, CellNumbers_P1, JointNumbers_P1,
		     GlobalCellNo_P0, GlobalCellNo_P1, Coll);
  }
  
  /// fe communicators
  velocity_comm = new TParFECommunicator3D (comm, velocity_space); 
  pressure_comm = new TParFECommunicator3D (comm, pressure_space); 
  grid_comm_p0 = new TParFECommunicator3D (comm, grid_space_P0, GlobalCellNo_P0);
  grid_comm_p1 = new TParFECommunicator3D (comm, grid_space_P1, GlobalCellNo_P1);

  /// count DOF's
  N_U = velocity_space->GetN_DegreesOfFreedom();
  N_P = pressure_space->GetN_DegreesOfFreedom();
  N_Grid_P0 = grid_space_P0->GetN_DegreesOfFreedom();
  N_Grid_P1 = grid_space_P1->GetN_DegreesOfFreedom();
  N_Grid = grid_space->GetN_DegreesOfFreedom();
  N_Unknowns = 3*N_U + N_P;
  ActiveBound = velocity_space->GetActiveBound();
  
  /// alloc arrays
  // arrays for navier stokes
  if ( velocity_comm->IsMember() )
  {
    sol = new double [N_Unknowns];
    rhs = new double [N_Unknowns];
    res = new double [N_Unknowns];
    grid_velo = new double [3*N_Grid];
    grid_pos = new double [3*N_Grid];
    aux_array = new double [N_U];
     
    memset(sol, 0, N_Unknowns*sizeof(sol[0]));
    memset(rhs, 0, N_Unknowns*sizeof(rhs[0]));
    memset(res, 0, N_Unknowns*sizeof(res[0]));
    memset(grid_velo, 0, 3*N_Grid*sizeof(grid_velo[0]));
    memset(grid_pos, 0, 3*N_Grid*sizeof(grid_pos[0]));
    memset(aux_array, 0, N_U*sizeof(aux_array[0]));
  }
  else
  {
    sol = rhs = res = grid_velo = grid_pos = NULL;
  }
  
  // arrays for grid equation phase 0
  if ( grid_comm_p0->IsMember() )
  {
    grid_sol_p0 = new double [3*N_Grid_P0];
    grid_rhs_p0 = new double [3*N_Grid_P0];
    
    memset(grid_rhs_p0, 0, 3*N_Grid_P0*sizeof(grid_rhs_p0[0]));
    memset(grid_sol_p0, 0, 3*N_Grid_P0*sizeof(grid_sol_p0[0]));  
  }
  else
  {
    grid_sol_p0 = grid_rhs_p0 = NULL;
  }
  
  // arrays for grid equation phase 1
  if ( grid_comm_p1->IsMember() )
  {
    grid_sol_p1 = new double [3*N_Grid_P1];
    grid_rhs_p1 = new double [3*N_Grid_P1]; 
    
    memset(grid_rhs_p1, 0, 3*N_Grid_P1*sizeof(grid_rhs_p1[0]));
    memset(grid_sol_p1, 0, 3*N_Grid_P1*sizeof(grid_sol_p1[0]));
  }
  else
  {
    grid_sol_p1 = grid_rhs_p1 = NULL;
  }
    
  /// par vectors;
  // velocity
  if ( velocity_comm->IsMember() )
  {
    Sol_ux = new TParVector3D (velocity_comm, N_U, sol);
    Sol_uy = new TParVector3D (velocity_comm, N_U, sol +   N_U);
    Sol_uz = new TParVector3D (velocity_comm, N_U, sol + 2*N_U);
    Sol_p  = new TParVector3D (pressure_comm, N_P, sol + 3*N_U);
  
    Rhs_ux = new TParVector3D (velocity_comm, N_U, rhs);
    Rhs_uy = new TParVector3D (velocity_comm, N_U, rhs +   N_U);
    Rhs_uz = new TParVector3D (velocity_comm, N_U, rhs + 2*N_U);
    Rhs_p  = new TParVector3D (pressure_comm, N_P, rhs + 3*N_U);
  
    Res_ux = new TParVector3D (velocity_comm, N_U, res);
    Res_uy = new TParVector3D (velocity_comm, N_U, res +   N_U);
    Res_uz = new TParVector3D (velocity_comm, N_U, res + 2*N_U);
    Res_p  = new TParVector3D (pressure_comm, N_P, res + 3*N_U);
  }
  else
  {
    Sol_ux = NULL;
    Sol_uy = NULL;
    Sol_uz = NULL;
    Sol_p  = NULL;  
    Rhs_ux = NULL;
    Rhs_uy = NULL;
    Rhs_uz = NULL;
    Rhs_p  = NULL;  
    Res_ux = NULL;
    Res_uy = NULL;
    Res_uz = NULL;
    Res_p  = NULL;
  }
  
  if ( grid_comm_p0->IsMember() )
  {
    GridRhs_P0_x = new TParVector3D (grid_comm_p0, N_Grid_P0, grid_rhs_p0);
    GridRhs_P0_y = new TParVector3D (grid_comm_p0, N_Grid_P0, grid_rhs_p0 +   N_Grid_P0);
    GridRhs_P0_z = new TParVector3D (grid_comm_p0, N_Grid_P0, grid_rhs_p0 + 2*N_Grid_P0);
  
    GridSol_P0_x = new TParVector3D (grid_comm_p0, N_Grid_P0, grid_sol_p0);
    GridSol_P0_y = new TParVector3D (grid_comm_p0, N_Grid_P0, grid_sol_p0 +   N_Grid_P0);
    GridSol_P0_z = new TParVector3D (grid_comm_p0, N_Grid_P0, grid_sol_p0 + 2*N_Grid_P0);
  }
  else
  {
    GridRhs_P0_x = NULL;
    GridRhs_P0_y = NULL;
    GridRhs_P0_z = NULL;  
    GridSol_P0_x = NULL;
    GridSol_P0_y = NULL;
    GridSol_P0_z = NULL;
  }
  
  if ( grid_comm_p1->IsMember() )
  {
    GridRhs_P1_x = new TParVector3D (grid_comm_p1, N_Grid_P1, grid_rhs_p1);
    GridRhs_P1_y = new TParVector3D (grid_comm_p1, N_Grid_P1, grid_rhs_p1 +   N_Grid_P1);
    GridRhs_P1_z = new TParVector3D (grid_comm_p1, N_Grid_P1, grid_rhs_p1 + 2*N_Grid_P1);
    
    GridSol_P1_x = new TParVector3D (grid_comm_p1, N_Grid_P1, grid_sol_p1);
    GridSol_P1_y = new TParVector3D (grid_comm_p1, N_Grid_P1, grid_sol_p1 +   N_Grid_P1);
    GridSol_P1_z = new TParVector3D (grid_comm_p1, N_Grid_P1, grid_sol_p1 + 2*N_Grid_P1);
  }
  else
  {
    GridRhs_P1_x = NULL;
    GridRhs_P1_y = NULL;
    GridRhs_P1_z = NULL;    
    GridSol_P1_x = NULL;
    GridSol_P1_y = NULL;
    GridSol_P1_z = NULL;
  }
  
  /// fe functions
  if ( velocity_comm->IsMember() )
  {
    // velocity
    u = new TFEVectFunct3D (velocity_space, UName, Desc, sol, N_U, 3);
    ux = u->GetComponent(0);
    uy = u->GetComponent(1);
    uz = u->GetComponent(2);
    
    // pressure
    p = new TFEFunction3D (pressure_space, PName, Desc, sol+3*N_U, N_P);
    
    // grid velocity
    g = new TFEVectFunct3D (grid_space, GName, Desc, grid_velo, N_Grid, 3);
    gx = g->GetComponent(0);
    gy = g->GetComponent(1);
    gz = g->GetComponent(2);
  }
  else
  {
    u = NULL;
    ux = NULL;
    uy = NULL;
    uz = NULL;
    p = NULL;
    g = NULL;
    gx = NULL;
    gy = NULL;
    gz = NULL;
  }
  
  if ( grid_comm_p0->IsMember() )
    g_p0 = new TFEVectFunct3D (grid_space_P0, Name, Desc, grid_sol_p0, N_Grid_P0, 3);
  else
    g_p0 = NULL;
  
  if ( grid_comm_p1->IsMember() )
    g_p1 = new TFEVectFunct3D (grid_space_P1, Name, Desc, grid_sol_p1, N_Grid_P1, 3);
  else
    g_p1 = NULL;
  
  // grid position
  if ( velocity_comm->IsMember() )
    grid = new TFEVectFunct3D(grid_space, GName, Desc, grid_pos, N_Grid, 3);
  else
    grid = NULL;
  
  // add to Output
  if ( velocity_comm->IsMember() )
  {
    Output->AddFEVectFunct(u);
    Output->AddFEFunction(p);
    Output->AddFEVectFunct(g);
  }
  
  /// structures and matrices
  // structures
  if ( velocity_comm->IsMember() )
  {
    sqstructureA = new TSquareStructure3D (velocity_space);
    structureB = new TStructure3D (pressure_space, velocity_space);
    structureBT = new TStructure3D (velocity_space, pressure_space);
    
    sqstructureA->Sort();
    structureB->Sort();
    structureBT->Sort();
    
    
  }
  else
  {
    sqstructureA = NULL;
    structureB = NULL;
    structureBT = NULL;
  }
  
  if ( grid_comm_p0->IsMember() )
  {
    sqstructureGrid_P0 = new TSquareStructure3D (grid_space_P0);
    sqstructureGrid_P0->Sort();
  }
  else
    sqstructureGrid_P0 = NULL;
  
  if ( grid_comm_p1->IsMember() )
  {
    sqstructureGrid_P1 = new TSquareStructure3D (grid_space_P1);
    sqstructureGrid_P1->Sort();
  }
  else
    sqstructureGrid_P1 = NULL; 
  
  // matrices
  if ( velocity_comm->IsWorker() )
  {
    A11 = new TSquareMatrix3D (sqstructureA);
    A12 = new TSquareMatrix3D (sqstructureA);
    A13 = new TSquareMatrix3D (sqstructureA);
    A21 = new TSquareMatrix3D (sqstructureA);
    A22 = new TSquareMatrix3D (sqstructureA);
    A23 = new TSquareMatrix3D (sqstructureA);
    A31 = new TSquareMatrix3D (sqstructureA);
    A32 = new TSquareMatrix3D (sqstructureA);
    A33 = new TSquareMatrix3D (sqstructureA);
    
    M11 = new TSquareMatrix3D (sqstructureA);
    M22 = new TSquareMatrix3D (sqstructureA);
    M33 = new TSquareMatrix3D (sqstructureA);
        
    B1 = new TMatrix3D (structureB);
    B2 = new TMatrix3D (structureB);
    B3 = new TMatrix3D (structureB);
    B1T = new TMatrix3D (structureBT);
    B2T = new TMatrix3D (structureBT);
    B3T = new TMatrix3D (structureBT);
  }
  else
  {
    A11 = A12 = A13 = NULL;
    A21 = A22 = A23 = NULL;
    A31 = A32 = A33 = NULL;
    M11 = M22 = M33 = NULL;
    
    B1 = B2 = B3 = NULL;
    B1T = B2T = B3T = NULL;
  }
  
  if ( grid_comm_p0->IsWorker() )
  {
    G11_P0 = new TSquareMatrix3D (sqstructureGrid_P0);
    G12_P0 = new TSquareMatrix3D (sqstructureGrid_P0);
    G13_P0 = new TSquareMatrix3D (sqstructureGrid_P0);
    G21_P0 = new TSquareMatrix3D (sqstructureGrid_P0);
    G22_P0 = new TSquareMatrix3D (sqstructureGrid_P0);
    G23_P0 = new TSquareMatrix3D (sqstructureGrid_P0);
    G31_P0 = new TSquareMatrix3D (sqstructureGrid_P0);
    G32_P0 = new TSquareMatrix3D (sqstructureGrid_P0);
    G33_P0 = new TSquareMatrix3D (sqstructureGrid_P0);
  }
  else
  {
    G11_P0 = G12_P0 = G13_P0 = NULL;
    G21_P0 = G22_P0 = G23_P0 = NULL;
    G31_P0 = G32_P0 = G33_P0 = NULL;
  }
  
  if ( grid_comm_p1->IsWorker() )
  {
    G11_P1 = new TSquareMatrix3D (sqstructureGrid_P1);
    G12_P1 = new TSquareMatrix3D (sqstructureGrid_P1);
    G13_P1 = new TSquareMatrix3D (sqstructureGrid_P1);
    G21_P1 = new TSquareMatrix3D (sqstructureGrid_P1);
    G22_P1 = new TSquareMatrix3D (sqstructureGrid_P1);
    G23_P1 = new TSquareMatrix3D (sqstructureGrid_P1);
    G31_P1 = new TSquareMatrix3D (sqstructureGrid_P1);
    G32_P1 = new TSquareMatrix3D (sqstructureGrid_P1);
    G33_P1 = new TSquareMatrix3D (sqstructureGrid_P1);
  }
  else
  {
    G11_P1 = G12_P1 = G13_P1 = NULL;
    G21_P1 = G22_P1 = G23_P1 = NULL;
    G31_P1 = G32_P1 = G33_P1 = NULL;
  }
  
  // init Solver
  if ( velocity_comm->IsMember() )
  {
    Solver = new TMumpsSolver (velocity_comm, pressure_comm, sqstructureA, structureB, structureBT);
  }
  else
    Solver = NULL;

  if ( 0 < (int) TDatabase::ParamDB->P1)
  {
    if ( grid_comm_p0->IsMember() )
      Solver_Grid_P0 = new TMumpsSolver (grid_comm_p0, sqstructureGrid_P0);
    else
      Solver_Grid_P0 = NULL;
    
    if ( grid_comm_p1->IsMember() )
      Solver_Grid_P1 = new TMumpsSolver (grid_comm_p1, sqstructureGrid_P1);
    else
      Solver_Grid_P1 = NULL;
  }
  else
  {
    Solver_Grid_P0 = NULL;
    Solver_Grid_P1 = NULL;
  }
  
  if ( velocity_comm->IsWorker() )
  {
    // aux object;
    FEFCT[0] = ux; FEFCT[1] = uy; FEFCT[2] = uz;
    FEFCT[3] = gx; FEFCT[4] = gy; FEFCT[5] = gz;
    aux = new TAuxParam3D (MovingTimeNSN_FESpacesVelo, MovingTimeNSN_FctVelo,
			    MovingTimeNSN_ParamFctVelo, MovingTimeNSN_FEValuesVelo,
			    NULL, FEFCT, MovingTimeNSFctVelo_TwoPhase, MovingTimeNSFEFctIndexVelo,
			    MovingTimeNSFEMultiIndexVelo, MovingTimeNSN_ParamsVelo,
			    MovingTimeNSBeginParamVelo);
  }
  else
    aux = NULL;
  
  // init grid
  if ( velocity_comm->IsMember() )
  {
    grid->GridToData();
  }
  
  if ( grid_comm_p1->IsMember() )
  {
    // write normals into vertices
    SetNormals(Coll_P1, N_SurfaceJoints_P1, CellNumbers_P1, JointNumbers_P1);
  }
  
  BEGIN_SEQ
  // print some informations
  OutPut (endl);
  OutPut ("Velocity DOF: " << 3*N_U << endl);
  OutPut ("Pressure DOF: " << N_P << endl);
  OutPut ("Grid phase 0 DOF: " << 3*N_Grid_P0 << endl);
  OutPut ("Grid phase 1 DOF: " << 3*N_Grid_P1 << endl);
  OutPut ("Number of unknowns: " << N_Unknowns);
  OutPut (" + " << 3*N_Grid_P0 << " + " << 3*N_Grid_P1 << endl);
  END_SEQ
}

void WriteSolution(const char *basename, int number, TOutput3D *Output)
{
  std::ostringstream os;
  
  os << basename << ".";
  os.width(5);
  os.fill('0');
  os << number << ".vtk" << ends;
  
  Output->WriteVtk(os.str().c_str());
}

void PrintMemory()
{
  int bmem = GetMemory();
  int kmem = bmem / 1024;
  int mmem = kmem / 1024;
  
  char mbyte[] = " MByte)";
  char kbyte[] = " KByte)";
  
  int print = mmem == 0 ? kmem : mmem;
  char *p   = mmem == 0 ? kbyte : mbyte;
  
  cout << "Memory used: " << bmem << " (" << print << p << endl;
}

void PrintMemory(int bmem)
{
//   int bmem = GetMemory();
  int kmem = bmem / 1024;
  int mmem = kmem / 1024;
  
  char mbyte[] = " MByte)";
  char kbyte[] = " KByte)";
  
  int print = mmem == 0 ? kmem : mmem;
  char *p   = mmem == 0 ? kbyte : mbyte;
  
  cout << "Memory used: " << bmem << " (" << print << p << endl;
}

void GetFacets(MUTABLE_DATA *Data, int &N_Points, double *&Points, int *&Facets)
{
  int MaxLen, len, CellNr, JointNr, ivar, counter;
  const int *TmpFV, *TmpLen, *indices;
  double x, y, z;
  TVertex *Vertex;
  TBaseCell *Cell;
  
// reset clipboard
  for (int i=0;i<Data->N_SurfaceJoints_P1;++i)
  {
    CellNr = Data->CellNumbers_P1[i];
    JointNr = Data->JointNumbers_P1[i];
    
    Cell = Data->Coll_P1->GetCell(CellNr);
    
    Cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
    
    len = TmpLen[JointNr];
    indices = TmpFV+JointNr*MaxLen;
    for (int j=0;j<len;++j)
    {
      Cell->GetVertex(indices[j])->SetClipBoard(-1);
    }
  } // end for i
  
  // number vertices on surface phase 1
  counter = 0;
  for (int i=0;i<Data->N_SurfaceJoints_P1;++i)
  {
    CellNr = Data->CellNumbers_P1[i];
    JointNr = Data->JointNumbers_P1[i];
    
    Cell = Data->Coll_P1->GetCell(CellNr);
    
    Cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
    
    len = TmpLen[JointNr];
    indices = TmpFV+JointNr*MaxLen;
    for (int j=0;j<len;++j)
    {
      Vertex = Cell->GetVertex(indices[j]);
      
      if ( Vertex->GetClipBoard() == -1 )
      {
	Vertex->SetClipBoard(counter++);
      }
    }
  } // end for i
  
  N_Points = counter;
  Points = new double [3*N_Points];
  Facets = new int [3*Data->N_SurfaceJoints_P1];
  
  for (int i=0;i<Data->N_SurfaceJoints_P1;++i)
  {
    CellNr = Data->CellNumbers_P1[i];
    JointNr = Data->JointNumbers_P1[i];
    
    Cell = Data->Coll_P1->GetCell(CellNr);
    
    Cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
    
    len = TmpLen[JointNr];
    indices = TmpFV+JointNr*MaxLen;
    for (int j=0;j<len;++j)
    {
      Vertex = Cell->GetVertex(indices[j]);
      Vertex->GetCoords(x, y, z);
      ivar = Vertex->GetClipBoard();
      
      Points[3*ivar  ] = x;
      Points[3*ivar+1] = y;
      Points[3*ivar+2] = z;
      
      Facets[3*i+j] = ivar;
    }
    
    ivar = Facets[3*i];
    Facets[3*i] = Facets[3*i+1];
    Facets[3*i+1] = ivar;
  } // end for i
}

void Remesh(MPI_Comm comm, TDomain *Domain, MUTABLE_DATA *Data)
{
  int rank, size, MaxCpV;
  int N_Points, N_Facets, *Facets;
  int iBuff[10];
  double Regions[10], x, y, z, *Points;
  double dBuff[10], StartX, StartY, StartZ, BoundX, BoundY, BoundZ;
  double mem;
  
  mem = GetMemory();
  
  TBaseCell *Cell;
  
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  
  ROOTOUT (endl);
  ROOTOUT ("=== REMESHING ==================== " << endl);
  
  // save old data
  MUTABLE_DATA Data_old;
  Data_old = *Data;
  
  // already free some memory
  delete Data_old.sqstructureA; delete Data_old.sqstructureGrid_P0;
  delete Data_old.sqstructureGrid_P1;
  delete Data_old.structureB; delete Data_old.structureBT;
  
  delete [] Data_old.rhs;
  delete [] Data_old.res;
  delete [] Data_old.grid_velo;
  delete [] Data_old.grid_pos;
  delete [] Data_old.grid_rhs_p0;
  delete [] Data_old.grid_rhs_p1;
  delete [] Data_old.grid_sol_p0;
  delete [] Data_old.grid_sol_p1;
  delete [] Data_old.aux_array;
  
  delete [] Data_old.phase_marker->GetValues();
  delete Data_old.phase_marker;
  
  if ( rank == 0 ) // root
  {
    // get surfaces
    GetFacets(Data, N_Points, Points, Facets);
    
    // get regions
    Cell = Data->Coll_P0->GetCell(0),
    Regions[0] = 0;
    Regions[1] = 0;
    Regions[2] = 0;
    Regions[3] = 0;
    Regions[4] = 0;
    
    for (int i=0;i<4;++i)
    {
      Cell->GetVertex(i)->GetCoords(x, y, z);
      
      Regions[0] += 0.25*x;
      Regions[1] += 0.25*y;
      Regions[2] += 0.25*z;
    }
    
    Cell = Data->Coll_P1->GetCell(0);
    
    Regions[5] = 0;
    Regions[6] = 0;
    Regions[7] = 0;
    Regions[8] = 1;
    Regions[9] = 0;
    
    for (int i=0;i<4;++i)
    {
      Cell->GetVertex(i)->GetCoords(x, y, z);
      
      Regions[5] += 0.25*x;
      Regions[6] += 0.25*y;
      Regions[7] += 0.25*z;
    }
    
    N_Facets = Data->N_SurfaceJoints_P1;
    
    Domain->GetBoundBox(dBuff[0], dBuff[1], dBuff[2],
			dBuff[3], dBuff[4], dBuff[5]);	
    
    iBuff[0] = N_Points;
    iBuff[1] = N_Facets;
    
    MPI_Bcast(iBuff, 2, MPI_INT, 0, comm);
    MPI_Bcast(Regions, 10, MPI_DOUBLE, 0, comm);
    MPI_Bcast(dBuff, 6, MPI_DOUBLE, 0, comm);
    
    MPI_Bcast(Points, 3*N_Points, MPI_DOUBLE, 0, comm);
    MPI_Bcast(Facets, 3*N_Facets, MPI_INT, 0, comm);
  }
  else // worker
  {
    // free some memory
    delete Data_old.A11; delete Data_old.A12; delete Data_old.A13;
    delete Data_old.A21; delete Data_old.A22; delete Data_old.A23;
    delete Data_old.A31; delete Data_old.A32; delete Data_old.A33;
    
    delete Data_old.G11_P0; delete Data_old.G12_P0; delete Data_old.G13_P0;
    delete Data_old.G21_P0; delete Data_old.G22_P0; delete Data_old.G23_P0;
    delete Data_old.G31_P0; delete Data_old.G32_P0; delete Data_old.G33_P0;
    
    delete Data_old.G11_P1; delete Data_old.G12_P1; delete Data_old.G13_P1;
    delete Data_old.G21_P1; delete Data_old.G22_P1; delete Data_old.G23_P1;
    delete Data_old.G31_P1; delete Data_old.G32_P1; delete Data_old.G33_P1;
    
    delete Data_old.B1; delete Data_old.B2; delete Data_old.B3;
    delete Data_old.B1T; delete Data_old.B2T; delete Data_old.B3T;
    
    delete Data_old.aux;
        
    MPI_Bcast(iBuff, 2, MPI_INT, 0, comm);
    MPI_Bcast(Regions, 10, MPI_DOUBLE, 0, comm);
    MPI_Bcast(dBuff, 6, MPI_DOUBLE, 0, comm);
    
    N_Points = iBuff[0];
    N_Facets = iBuff[1];
    
    Points = new double [3*N_Points];
    Facets = new int [3*N_Facets];
    
    MPI_Bcast(Points, 3*N_Points, MPI_DOUBLE, 0, comm);
    MPI_Bcast(Facets, 3*N_Facets, MPI_INT, 0, comm);
    
    Domain->SetBoundBox(dBuff[0], dBuff[1], dBuff[2],
			dBuff[3], dBuff[4], dBuff[5]);
  }
  
  TTetGenMeshLoader::TDummyDomain DummyDomain;
  TTetGenMeshLoader MeshLoader;
    
  // create new grid on all ranks
  MeshLoader.Generate(N_Points, Points, N_Facets, Facets, 2, Regions,
		      Domain, &DummyDomain);  
		      
  // mesh partition on all ranks
  Partition_Mesh(comm, Domain, MaxCpV);
  
  // create new data
  Data->Create(Domain, comm);
  
  /// interpolate from old to new mesh on root
  if ( rank == 0 ) // root
  {
    double t1, t2;
    
    t1 = GetTime();
    
    TOctTree::TBoundBox BoundBox;
    BoundBox.StartX = DummyDomain.StartX;
    BoundBox.BoundX = DummyDomain.BoundX;
    BoundBox.StartY = DummyDomain.StartY;
    BoundBox.BoundY = DummyDomain.BoundY;
    BoundBox.StartZ = DummyDomain.StartZ;
    BoundBox.BoundZ = DummyDomain.BoundZ;
    
    TOctTree OctTree(Data_old.Coll, &BoundBox);
    
    Data->ux->Interpolate(Data_old.ux, &OctTree);
    Data->uy->Interpolate(Data_old.uy, &OctTree);
    Data->uz->Interpolate(Data_old.uz, &OctTree);
    Data->p->Interpolate(Data_old.p, &OctTree);
    t2 = GetTime();
  
    OutPut("time for interpolation: " << t2-t1 << endl);
  }
  
  // distribute interpolated solution to worker
  Data->Sol_ux->ScatterFromRoot();
  Data->Sol_uy->ScatterFromRoot();
  Data->Sol_uz->ScatterFromRoot();
  Data->Sol_p->ScatterFromRoot();
  
  // free rest of memory
  FreeDomain(&DummyDomain);
  
  delete Data_old.Coll;
  delete Data_old.Coll_P0;
  delete Data_old.Coll_P1;
  
  delete [] Data_old.GlobalCellNo_P0;
  delete [] Data_old.GlobalCellNo_P1;
  
  delete [] Data_old.CellNumbers_P0;
  delete [] Data_old.CellNumbers_P1;
  delete [] Data_old.JointNumbers_P0;
  delete [] Data_old.JointNumbers_P1;
  
  delete Data_old.velocity_space;
  delete Data_old.pressure_space;
  delete Data_old.grid_space;
  delete Data_old.grid_space_P0;
  delete Data_old.grid_space_P1;
  
  delete Data_old.velocity_comm;
  delete Data_old.pressure_comm;
  delete Data_old.grid_comm_p0;
  delete Data_old.grid_comm_p1;
  
  delete Data_old.Sol_ux; delete Data_old.Sol_uy; delete Data_old.Sol_uz;
  delete Data_old.Sol_p;
  delete Data_old.Rhs_ux; delete Data_old.Rhs_uy; delete Data_old.Rhs_uz;
  delete Data_old.Rhs_p;
  delete Data_old.Res_ux; delete Data_old.Res_uy; delete Data_old.Res_uz;
  delete Data_old.Res_p;
  delete Data_old.GridRhs_P0_x; delete Data_old.GridRhs_P0_y; delete Data_old.GridRhs_P0_z;
  delete Data_old.GridRhs_P1_x; delete Data_old.GridRhs_P1_y; delete Data_old.GridRhs_P1_z;
  delete Data_old.GridSol_P0_x; delete Data_old.GridSol_P0_y; delete Data_old.GridSol_P0_z;
  delete Data_old.GridSol_P1_x; delete Data_old.GridSol_P1_y; delete Data_old.GridSol_P1_z;
  
  delete Data_old.u; delete Data_old.ux; delete Data_old.uy; delete Data_old.uz;
  delete Data_old.g; delete Data_old.gx; delete Data_old.gy; delete Data_old.gz;
  delete Data_old.p;
  delete Data_old.g_p0; delete Data_old.g_p1;
  delete Data_old.grid;
  
  delete Data_old.Output;
  
  delete Data_old.Solver;
  if ( 0 < (int) TDatabase::ParamDB->P1 )
  {
    delete Data_old.Solver_Grid_P0;
    delete Data_old.Solver_Grid_P1;
  }
  
  ROOTOUT ("=== REMESH DONE ================== " << endl);
  

  delete [] Facets;
  delete [] Points;
  
  BEGIN_SEQ
  OutPut ("memory before: ");
  PrintMemory(mem);
  OutPut ("memory after: ");
  PrintMemory();
  END_SEQ
}

void FreeDomain(TTetGenMeshLoader::TDummyDomain *Domain)
{
  int N_Cells, NeibJoint, N_Vertices, N_, N_BoundParts, ivar;
  TBaseCell **Cells, *Cell, *Neib;
  TVertex *Vertex, **Vertices;
  TBoundPart **BdParts, *BoundPart;
  TBoundComp3D *BoundComp;
  
  N_Cells = Domain->N_RootCells;
  Cells = Domain->CellTree;
  
  // reset clipboard
  for (int i=0;i<N_Cells;++i)
  {
    Cell = Cells[i];
    
    for (int j=0;j<4;++j)
    {
      Cell->GetVertex(j)->SetClipBoard(-1);
    }
  } // end for i
  
  Vertices = new TVertex* [4*N_Cells];
  
  N_Vertices = 0;
  for (int i=0;i<N_Cells;++i)
  {
    Cell = Cells[i];
    
    for (int j=0;j<4;++j)
    {
      Vertex = Cell->GetVertex(j);
      if ( Vertex->GetClipBoard() == -1 )
      {
	Vertex->SetClipBoard(N_Vertices);
	Vertices[N_Vertices] = Vertex;
	
	++N_Vertices;	
      }
    }
    
    delete Cell;
  }
  
  for (int i=0;i<N_Vertices;++i)
  {
    delete Vertices[i];
  }
  
  delete [] Vertices;
  delete [] Cells;
  
  N_BoundParts = Domain->N_BoundParts;
  BdParts = Domain->BdParts;  
  
  for (int i=0;i<N_BoundParts;++i)
  {
    BoundPart = BdParts[i];
    N_ = BoundPart->GetN_BdComps();
    
    for (int j=0;j<N_;++j)
    {
      BoundComp = BoundPart->GetBdComp(j);
      
      delete BoundComp;
    }
    
    delete BoundPart;
  }
  
  delete [] BdParts;
  delete [] Domain->StartBdCompID;
//   delete [] Domain->Interfaces;
}

double CheckForRemesh(TCollection *Coll)
{
  int N_Cells;
  double diameter, ratio, maxratio = 0, meanratio = 0, vol, r;
  double baryx, baryy, baryz, x, y, z, ax, ay, az, bx, by, bz, cx, cy, cz, rx, ry, rz;
  double ox, oy, oz, a, b, c;
  TBaseCell *Cell;
  
  N_Cells = Coll->GetN_Cells();
  
  for (int i=0;i<N_Cells;++i)
  {
    Cell = Coll->GetCell(i);
    
    Cell->GetVertex(0)->GetCoords(ox, oy, oz);
    
    Cell->GetVertex(1)->GetCoords(x, y, z);
    ax = x - ox;
    ay = y - oy;
    az = z - oz;
    
    Cell->GetVertex(2)->GetCoords(x, y, z);
    bx = x - ox;
    by = y - oy;
    bz = z - oz;
    
    Cell->GetVertex(3)->GetCoords(x, y, z);
    cx = x - ox;
    cy = y - oy;
    cz = z - oz;
    
    x = ax * (by*cz - bz*cy);
    y = ay * (bz*cx - bx*cz);
    z = az * (bx*cy - by*cx);
    
    vol = fabs(x+y+z);
    vol /= 6;
    
    a = ax*ax + ay*ay + az*az;
    b = bx*bx + by*by + bz*bz;
    c = cx*cx + cy*cy + cz*cz;
    
    rx = a * (by*cz - bz*cy) + b * (cy*az - cz*ay) + c * (ay*bz - az*by);
    ry = a * (bz*cx - bx*cz) + b * (cz*ax - cx*az) + c * (az*bx - ax*bz);
    rz = a * (bx*cy - by*cx) + b * (cx*ay - cy*ax) + c * (ax*by - ay*bx);
        
    r = sqrt(rx*rx + ry*ry + rz*rz);    
    r /= 12*vol;
    
    ratio = r / Cell->GetShortestEdge();
    
    if ( maxratio < ratio ) maxratio = ratio;
    
    meanratio += ratio;
  }
  
  meanratio /= N_Cells;
  
  OutPut("Max ratio: " << maxratio << endl);
  OutPut("Mean ratio: " << meanratio << endl);
  
//   exit(0);
  
  return maxratio;
}

