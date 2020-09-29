#include <MooNMD.h>

#include <stdlib.h>

// #include "../Examples/TNSE_3D/AnsatzLinConst.h"
#include "../Examples/TNSE_3D/TwoPhase.h"

void PrintMemory();
void MatReset(int N_SQMAT, TSquareMatrix3D **SQMATRICES, int N_MAT, TMatrix3D **MATRICES);
void WriteSolution(const char *basename, int number, TOutput3D *Output);
void CheckMesh(TFESpace3D *fesp);
void SplitGrid_TwoPhase(TDomain *Domain, TCollection *&Coll, TCollection *&Coll_P0,
			TCollection *&Coll_P1, TFEFunction3D *&phase_marker,
			int *&GlobalCellNo_P0, int *&GlobalCellNo_P1,
			int &N_SurfaceJoints_P0, int *&CellNumbers_P0, int *&JointNumbers_P0,
			int &N_SurfaceJoints_P1, int *&CellNumbers_P1, int *&JointNumbers_P1);

void WriteGrid(const char *basename, int number, TDomain *Domain, TCollection *Phase);

void Remesh(TCollection *Coll, TDomain *Domain);
			
int main(int argc, char **argv)
{
  TDatabase Database;
  TFEDatabase3D FEDatabase;
  TDomain Domain;
  
  TCollection *Coll, *Coll_P0, *Coll_P1;
  int *GlobalCellNo_P0, *GlobalCellNo_P1;
  
  TFESpace3D *velocity_space, *pressure_space, *grid_space, *grid_space_P0, *grid_space_P1;
  int pressure_space_code;
  
  TDiscreteForm3D *DiscreteFormGalerkin, *DiscreteFormNLGalerkin;
  TDiscreteForm3D *DiscreteFormRhs, *DiscreteFormGrid;
  
  TFEFunction3D *phase_marker;
  
  if ( argc > 1 )
    Domain.ReadParam(argv[1]);
  else
  {
    cerr << "no readin file given" << endl;
    exit(0);
  }
  
  char Name[] = "name";
  char Desc[] = "desc";
  char UName[] = "u";
  char PName[] = "p";
  char GName[] = "g";
  
  char *SMESH = TDatabase::ParamDB->SMESHFILE;
  char *VTK = TDatabase::ParamDB->VTKBASENAME;
  int NSTYPE = TDatabase::ParamDB->NSTYPE;
  
  if ( NSTYPE != 4 )
  {
    cerr << "Use NSTYPE = 4 !" << endl;
    exit(0);
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
  
  Domain.Tetgen(SMESH);
  
  int N_SurfaceJoints_P0, *CellNumbers_P0, *JointNumbers_P0;
  int N_SurfaceJoints_P1, *CellNumbers_P1, *JointNumbers_P1;
  
  SplitGrid_TwoPhase(&Domain, Coll, Coll_P0, Coll_P1, phase_marker,
		     GlobalCellNo_P0, GlobalCellNo_P1,
		     N_SurfaceJoints_P0, CellNumbers_P0, JointNumbers_P0,
		     N_SurfaceJoints_P1, CellNumbers_P1, JointNumbers_P1);
  
  TOutput3D Output (3, 1, 2, 0, &Domain);
  
  Output.AddFEFunction(phase_marker);
  Output.WriteVtk("phase.vtk");
  
  // discreteforms
  InitializeDiscreteForms(DiscreteFormGalerkin, DiscreteFormNLGalerkin,
			  DiscreteFormRhs, LinCoeffs, NSTYPE);
			  
  DiscreteFormGrid = new TDiscreteForm3D (Name, Desc, ESGridN_Terms, ESGridDerivatives,
					  ESGridSpaceNumbers, ESGridN_Matrices,
					  ESGridN_Rhs, ESGridRowSpace, ESGridColumnSpace,
					  ESGridRhsSpace, ESGridAssemble, GridCoeffs, NULL);
  
  // velocity and pressure space
  GetVelocityAndPressureSpace3D (Coll, BoundCondition, velocity_space, pressure_space,
				 &pressure_space_code,
				 TDatabase::ParamDB->VELOCITY_SPACE,
				 TDatabase::ParamDB->PRESSURE_SPACE);
				 
  // grid spaces
  grid_space = new TFESpace3D (Coll, Name, Desc, GridBoundCond, 1);
  grid_space_P0 = new TFESpace3D (Coll_P0, Name, Desc, GridBoundCond, 1);
  grid_space_P1 = new TFESpace3D (Coll_P1, Name, Desc, GridBoundCond, 1);
  
  CheckMesh(velocity_space);
  
  // count DOF's
  int N_U = velocity_space->GetN_DegreesOfFreedom();
  int N_P = pressure_space->GetN_DegreesOfFreedom();
  int N_Grid_P0 = grid_space_P0->GetN_DegreesOfFreedom();
  int N_Grid_P1 = grid_space_P1->GetN_DegreesOfFreedom();
  int N_Grid = grid_space->GetN_DegreesOfFreedom();
  int N_Unknowns = 3*N_U + N_P;
  int ActiveBound = velocity_space->GetActiveBound();
  
  // alloc arrays
  double *sol = new double [N_Unknowns];
  double *rhs = new double [N_Unknowns];
  double *res = new double [N_Unknowns];
  double *grid_velo = new double [3*N_Grid];
  double *grid_pos = new double [3*N_Grid];
  double *grid_rhs_p0 = new double [3*N_Grid_P0];
  double *grid_rhs_p1 = new double [3*N_Grid_P1];
  double *grid_sol_p0 = new double [3*N_Grid_P0];
  double *grid_sol_p1 = new double [3*N_Grid_P1];
  double *aux_array = new double [N_U];
  double *dptr1, *dptr2;
  
  memset(sol, 0, N_Unknowns*sizeof(sol[0]));
  memset(rhs, 0, N_Unknowns*sizeof(rhs[0]));
  memset(grid_velo, 0, 3*N_Grid*sizeof(grid_velo[0]));
  memset(grid_pos, 0, 3*N_Grid*sizeof(grid_pos[0]));
  memset(grid_rhs_p0, 0, 3*N_Grid_P0*sizeof(grid_rhs_p0[0]));
  memset(grid_rhs_p1, 0, 3*N_Grid_P1*sizeof(grid_rhs_p1[0]));
  
  /// fe functions
  // velocity
  TFEVectFunct3D *u = new TFEVectFunct3D (velocity_space, UName, Desc, sol, N_U, 3);
  TFEFunction3D *ux = u->GetComponent(0);
  TFEFunction3D *uy = u->GetComponent(1);
  TFEFunction3D *uz = u->GetComponent(2);
  
  // pressure
  TFEFunction3D *p = new TFEFunction3D (pressure_space, PName, Desc, sol+3*N_U, N_P);
  
  // grid velocity
  TFEVectFunct3D *g = new TFEVectFunct3D (grid_space, GName, Desc, grid_velo, N_Grid, 3);
  TFEFunction3D *gx = g->GetComponent(0);
  TFEFunction3D *gy = g->GetComponent(1);
  TFEFunction3D *gz = g->GetComponent(2);
  
  TFEVectFunct3D *g_p0 = new TFEVectFunct3D (grid_space_P0, Name, Desc, grid_sol_p0, N_Grid_P0, 3);
  TFEVectFunct3D *g_p1 = new TFEVectFunct3D (grid_space_P1, Name, Desc, grid_sol_p1, N_Grid_P1, 3);
  
  // grid position
  TFEVectFunct3D *grid = new TFEVectFunct3D(grid_space, GName, Desc, grid_pos, N_Grid, 3);
  
  // add to Output
  Output.AddFEVectFunct(u);
  Output.AddFEFunction(p);
  Output.AddFEVectFunct(g);
  
  /// structures and matrices
  // structures
  TSquareStructure3D sqstructureA (velocity_space);
  TSquareStructure3D sqstructureGrid_P0 (grid_space_P0);
  TSquareStructure3D sqstructureGrid_P1 (grid_space_P1);
  TStructure3D structureB (pressure_space, velocity_space);
  TStructure3D structureBT (velocity_space, pressure_space);
  sqstructureA.Sort();
  sqstructureGrid_P0.Sort();
  sqstructureGrid_P1.Sort();
  structureB.Sort();
  structureBT.Sort();
  
  // matrices
  TSquareMatrix3D A11 (&sqstructureA);
  TSquareMatrix3D A12 (&sqstructureA);
  TSquareMatrix3D A13 (&sqstructureA);
  TSquareMatrix3D A21 (&sqstructureA);
  TSquareMatrix3D A22 (&sqstructureA);
  TSquareMatrix3D A23 (&sqstructureA);
  TSquareMatrix3D A31 (&sqstructureA);
  TSquareMatrix3D A32 (&sqstructureA);
  TSquareMatrix3D A33 (&sqstructureA);
  
  TSquareMatrix3D M11 (&sqstructureA);
  TSquareMatrix3D M22 (&sqstructureA);
  TSquareMatrix3D M33 (&sqstructureA);
  
  TSquareMatrix3D G11_P0 (&sqstructureGrid_P0);
  TSquareMatrix3D G12_P0 (&sqstructureGrid_P0);
  TSquareMatrix3D G13_P0 (&sqstructureGrid_P0);
  TSquareMatrix3D G21_P0 (&sqstructureGrid_P0);
  TSquareMatrix3D G22_P0 (&sqstructureGrid_P0);
  TSquareMatrix3D G23_P0 (&sqstructureGrid_P0);
  TSquareMatrix3D G31_P0 (&sqstructureGrid_P0);
  TSquareMatrix3D G32_P0 (&sqstructureGrid_P0);
  TSquareMatrix3D G33_P0 (&sqstructureGrid_P0);
  
  TSquareMatrix3D G11_P1 (&sqstructureGrid_P1);
  TSquareMatrix3D G12_P1 (&sqstructureGrid_P1);
  TSquareMatrix3D G13_P1 (&sqstructureGrid_P1);
  TSquareMatrix3D G21_P1 (&sqstructureGrid_P1);
  TSquareMatrix3D G22_P1 (&sqstructureGrid_P1);
  TSquareMatrix3D G23_P1 (&sqstructureGrid_P1);
  TSquareMatrix3D G31_P1 (&sqstructureGrid_P1);
  TSquareMatrix3D G32_P1 (&sqstructureGrid_P1);
  TSquareMatrix3D G33_P1 (&sqstructureGrid_P1);
  
  TMatrix3D B1 (&structureB);
  TMatrix3D B2 (&structureB);
  TMatrix3D B3 (&structureB);
  TMatrix3D B1T (&structureBT);
  TMatrix3D B2T (&structureBT);
  TMatrix3D B3T (&structureBT);
  
  /// prepare assemble
  TSquareMatrix3D *SQMATRICES[12];
  TSquareMatrix3D *SQMATRICES_GRID[12];
  TSquareMatrix3D *sqmptr1;
  TMatrix3D *MATRICES[6];
  TFESpace3D *FESP[4], *FESPRHS[4];
  TFEFunction3D *FEFCT[6];
  int N_MAT, N_SQMAT, N_RHS, N_FESP;
  double *RHS[3], *RHS_GRID_P0[3], *RHS_GRID_P1[3];
  
  RHS[0] = rhs        ;
  RHS[1] = rhs +   N_U;
  RHS[2] = rhs + 2*N_U;
  RHS_GRID_P0[0] = grid_rhs_p0              ;
  RHS_GRID_P0[1] = grid_rhs_p0 +   N_Grid_P0;
  RHS_GRID_P0[2] = grid_rhs_p0 + 2*N_Grid_P0;
  RHS_GRID_P1[0] = grid_rhs_p1              ;
  RHS_GRID_P1[1] = grid_rhs_p1 +   N_Grid_P1;
  RHS_GRID_P1[2] = grid_rhs_p1 + 2*N_Grid_P1;
  N_RHS = 3;
  
  FESPRHS[0] = velocity_space;
  FESPRHS[1] = velocity_space;
  FESPRHS[2] = velocity_space;
  
  // aux object;
  FEFCT[0] = ux; FEFCT[1] = uy; FEFCT[2] = uz;
  FEFCT[3] = gx; FEFCT[4] = gy; FEFCT[5] = gz;
  TAuxParam3D aux (MovingTimeNSN_FESpacesVelo, MovingTimeNSN_FctVelo,
		   MovingTimeNSN_ParamFctVelo, MovingTimeNSN_FEValuesVelo,
		   FESP, FEFCT, MovingTimeNSFctVelo_TwoPhase, MovingTimeNSFEFctIndexVelo,
		   MovingTimeNSFEMultiIndexVelo, MovingTimeNSN_ParamsVelo,
		   MovingTimeNSBeginParamVelo);
		   
  TAuxParam3D aux_dummy (0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
  
  Remesh(Coll, &Domain);
exit(0);
  
  /// time loop (back euler);
  int img = 0;
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  
  ux->Interpolate(InitialU1);
  uy->Interpolate(InitialU2);
  uz->Interpolate(InitialU3);
  p->Interpolate(InitialP);
  
  // print some informations
  OutPut (endl);
  OutPut ("Velocity DOF: " << 3*N_U << endl);
  OutPut ("Pressure DOF: " << N_P << endl);
  OutPut ("Grid phase 0 DOF: " << 3*N_Grid_P0 << endl);
  OutPut ("Grid phase 1 DOF: " << 3*N_Grid_P1 << endl);
  OutPut ("Number of unknowns: " << N_Unknowns);
  OutPut (" + " << 3*N_Grid_P0 << " + " << 3*N_Grid_P1 << endl);
 
  OutPut(endl);
  OutPut("RE_NR: " << TDatabase::ParamDB->RE_NR << endl);
  OutPut("WE_NR: " << TDatabase::ParamDB->WB_NR << endl);
  
  OutPut(endl);
  OutPut("STARTTIME: " << TDatabase::TimeDB->CURRENTTIME << endl);
  OutPut("ENDTIME: " << TDatabase::TimeDB->ENDTIME << endl);
  OutPut("TIMESTEPLENGTH: " << tau << endl);
  double Volumen = GetVolume(Coll);
  double Volumen_P0 = GetVolume(Coll_P0);
  double Volumen_P1 = GetVolume(Coll_P1);
  OutPut("Volume: " << Volumen << endl);
  OutPut("Volume phase 0: " << Volumen_P0 << endl);
  OutPut("Volume phase 1: " << Volumen_P1 << endl);
  OutPut("Density ratio: " << TDatabase::ParamDB->P6 << endl);
  OutPut("Viscosity ratio: " << TDatabase::ParamDB->P5 << endl << endl);
  
  WriteSolution(VTK, img++, &Output);
  
  // init grid
  grid->GridToData();
  
  // initialise directsolver
  TPardisoSolver Solver;
  
  PrintMemory();
  
  /// time loop begin
  for(int timestep=0;TDatabase::TimeDB->CURRENTTIME <= TDatabase::TimeDB->ENDTIME;++timestep)
  {
    TDatabase::TimeDB->CURRENTTIME += tau;
    
    OutPut (endl << "CURRENTTIME: " << TDatabase::TimeDB->CURRENTTIME << endl);
    
    // reset rhs
    memset(rhs, 0, N_Unknowns*sizeof(rhs[0]));
    memset(grid_rhs_p0, 0, 3*N_Grid_P0*sizeof(grid_rhs_p0[0]));
    memset(grid_rhs_p1, 0, 3*N_Grid_P1*sizeof(grid_rhs_p1[0]));
    
    /// assemble systems
    FESP[0] = velocity_space;
    FESP[1] = pressure_space;
    FESP[2] = grid_space;
    N_FESP = 3;
    
    SQMATRICES[0] = &A11;
    SQMATRICES[1] = &A12;
    SQMATRICES[2] = &A13;
    SQMATRICES[3] = &A21;
    SQMATRICES[4] = &A22;
    SQMATRICES[5] = &A23;
    SQMATRICES[6] = &A31;
    SQMATRICES[7] = &A32;
    SQMATRICES[8] = &A33;
    SQMATRICES[9] = &M11;
    SQMATRICES[10] = &M22;
    SQMATRICES[11] = &M33;
    N_SQMAT = 12;
    
    MATRICES[0] = &B1;
    MATRICES[1] = &B2;
    MATRICES[2] = &B3;
    MATRICES[3] = &B1T;
    MATRICES[4] = &B2T;
    MATRICES[5] = &B3T;
    N_MAT = 6;
    
    MatReset(N_SQMAT, SQMATRICES, N_MAT, MATRICES);
    
    // assemble
    Assemble3D (N_FESP, FESP,
		N_SQMAT, SQMATRICES,
		N_MAT, MATRICES,
		N_RHS, RHS, FESPRHS,
		DiscreteFormGalerkin,
		BOUNDCOND, BOUNDVALUE, &aux);
    
    // grid phase 0;
    FESP[0] = grid_space_P0;
    N_FESP = 1;
		
    SQMATRICES[0] = &G11_P0;
    SQMATRICES[1] = &G12_P0;
    SQMATRICES[2] = &G13_P0;
    SQMATRICES[3] = &G21_P0;
    SQMATRICES[4] = &G22_P0;
    SQMATRICES[5] = &G23_P0;
    SQMATRICES[6] = &G31_P0;
    SQMATRICES[7] = &G32_P0;
    SQMATRICES[8] = &G33_P0;
    N_SQMAT = 9;
    N_MAT = 0;
    
    MatReset(N_SQMAT, SQMATRICES, N_MAT, NULL);
    
    Assemble3D (N_FESP, FESP,
		N_SQMAT, SQMATRICES,
		N_MAT, NULL,
		0, NULL, NULL,
		DiscreteFormGrid,
		GRIDBOUNDCOND, GRIDBOUNDVALUE, &aux_dummy);
    
    // grid phase 1
    FESP[0] = grid_space_P1;
    N_FESP = 1;
		
    SQMATRICES[0] = &G11_P1;
    SQMATRICES[1] = &G12_P1;
    SQMATRICES[2] = &G13_P1;
    SQMATRICES[3] = &G21_P1;
    SQMATRICES[4] = &G22_P1;
    SQMATRICES[5] = &G23_P1;
    SQMATRICES[6] = &G31_P1;
    SQMATRICES[7] = &G32_P1;
    SQMATRICES[8] = &G33_P1;
    N_SQMAT = 9;
    N_MAT = 0;
    
    MatReset(N_SQMAT, SQMATRICES, N_MAT, NULL);
    
    Assemble3D (N_FESP, FESP,
		N_SQMAT, SQMATRICES,
		N_MAT, NULL,
		0, NULL, NULL,
		DiscreteFormGrid,
		GRIDBOUNDCOND, GRIDBOUNDVALUE, &aux_dummy);
    
    // scale BT's
    Dscal(B1T.GetN_Entries(), tau, B1T.GetEntries());
    Dscal(B2T.GetN_Entries(), tau, B2T.GetEntries());
    Dscal(B3T.GetN_Entries(), tau, B3T.GetEntries());
    
    // add M*u_old to rhs and scale rhs
    SQMATRICES[0] = &M11;
    SQMATRICES[1] = &M22;
    SQMATRICES[2] = &M33;
    for (int i=0;i<3;++i)
    {
      dptr1 = sol + i*N_U;
      dptr2 = rhs + i*N_U;
      
      MatVectActive(SQMATRICES[i], dptr1, aux_array);
      Dscal(ActiveBound, tau, dptr2);
      Daxpy(ActiveBound, 1.0, aux_array, dptr2);
    }
    
    // TODO: add surface integrals to M here
    
    // add M block to A blocks and scale A blocks
    MatAdd2(&A11, &M11, tau);
    MatAdd2(&A22, &M22, tau);
    MatAdd2(&A33, &M33, tau);
    
    SQMATRICES[0] = &A12;
    SQMATRICES[1] = &A13;
    SQMATRICES[2] = &A21;
    SQMATRICES[3] = &A23;
    SQMATRICES[4] = &A31;
    SQMATRICES[5] = &A32;
    
    for (int i=0;i<6;++i)
      Dscal(SQMATRICES[i]->GetN_Entries(), tau, SQMATRICES[i]->GetEntries());
    
    /// solve system
    Solver.SetMatrix(&A11, &A12, &A13, &A21, &A22, &A23, &A31, &A32, &A33,
		     &B1T, &B2T, &B3T, &B1, &B2, &B3);
		     
    if ( timestep % 10 == 0 )
    {
      Solver.Analyse();
      Solver.Factorize();
      Solver.Solve(sol, rhs);
    }
    else
    {
      Solver.FactorizeSolve(sol, rhs);
    }
//     DirectSolver(&A11, &A12, &A13, &A21, &A22, &A23, &A31, &A32, &A33,
// 		 &B1T, &B2T, &B3T, &B1, &B2, &B3, rhs, sol, 3);
    
    // write normals into vertices
    SetNormals(Coll_P1, N_SurfaceJoints_P1, CellNumbers_P1, JointNumbers_P1);
    
    /// start nonlinear loop
    for (int iter=0;iter<TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;++iter)
    {
      // solve grid equation
      memset(grid_rhs_p0, 0, 3*N_Grid_P0*sizeof(grid_rhs_p0[0]));
      memset(grid_rhs_p1, 0, 3*N_Grid_P1*sizeof(grid_rhs_p1[0]));
      
      SetGridRhs(u, grid_space_P0, N_SurfaceJoints_P0, CellNumbers_P0, JointNumbers_P0,
		 GlobalCellNo_P0, tau, grid_rhs_p0, grid_rhs_p0+N_Grid_P0, grid_rhs_p0+2*N_Grid_P0);
      SetGridRhs(u, grid_space_P1, N_SurfaceJoints_P1, CellNumbers_P1, JointNumbers_P1,
		 GlobalCellNo_P1, tau, grid_rhs_p1, grid_rhs_p1+N_Grid_P1, grid_rhs_p1+2*N_Grid_P1);
		 
      SQMATRICES[0] = &G11_P0;
      SQMATRICES[1] = &G12_P0;
      SQMATRICES[2] = &G13_P0;
      SQMATRICES[3] = &G21_P0;
      SQMATRICES[4] = &G22_P0;
      SQMATRICES[5] = &G23_P0;
      SQMATRICES[6] = &G31_P0;
      SQMATRICES[7] = &G32_P0;
      SQMATRICES[8] = &G33_P0;
      DirectSolver(SQMATRICES, 3, 3, grid_sol_p0, grid_rhs_p0);
      
      SQMATRICES[0] = &G11_P1;
      SQMATRICES[1] = &G12_P1;
      SQMATRICES[2] = &G13_P1;
      SQMATRICES[3] = &G21_P1;
      SQMATRICES[4] = &G22_P1;
      SQMATRICES[5] = &G23_P1;
      SQMATRICES[6] = &G31_P1;
      SQMATRICES[7] = &G32_P1;
      SQMATRICES[8] = &G33_P1;
      DirectSolver(SQMATRICES, 3, 3, grid_sol_p1, grid_rhs_p1);
      
//       for (int i=0;i<3*N_Grid_P0;++i)
// 	grid_sol_p0[i] = 1.;
//       
//       for (int i=0;i<3*N_Grid_P1;++i)
// 	grid_sol_p1[i] = 1.;

//       TFEFunction3D *func;
//       
//       func = g_p0->GetComponent(0);
//       func->Interpolate(ExactU1);
//       func = g_p0->GetComponent(1);
//       func->Interpolate(ExactU2);
//       func = g_p0->GetComponent(2);
//       func->Interpolate(ExactU3);
//       
//       func = g_p1->GetComponent(0);
//       func->Interpolate(ExactU1);
//       func = g_p1->GetComponent(1);
//       func->Interpolate(ExactU2);
//       func = g_p1->GetComponent(2);
//       func->Interpolate(ExactU3);
      
      MapGridVelo(grid_space_P0, g, GlobalCellNo_P0, grid_sol_p0);
      MapGridVelo(grid_space_P1, g, GlobalCellNo_P1, grid_sol_p1);
      Dscal(3*N_Grid, 1/tau, grid_velo);
      
      /// assemble nonlinear part
      FESP[0] = velocity_space;
      FESP[1] = grid_space;
      N_FESP = 2;
      
      SQMATRICES[0] = &A11;
      SQMATRICES[1] = &A22;
      SQMATRICES[2] = &A33;
      N_SQMAT = 3;
      N_MAT = 0,
      
      MatReset(N_SQMAT, SQMATRICES, N_MAT, NULL);
      
      Assemble3D (N_FESP, FESP,
		  N_SQMAT, SQMATRICES,
		  0, NULL,
		  0, NULL, NULL,
		  DiscreteFormNLGalerkin,
		  BOUNDCOND, BOUNDVALUE, &aux);
		  
      // add M blocks to A blocks 
      MatAdd2(&A11, &M11, tau);
      MatAdd2(&A22, &M22, tau);
      MatAdd2(&A33, &M33, tau);
      
      // calculate defect
      CoupledDefect(&A11, &A12, &A13, &A21, &A22, &A23, &A31, &A32, &A33,
		    &B1, &B2, &B3, &B1T, &B2T, &B3T, sol, rhs, res);
		  
      if ( TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE )
	IntoL20Vector3D(res+3*N_U, N_P, TDatabase::ParamDB->PRESSURE_SPACE);
      
      double norm = Dnorm(N_Unknowns, res);
      
      OutPut("Nonlinear iterationstep: " << iter);
      OutPut("  ; RES: " << norm << endl);
      
      if ( norm < TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE )
	break;
      
      // solve linear system
      Solver.SetMatrix(&A11, &A12, &A13, &A21, &A22, &A23, &A31, &A32, &A33,
		     &B1T, &B2T, &B3T, &B1, &B2, &B3);
      Solver.FactorizeSolve(sol, rhs);
//       DirectSolver(&A11, &A12, &A13, &A21, &A22, &A23, &A31, &A32, &A33,
// 		   &B1T, &B2T, &B3T, &B1, &B2, &B3, rhs, sol, 3);
      
    } // end nonlinear loop
    
    OutPut(endl);
    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    {
      OutPut("PROJECTING PRESSURE" << endl);
      IntoL20FEFunction3D(p->GetValues(), p->GetLength(), p->GetFESpace3D(),
			  TDatabase::ParamDB->VELOCITY_SPACE, TDatabase::ParamDB->PRESSURE_SPACE);
    }
    
    // Move grid
    Daxpy(3*N_Grid, tau, grid_velo, grid_pos);
    grid->DataToGrid();
    
    if ( timestep % TDatabase::TimeDB->STEPS_PER_IMAGE == 0 )
    {
      WriteSolution (VTK, img, &Output);
      WriteGrid("grid_p1", img, &Domain, Coll_P1);
      ++img;
    }

    Volumen_P1 = GetVolume(Coll_P1);
    OutPut("Volume phase 1: " << Volumen_P1 << endl);
    
    PrintMemory();    
  } // end time loop
  
  return 0;
}

void SplitGrid_TwoPhase(TDomain *Domain, TCollection *&Coll, TCollection *&Coll_P0,
			TCollection *&Coll_P1, TFEFunction3D *&phase_marker,
			int *&GlobalCellNo_P0, int *&GlobalCellNo_P1,
			int &N_SurfaceJoints_P0, int *&CellNumbers_P0, int *&JointNumbers_P0,
			int &N_SurfaceJoints_P1, int *&CellNumbers_P1, int *&JointNumbers_P1)
{
  int N_Cells, N_Cells_P0, N_Cells_P1, Phase, p0, p1;
  int N_DOF, *DOF, *GlobalNumbers, *BeginIndex, *JointDOF, N_JointDOF;
  TBaseCell *Cell, **Cells_P0, **Cells_P1, *Neib;
  TJoint *Joint;
  FE3D FeID;
  
  Coll = Domain->GetCollection(It_Finest, 0);
  
  N_Cells = Coll->GetN_Cells();
  
  N_Cells_P0 = N_Cells_P1 = 0;
  N_SurfaceJoints_P0 = 0;
  N_SurfaceJoints_P1 = 0;
  for (int i=0;i<N_Cells;++i)
  {
    Cell = Coll->GetCell(i);
    
    Phase = Cell->GetPhase_ID();
    
    switch (Phase)
    {
      case 0:
	++N_Cells_P0;
	break;	
      case 1:
	++N_Cells_P1;
	break;
	
      default:
	cerr << "invalid phase id: " << Phase << endl;
	exit(0);
    } 
    
    for (int j=0;j<4;++j)
    {
      Joint = Cell->GetJoint(j);
      Neib  = Joint->GetNeighbour(Cell);
      
      if ( Neib )
      {
	if ( Neib->GetPhase_ID() != Phase )
	{
	  switch (Phase)
	  {
	    case 0:    
	      ++N_SurfaceJoints_P0;
	      break;
	    case 1:
	      ++N_SurfaceJoints_P1;
	      break;
	  }
	}
      }
      
//       if ( Joint->GetType() == BoundaryFace )
//       {
// 	  switch (Phase)
// 	  {
// 	    case 0:    
// 	      ++N_SurfaceJoints_P0;
// 	      break;
// 	    case 1:
// 	      ++N_SurfaceJoints_P1;
// 	      break;
// 	  }
//       }
    }
  }
  
  OutPut("N_Cells phase 0: " << N_Cells_P0 << endl);
  OutPut("N_Cells phase 1: " << N_Cells_P1 << endl);
  
  Cells_P0 = new TBaseCell* [N_Cells_P0];
  Cells_P1 = new TBaseCell* [N_Cells_P1];
  GlobalCellNo_P0 = new int [N_Cells_P0];
  GlobalCellNo_P1 = new int [N_Cells_P1];
  CellNumbers_P0 = new int [N_SurfaceJoints_P0];
  CellNumbers_P1 = new int [N_SurfaceJoints_P1];
  JointNumbers_P0 = new int [N_SurfaceJoints_P0];
  JointNumbers_P1 = new int [N_SurfaceJoints_P1];
  
  p0 = p1 = 0;
  N_SurfaceJoints_P0 = 0;
  N_SurfaceJoints_P1 = 0;
  for (int i=0;i<N_Cells;++i)
  {
    Cell = Coll->GetCell(i);
    
    Phase = Cell->GetPhase_ID();
    
    switch (Phase)
    {
      case 0:
	Cells_P0[p0] = Cell;
	GlobalCellNo_P0[p0] = i;
	++p0;
	break;
	
      case 1:
	Cells_P1[p1] = Cell;
	GlobalCellNo_P1[p1] = i;
	++p1;
	break;
    }
    
    for (int j=0;j<4;++j)
    {
      Joint = Cell->GetJoint(j);
      Neib = Joint->GetNeighbour(Cell);
      
      if ( Neib )
      {
	if ( Neib->GetPhase_ID() != Phase )
	{
	  switch (Phase)
	  {
	    case 0:
	      CellNumbers_P0[N_SurfaceJoints_P0] = p0-1;
	      JointNumbers_P0[N_SurfaceJoints_P0] = j;
	      ++N_SurfaceJoints_P0;
	      break;
	    case 1:
	      CellNumbers_P1[N_SurfaceJoints_P1] = p1-1;
	      JointNumbers_P1[N_SurfaceJoints_P1] = j;
	      ++N_SurfaceJoints_P1;
	      break;
	  }
	}
      }
      
//       if ( Joint->GetType() == BoundaryFace )
//       {
// 	  switch (Phase)
// 	  {
// 	    case 0:
// 	      CellNumbers_P0[N_SurfaceJoints_P0] = p0-1;
// 	      JointNumbers_P0[N_SurfaceJoints_P0] = j;
// 	      ++N_SurfaceJoints_P0;
// 	      break;
// 	    case 1:
// 	      CellNumbers_P1[N_SurfaceJoints_P1] = p1-1;
// 	      JointNumbers_P1[N_SurfaceJoints_P1] = j;
// 	      ++N_SurfaceJoints_P1;
// 	      break;
// 	  }
//       }
    }
  } 
  
  Coll_P0 = new TCollection (N_Cells_P0, Cells_P0);
  Coll_P1 = new TCollection (N_Cells_P1, Cells_P1);
  
  TOutput3D Output(0, 0, 0, 0, Domain, Coll);
  TOutput3D Output_P0 (0, 0, 0, 0, Domain, Coll_P0);
  TOutput3D Output_P1 (0, 0, 0, 0, Domain, Coll_P1);
  
  Output.WriteVtk("grid.vtk");
  Output_P0.WriteVtk("grid.p0.vtk");
  Output_P1.WriteVtk("grid.p1.vtk");
  
  char name[] = "name";
  char desc[] = "desc";
  char marker[] = "marker";
  
  TFESpace3D *phase_space = new TFESpace3D (Coll, name, desc, GridBoundCond, 1);
  
  N_DOF = phase_space->GetN_DegreesOfFreedom();
  
  double *phase_space_values = new double [N_DOF];
  memset(phase_space_values, 0, N_DOF*sizeof(phase_space_values[0]));
  phase_marker = new TFEFunction3D (phase_space, marker, desc, phase_space_values, N_DOF);
  
  GlobalNumbers = phase_space->GetGlobalNumbers();
  BeginIndex    = phase_space->GetBeginIndex();
  
  for (int i=0;i<N_Cells;++i)
  {
    Cell = Coll->GetCell(i);
    
    FeID = phase_space->GetFE3D(i, Cell);
    
    DOF = GlobalNumbers + BeginIndex[i];
    
    Phase = Cell->GetPhase_ID();
    
    N_DOF = TFEDatabase3D::GetFEDesc3DFromFE3D(FeID)->GetN_DOF();
    
    for (int j=0;j<N_DOF;++j)
    {
      phase_space_values[DOF[j]] = 2*Phase - 1;
    }
  }
  
  for (int i=0;i<N_Cells;++i)
  {
    Cell = Coll->GetCell(i);
    
    FeID = phase_space->GetFE3D(i, Cell);
    
    DOF = GlobalNumbers + BeginIndex[i];
    
    Phase = Cell->GetPhase_ID();
    
    for (int j=0;j<4;++j)
    {
      Joint = Cell->GetJoint(j);
      
      Neib = Joint->GetNeighbour(Cell);
      
      p0 = Phase;
      
      if (Neib) p0 = Neib->GetPhase_ID();
      
      if ( Phase != p0 )
      {
	N_JointDOF = TFEDatabase3D::GetFEDesc3DFromFE3D(FeID)->GetN_JointDOF();
	JointDOF   = TFEDatabase3D::GetFEDesc3DFromFE3D(FeID)->GetJointDOF(j);
	
	for (int k=0;k<N_JointDOF;++k)
	{
	  phase_space_values[DOF[JointDOF[k]]] = 0.;
	}
      }
    }
  }
  
  {
    TBaseCell **SurfCells = new TBaseCell* [N_SurfaceJoints_P0];
    int N_SurfCells = N_SurfaceJoints_P0;
    for (int i=0;i<N_SurfaceJoints_P0;++i)
    {
      SurfCells[i] = Coll_P0->GetCell(CellNumbers_P0[i]);
    }
    
    TCollection *SurfColl = new TCollection (N_SurfCells, SurfCells);
    TOutput3D Output_Surf(0, 0, 0, 0, Domain, SurfColl);
    Output_Surf.WriteVtk("surf.p0.vtk");
  }
  
  {
    TBaseCell **SurfCells = new TBaseCell* [N_SurfaceJoints_P1];
    int N_SurfCells = N_SurfaceJoints_P1;
    for (int i=0;i<N_SurfaceJoints_P1;++i)
    {
      SurfCells[i] = Coll_P1->GetCell(CellNumbers_P1[i]);
    }
    
    TCollection *SurfColl = new TCollection (N_SurfCells, SurfCells);
    TOutput3D Output_Surf(0, 0, 0, 0, Domain, SurfColl);
    Output_Surf.WriteVtk("surf.p1.vtk");
  }
}

void CheckMesh(TFESpace3D *fesp)
{
  int N_Cells, N_DOF, found;
  TBaseCell *Cell;
  TCollection *Coll;
  int *BeginIndex, *GlobalNumbers, *DOF;
  int ActiveBound;
  FE3D feid;
  TFE3D *fe;
  
  Coll = fesp->GetCollection();
  N_Cells = Coll->GetN_Cells();
  ActiveBound = fesp->GetActiveBound();
  
  GlobalNumbers = fesp->GetGlobalNumbers();
  BeginIndex    = fesp->GetBeginIndex();
  
  for (int i=0;i<N_Cells;++i)
  {
    Cell = Coll->GetCell(i);

    DOF = GlobalNumbers + BeginIndex[i];
    
    feid = fesp->GetFE3D(i, Cell);
    fe   = TFEDatabase3D::GetFE3D(feid);
    N_DOF = fe->GetN_DOF();
    
    found = 0;
    for (int j=0;j<N_DOF;++j)
    {
      if ( DOF[j] < ActiveBound )
      {
	found = 1;
	break;
      }
    }
    
    if ( found == 0 )
    {
      cerr << "Mesh check failed! Exit" << endl;
      exit(0);
    }
  }
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

void MatReset(int N_SQMAT, TSquareMatrix3D **SQMATRICES, int N_MAT, TMatrix3D **MATRICES)
{
  for (int i=0;i<N_SQMAT;++i)
    SQMATRICES[i]->Reset();
  
  for (int i=0;i<N_MAT;++i)
    MATRICES[i]->Reset();
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

void WriteGrid(const char *basename, int number, TDomain *Domain, TCollection *Phase)
{
  std::ostringstream os;
  
  os << basename << ".";
  os.width(5);
  os.fill('0');
  os << number << ".vtk" << ends;
  
  TOutput3D Output(0,0,0,0, Domain, Phase);
  Output.WriteVtk(os.str().c_str());
}

void Remesh(TCollection *Coll, TDomain *Domain)
{
  int N_Cells, counter, counter2, MaxLen, Phase, N_Points;
  int *CellNumbers, *JointNumbers, CellNr, JointNr, N_SurfaceJoints;
  int *Facets, index;
  const int *TmpFV, *TmpLen;
  double *Points, x, y, z, P0x, P0y, P0z, P1x, P1y, P1z, Regions[10];
  JointType jtype;
  TBaseCell *Cell, *Neib, *Cell_inner, *Cell_outer;
  TJoint *Joint;
  TVertex *Vertex;
  
  N_Cells = Coll->GetN_Cells();
  
  /// reset and count surface joints
  N_SurfaceJoints = 0;
  index = -1;
  for (int i=0;i<N_Cells;++i)
  {
    Cell = Coll->GetCell(i);
    
    Phase = Cell->GetPhase_ID();
    if ( Phase == 1 )
    {
      if ( index  == -1 )
      {
	Cell_inner = Cell;
	index = 6;
      }
      continue;
    }
    
    Cell_outer = Cell;
    
    for (int j=0;j<4;++j)
    {
      Vertex = Cell->GetVertex(j);
      Vertex->SetClipBoard(-1);
      
      Joint = Cell->GetJoint(j);
      
      if ( jtype == BoundaryFace || IsoBoundFace )
	++N_SurfaceJoints;
      
      Neib = Joint->GetNeighbour(Cell);
      
      if ( Neib && Neib->GetPhase_ID() != Phase )
	++N_SurfaceJoints;
    }
  } // end for i
  
  CellNumbers = new int [N_SurfaceJoints];
  JointNumbers = new int [N_SurfaceJoints];
  
  /// count vertices on surfaces and surface cells
  counter = 0;
  counter2 = 0;
  for (int i=0;i<N_Cells;++i)
  {
    Cell = Coll->GetCell(i);
    
    Cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
    
    Phase = Cell->GetPhase_ID();
    if ( Phase == 1 ) continue;
    
    for (int j=0;j<4;++j)
    {
      Joint = Cell->GetJoint(j);
      jtype = Joint->GetType();
      
      // boundary
      if ( jtype == BoundaryFace || IsoBoundFace )
      {
	for (int k=0;k<TmpLen[j];++k)
	{
	  Vertex = Cell->GetVertex(TmpFV[j*MaxLen+k]);
	  if ( Vertex->GetClipBoard() == -1 )
	  {
	    Vertex->SetClipBoard(counter++);
	  }
	}
	
	CellNumbers[counter2] = i;
	JointNumbers[counter2] = j;
	++counter2;
      }
      
      //interface
      Neib = Joint->GetNeighbour(Cell);
      if ( Neib )
      {
	if ( Neib->GetPhase_ID() != Phase )
	{
	  for (int k=0;k<TmpLen[j];++k)
	  {
	    Vertex = Cell->GetVertex(TmpFV[j*MaxLen+k]);
	    if ( Vertex->GetClipBoard() == -1)
	    {
	      Vertex->SetClipBoard(counter++);
	    }
	  }
	  
	  CellNumbers[counter2] = i;
	  JointNumbers[counter2] = j;
	  ++counter2;
	}
      }
    } // end for j
  } // end for i
  
  N_Points = counter;
  Points = new double [3*N_Points];
  Facets = new int [3*N_SurfaceJoints];
  
  for (int i=0;i<N_SurfaceJoints;++i)
  {
    CellNr = CellNumbers[i];
    JointNr = JointNumbers[i];
    
    Cell = Coll->GetCell(CellNr);
    Joint = Cell->GetJoint(JointNr);
    
    Cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
    
    for (int j=0;j<TmpLen[JointNr];++j)
    {
      Vertex = Cell->GetVertex(TmpFV[JointNr*MaxLen+j]);
      counter = Vertex->GetClipBoard();
      Vertex->GetCoords(x, y, z);
      
      Points[3*counter+0] = x;
      Points[3*counter+1] = y;
      Points[3*counter+2] = z;
      
      Facets[3*i+j] = counter;
    }
  }
  
  P0x = P0y = P0z = 0;
  P1x = P1y = P1z = 0;
  for (int i=0;i<4;++i)
  {
    Cell_outer->GetVertex(i)->GetCoords(x, y, z);
    P0x += 0.25*x;
    P0y += 0.25*y;
    P0z += 0.25*z;
    
    Cell_inner->GetVertex(i)->GetCoords(x, y, z);
    P1x += 0.25*x;
    P1y += 0.25*y;
    P1z += 0.25*z;
  }
  
  Regions[0] = P0x;
  Regions[1] = P0y;
  Regions[2] = P0z;
  Regions[3] = 0;
  Regions[4] = 0;
  Regions[5] = P1x;
  Regions[6] = P1y;
  Regions[7] = P1z;
  Regions[8] = 1;
  Regions[9] = 0;
  
  TTetGenMeshLoader::TDummyDomain DummyDomain;
  TTetGenMeshLoader MeshLoader;
  
  MeshLoader.Generate(N_Points, Points, N_SurfaceJoints, Facets, 2, Regions, &DummyDomain);
  
  delete [] CellNumbers;
  delete [] JointNumbers;
  delete [] Points;
}

