#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <Matrix3D.h>
#include <SquareMatrix3D.h>
#include <LinAlg.h>
#include <DirectSolver.h>
#include <Output3D.h>
#include <DiscreteForm3D.h>
#include <MainUtilities.h>
#include <MooNMD_Io.h>
#include <AuxParam3D.h>
#include <TNSE3D_ParamRout.h>
#include <Assemble3D.h>
#include <PardisoSolver.h>
#include <MovingTNSE3D.h>
#include <FreeSurface3D.h>
#include <TetraAffin.h>
#include <BoundFace.h>
#include <IsoBoundFace.h>
#include <TetraIsoparametric.h>

#include <stdlib.h>
#include <string.h>
#include <sstream>

// #define UMFPACK
// #define EXACTCURV_SPHERE

// #include "../Examples/TNSE_3D/AnsatzLinConst.h"
// #include "../Examples/TNSE_3D/Benchmark3D.h"
#include "../Examples/TNSE_3D/FreeSurf3D.h"

void PrintMemory();
void WriteGrid(const char *name, TDomain *Domain);
void WriteGrid(const char *name, TDomain *Domain, TCollection *Coll);
void WriteSolution(const char *basename, int number, TOutput3D *Output);
void WriteSurface(TDomain *Domain, TCollection *Coll, int N_BoundFaces, int *CellNumbers);

void PrintRhs(TCollection *Coll, int N_BoundFaces, int *CellNumbers, int *JointNumbers, TFESpace3D *fesp,
	      double *rhs, double *rhs_ref);

void MakeIsoparametric_Sphere(TCollection *Coll, int N_BoundFaces, int *CellNumbers,
			      int *JointNumbers, TFESpace3D *fespace,
			      TVertex **&AllVertices, int *&VertexCell, int &N_Vertices);
			      
void MoveIsoVertices(TFEVectFunct3D *velocity, int N_BoundFaces, int *CellNumbers, int *JointNumbers,
		     double tau);
			      
void CheckIso(TCollection *Coll, int N_BoundFaces, int *CellNumbers, int *JointNumbers, int img);
	      
void CheckMesh(TFESpace3D *fesp);

typedef struct
{
  TSquareMatrix3D *A11;
  TSquareMatrix3D *A12;
  TSquareMatrix3D *A13;
  TSquareMatrix3D *A21;
  TSquareMatrix3D *A22;
  TSquareMatrix3D *A23;
  TSquareMatrix3D *A31;
  TSquareMatrix3D *A32;
  TSquareMatrix3D *A33;
  
  TDirectSolver *Solver;
  
  void Init (TSquareStructure3D *structure);
  
  void PrepareAssemble(int &N_FESP, TFESpace3D **FESP,
		       int &N_SQMAT, TSquareMatrix3D **SQMATRICES,
		       int &N_MAT, TMatrix3D **MATRICES);
  
  void GetMatrices(TSquareMatrix3D **SQMATRICES);
  void Solve(double *sol, double *rhs);
		       
} MATRICES_GRID;

typedef struct
{
  TSquareMatrix3D *sqmatrixM;
  TSquareMatrix3D *sqmatrixA;
  TSquareMatrix3D *sqmatrixA11;
  TSquareMatrix3D *sqmatrixA12;
  TSquareMatrix3D *sqmatrixA13;
  TSquareMatrix3D *sqmatrixA21;
  TSquareMatrix3D *sqmatrixA22;
  TSquareMatrix3D *sqmatrixA23;
  TSquareMatrix3D *sqmatrixA31;
  TSquareMatrix3D *sqmatrixA32;
  TSquareMatrix3D *sqmatrixA33;
  TSquareMatrix3D *sqmatrixM11;
  TSquareMatrix3D *sqmatrixM22;
  TSquareMatrix3D *sqmatrixM33;
  TMatrix3D       *matrixB1;
  TMatrix3D       *matrixB2;
  TMatrix3D       *matrixB3;
  TMatrix3D       *matrixB1T;
  TMatrix3D       *matrixB2T;
  TMatrix3D       *matrixB3T;
  
  TDirectSolver *Solver;
  
  int NSTYPE;
  
  void Init(TSquareStructure3D *sqstructureA, TStructure3D *structureB,
	    TStructure3D *structureBT);
	    
  int PrepareAssemble(int &N_FESP, TFESpace3D **FESP,
		       int &N_SQMAT, TSquareMatrix3D **SQMATRICES,
		       int &N_MAT, TMatrix3D **MATRICES);
		       
  int PrepareAssembleNL(int &N_FESP, TFESpace3D **FESP,
		         int &N_SQMAT, TSquareMatrix3D **SQMATRICES,
		         int &N_MAT, TMatrix3D **MATRICES);
			 
  double GetDefect(int N_, double *sol, double *rhs);
		       
  void Solve(double *sol, double *rhs);
  
  void BuildSystem();
  
  void AddToRhs(double tau, double *sol, double *rhs);
  void AddMToA(double tau);
  void AddMToA_diag(double tau);
  void Scale(double factor);
		       
  void Print();
  
} MATRICES_NSTYPE;

int main(int argc, char **argv)
{
  TDomain Domain;
  TDatabase Database;
  TFEDatabase3D FEDatabase;
  TCollection *Coll;
  
  TFESpace3D *velocity_space, *pressure_space;
  
  TDiscreteForm3D *DiscreteFormGalerkin, *DiscreteFormNLGalerkin;
  TDiscreteForm3D *DiscreteFormRHS, *DiscreteFormGrid;
  
  MATRICES_NSTYPE Matrices;
  MATRICES_GRID   Matrices_Grid;
  
  double Volumen;
  
  int ret, NSTYPE, pressure_space_code, index;
  int img=0;
  
  BoundCondFunct3D      *BOUNDCOND[4], *GRIDBOUNDCOND[4];
  BoundValueFunct3D     *BOUNDVALUE[4], *GRIDBOUNDVALUE[4];
  
  char *GEO, *BND, *SMESH, *VTK;
  
  char Name[]  = "name";
  char UName[] = "u";
  char PName[] = "p";
  char GName[] = "g";
  char EUName[] = "exactU";
  char EPName[] = "exactP";
  char Desc[]  = "desc";
  
  if (argc > 1)
    ret = Domain.ReadParam(argv[1]);

  else 
  {
    cerr << "no data file given !" << endl;
    return 0;
  }
  
  ExampleFile();
  
  GEO = TDatabase::ParamDB->GEOFILE;
  BND = TDatabase::ParamDB->BNDFILE;
  SMESH = TDatabase::ParamDB->SMESHFILE;
  NSTYPE = TDatabase::ParamDB->NSTYPE;
  VTK = TDatabase::ParamDB->VTKBASENAME;
  
  //
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
  
  // create mesh
//   Domain.Init(BND, GEO);
  Domain.Tetgen(SMESH);

  WriteGrid("grid_coarse.vtk", &Domain);
  
  for (int i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;++i)
    Domain.RegRefineAll();
  
  WriteGrid("grid_fine.vtk", &Domain);

  // discrete forms
  InitializeDiscreteForms(DiscreteFormGalerkin, DiscreteFormNLGalerkin,
			  DiscreteFormRHS, LinCoeffs, NSTYPE);
			  
  DiscreteFormGrid = new TDiscreteForm3D (Name, Desc, ESGridN_Terms, ESGridDerivatives,
					  ESGridSpaceNumbers, ESGridN_Matrices,
					  ESGridN_Rhs, ESGridRowSpace, ESGridColumnSpace,
					  ESGridRhsSpace, ESGridAssemble, GridCoeffs, NULL);
			  
  // collection 
  Coll = Domain.GetCollection(It_Finest, 0);
  
  // init fespaces for velocity and pressure
  GetVelocityAndPressureSpace3D (Coll, BoundCondition, velocity_space,
				 pressure_space, &pressure_space_code,
				 TDatabase::ParamDB->VELOCITY_SPACE,
				 TDatabase::ParamDB->PRESSURE_SPACE);
			
//   TDatabase::ParamDB->PRESSURE_SPACE = pressure_space_code;
  
  // find free surface
  int N_SurfaceJoints, *CellNumbers, *JointNumbers, *MoveVertexCell=NULL, N_MoveVertices;
  TVertex **MoveVertices=NULL;
  FindFreeSurfaceFromJointType(Coll, BoundaryFace, N_SurfaceJoints,
		  CellNumbers, JointNumbers);
	
  if ( TDatabase::ParamDB->USE_ISOPARAMETRIC )
  {
    MakeIsoparametric_Sphere(Coll, N_SurfaceJoints, CellNumbers, JointNumbers, velocity_space,
			     MoveVertices, MoveVertexCell, N_MoveVertices);
  }
// exit(0);  
  OutPut("number of surface joints: " << N_SurfaceJoints << endl);
  WriteSurface(&Domain, Coll, N_SurfaceJoints, CellNumbers);

  TFESpace3D *grid_space = new TFESpace3D (Coll, Name, Desc, GridBoundCond, 1);
  
  CheckMesh(velocity_space);
  
  // count DOF's
  int N_U = velocity_space->GetN_DegreesOfFreedom();
  int N_P = pressure_space->GetN_DegreesOfFreedom();
  int N_Grid = grid_space->GetN_DegreesOfFreedom();
  int N_Unknowns = 3*N_U + N_P;
  int InnerBound = velocity_space->GetInnerBound();
  
  // arrays for solution and right hand side
  double *sol = new double [N_Unknowns];
  double *rhs = new double [N_Unknowns];
  double *exact = new double [N_Unknowns];
  double *rhs_test = new double [N_Unknowns];
//   double *rhs_tmp = new double [N_U];
  double *grid_pos  = new double [3*N_Grid];
  double *grid_velo = new double [3*N_Grid]; 
  double *grid_rhs  = new double [3*N_Grid];
  memset(sol, 0, N_Unknowns*sizeof(sol[0]));
  memset(rhs, 0, N_Unknowns*sizeof(rhs[0]));
  memset(exact, 0, N_Unknowns*sizeof(exact[0]));
  memset(grid_velo, 0, sizeof(grid_velo[0])*3*N_Grid);
  memset(grid_rhs, 0, sizeof(grid_rhs[0])*3*N_Grid);
  memset(grid_pos, 0, sizeof(grid_pos[0])*3*N_Grid);
  
  // fe function for velocity and pressure
  TFEVectFunct3D *u = new TFEVectFunct3D (velocity_space, UName, Desc, sol, N_U, 3);
  TFEFunction3D *ux = u->GetComponent(0);
  TFEFunction3D *uy = u->GetComponent(1);
  TFEFunction3D *uz = u->GetComponent(2);
  
  TFEFunction3D *p = new TFEFunction3D (pressure_space, PName, Desc, sol+3*N_U, N_P);
  
  TFEVectFunct3D *g = new TFEVectFunct3D (grid_space, GName, Desc, grid_velo, N_Grid, 3);
  TFEFunction3D *gx = g->GetComponent(0);
  TFEFunction3D *gy = g->GetComponent(1);
  TFEFunction3D *gz = g->GetComponent(2);
  
  TFEVectFunct3D *grid = new TFEVectFunct3D (grid_space, GName, Desc, grid_pos, N_Grid, 3);
  
  TFEVectFunct3D *fExactU = new TFEVectFunct3D (velocity_space, EUName, Desc, exact, N_U, 3);
  TFEFunction3D *fExactU1 = fExactU->GetComponent(0);
  TFEFunction3D *fExactU2 = fExactU->GetComponent(1);
  TFEFunction3D *fExactU3 = fExactU->GetComponent(2);
  TFEFunction3D *fExactP  = new TFEFunction3D (pressure_space, EPName, Desc, exact+3*N_U, N_P);
  
  
  TFEVectFunct3D *test = new TFEVectFunct3D (velocity_space, Name, Desc, rhs, N_U, 3);
  
  // output
  TOutput3D Output(2, 1, 1, 0, &Domain);
  Output.AddFEVectFunct(u);
  Output.AddFEVectFunct(g);
  Output.AddFEFunction(p);
  Output.AddFEVectFunct(test);
  Output.AddFEVectFunct(fExactU);
  Output.AddFEFunction(fExactP);
  
  // matrices structures
  TSquareStructure3D sqstructureA (velocity_space);
  TSquareStructure3D sqstructureGrid (grid_space);
  TStructure3D structureB  (pressure_space, velocity_space);
  TStructure3D structureBT (velocity_space, pressure_space);
  sqstructureA.Sort();
  structureB.Sort();
  structureBT.Sort();
  
  // init matrices
  Matrices.Init(&sqstructureA, &structureB, &structureBT);
  Matrices_Grid.Init(&sqstructureGrid);
  
  
  // prepare assemble
  int N_RHS, N_SQMAT, N_MAT, N_FESP;
  double *RHS[3], *RHS_grid[3];
  TFESpace3D *FESP[4], *FESPRHS[4];
  TFEFunction3D *FEFCT[6];
  TSquareMatrix3D *SQMATRICES[9], *FREESURF[3];
  TMatrix3D *MATRICES[6];
  
  RHS[0] = rhs        ;
  RHS[1] = rhs +   N_U;
  RHS[2] = rhs + 2*N_U;
  RHS_grid[0] = grid_rhs           ;
  RHS_grid[1] = grid_rhs +   N_Grid;
  RHS_grid[2] = grid_rhs + 2*N_Grid;
  N_RHS = 3;
  
  FESPRHS[0] = velocity_space;
  FESPRHS[1] = velocity_space;
  FESPRHS[2] = velocity_space;
			   
  // aux object
  FEFCT[0] = ux; FEFCT[1] = uy; FEFCT[2] = uz;
  FEFCT[3] = gx; FEFCT[4] = gy; FEFCT[5] = gz;
  TAuxParam3D aux (MovingTimeNSN_FESpacesVelo, MovingTimeNSN_FctVelo,
		   MovingTimeNSN_ParamFctVelo, MovingTimeNSN_FEValuesVelo,
		   FESP, FEFCT, MovingTimeNSFctVelo, MovingTimeNSFEFctIndexVelo,
		   MovingTimeNSFEMultiIndexVelo, MovingTimeNSN_ParamsVelo,
		   MovingTimeNSBeginParamVelo);

//   TAuxParam3D aux (TimeNSN_FESpacesVelo, TimeNSN_FctVelo,
// 		   TimeNSN_ParamFctVelo, TimeNSN_FEValuesVelo,
// 		   FESP, FEFCT, TimeNSFctVelo, TimeNSFEFctIndexVelo,
// 		   TimeNSFEMultiIndexVelo, TimeNSN_ParamsVelo,
// 		   TimeNSBeginParamVelo);
		   
  TAuxParam3D aux_dummy (0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
	      
  // time loop (back euler)
  double tau;
  tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  
  ux->Interpolate(InitialU1);
  uy->Interpolate(InitialU2);
  uz->Interpolate(InitialU3);
  p->Interpolate(InitialP);
  
  OutPut ("Velocity DOF: " << 3*N_U << endl);
  OutPut ("Pressure DOF: " << N_P << endl);
  OutPut ("Grid DOF: " << 3*N_Grid << endl);
  OutPut ("Number of unknowns: " << N_Unknowns);
  OutPut (" + " << 3*N_Grid << endl);
  
  OutPut(endl);
  OutPut("RE_NR: " << TDatabase::ParamDB->RE_NR << endl);
  OutPut("WE_NR: " << TDatabase::ParamDB->WB_NR << endl);
  
  OutPut(endl);
  OutPut("STARTTIME: " << TDatabase::TimeDB->CURRENTTIME << endl);
  OutPut("ENDTIME: " << TDatabase::TimeDB->ENDTIME << endl);
  OutPut("TIMESTEPLENGTH: " << tau << endl);
  Volumen = GetVolume(Coll);
  OutPut("Volume: " << Volumen << endl << endl);
  
  fExactU1->Interpolate(ExactU1);
  fExactU2->Interpolate(ExactU2);
  fExactU3->Interpolate(ExactU3);
  fExactP->Interpolate(ExactP);
  WriteSolution(VTK, img, &Output);
  if ( TDatabase::ParamDB->USE_ISOPARAMETRIC )
  {
     CheckIso(Coll, N_SurfaceJoints, CellNumbers, JointNumbers, img);
  }
  ++img;
  PrintMemory();
  
  // init grid
  grid->GridToData();
  
  for(int timestep=0;TDatabase::TimeDB->CURRENTTIME <= TDatabase::TimeDB->ENDTIME;++timestep)
  {
    TDatabase::TimeDB->CURRENTTIME += tau;
   
    OutPut (endl << "CURRENTTIME: " << TDatabase::TimeDB->CURRENTTIME << endl);
    
    // assemble complete system with rhs
    memset(rhs, 0, N_Unknowns*sizeof(rhs[0]));
    memset(rhs_test, 0, N_Unknowns*sizeof(rhs[0]));
    
    index = Matrices.PrepareAssemble (N_FESP, FESP,
				      N_SQMAT, SQMATRICES,
				      N_MAT, MATRICES);
				      
    FESP[index] = grid_space;
    ++N_FESP;
			      
    Assemble3D (N_FESP, FESP,
		N_SQMAT, SQMATRICES,
		N_MAT, MATRICES,
		N_RHS, RHS, FESPRHS,
		DiscreteFormGalerkin,
		BOUNDCOND, BOUNDVALUE, &aux);
		
    // assemble grid
    Matrices_Grid.PrepareAssemble(N_FESP, FESP,
				N_SQMAT, SQMATRICES,
				N_MAT, MATRICES);
				
    Assemble3D (N_FESP, FESP,
	      N_SQMAT, SQMATRICES,
	      N_MAT, MATRICES,
	      0, NULL, NULL,
	      DiscreteFormGrid,
	      GRIDBOUNDCOND, GRIDBOUNDVALUE, &aux_dummy);
	      
    // scale B's and BT's
    Matrices.Scale(tau);
    
//     #ifndef EXACTCURV_SPHERE
//     FreeSurfInt(Coll, N_SurfaceJoints, CellNumbers, JointNumbers, &aux, tau,
// 		  FREESURF, RHS[0], RHS[1], RHS[2]);
//     #else
//     FreeSurfInt_Sphere(velocity_space, tau, RHS[0], RHS[1], RHS[2]);
//     #endif
    
    Matrices.AddToRhs(tau, sol, rhs);
    
    switch (NSTYPE)
    {
      case 4:
	FREESURF[0] = Matrices.sqmatrixM11;
	FREESURF[1] = Matrices.sqmatrixM22;
	FREESURF[2] = Matrices.sqmatrixM33;
	break;
      case 2:
	FREESURF[0] = Matrices.sqmatrixM;
    }
    
    #ifndef EXACTCURV_SPHERE
    FreeSurfInt(Coll, N_SurfaceJoints, CellNumbers, JointNumbers, &aux, tau,
		  FREESURF, RHS[0], RHS[1], RHS[2]);

//     FreeSurfInt_new(Coll, N_SurfaceJoints, CellNumbers, JointNumbers, tau,
// 		    FREESURF, RHS[0], RHS[1], RHS[2]);
		  
    FreeSurfInt_Sphere(velocity_space, tau, rhs_test, rhs_test+N_U, rhs_test+2*N_U);
    
    PrintRhs(Coll, N_SurfaceJoints, CellNumbers, JointNumbers, velocity_space, rhs, rhs_test);
    #else
    FreeSurfInt_Sphere(velocity_space, tau, RHS[0], RHS[1], RHS[2]);
    FreeSurfInt(Coll, N_SurfaceJoints, CellNumbers, JointNumbers, &aux, tau,
		  FREESURF, rhs_test, rhs_test+N_U, rhs_test+2*N_U);
		  
    PrintRhs(Coll, N_SurfaceJoints, CellNumbers, JointNumbers, velocity_space, rhs, rhs_test);
    #endif
    
//     OutPut(Dnorm(3*N_U, rhs) << endl);
//     for (int ii=0;ii<3*N_U;++ii)
//     {
//       OutPut(rhs[ii] << " ");
//     }
    
    Matrices.AddMToA(tau);
    
    #ifndef UMFPACK
    if ( timestep == 0 || 1)
    {
//     Matrices.Solve(sol, rhs);
      Matrices.BuildSystem();
      Matrices.Solver->Analyse();
      Matrices.Solver->Factorize();
      Matrices.Solver->Solve(sol, rhs);
    }
    else 
      Matrices.Solve(sol, rhs);
    #endif
    
    #ifdef UMFPACK
    DirectSolver(Matrices.sqmatrixA, Matrices.matrixB1T, Matrices.matrixB2T,
		 Matrices.matrixB3T, Matrices.matrixB1, Matrices.matrixB2,
		 Matrices.matrixB3, rhs, sol);
    #endif
    
//     if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
//       IntoL20FEFunction3D(p->GetValues(), p->GetLength(), p->GetFESpace3D());
	       
    for (int iter=0;iter<TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;++iter)
    {
      double res;
      
      // solve grid equation
      memset(grid_rhs, 0, 3*N_Grid*sizeof(grid_rhs[0]));
      SetGridRhs(u, grid_space, N_SurfaceJoints, CellNumbers, JointNumbers, tau,
	       RHS_grid[0], RHS_grid[1], RHS_grid[2]);
    
      Matrices_Grid.Solve(grid_velo, grid_rhs);
    
      // scale to get grid velocity
      Dscal(3*N_Grid, 1/tau, grid_velo);
      
      // assemble nonlinear part
      index  = Matrices.PrepareAssembleNL(N_FESP, FESP,
					  N_SQMAT, SQMATRICES,
					  N_MAT, MATRICES);
					  
      FESP[index] = grid_space;
      ++N_FESP;
			
      Assemble3D (N_FESP, FESP,
		  N_SQMAT, SQMATRICES,
		  N_MAT, MATRICES,
		  0, NULL, NULL,
		  DiscreteFormNLGalerkin,
		  BOUNDCOND, BOUNDVALUE, &aux);
		  
      Matrices.AddMToA_diag(tau);
      
      res = Matrices.GetDefect(N_Unknowns, sol, rhs);
      OutPut ("Nonlinear iterationstep: " << iter);
      OutPut ("  ; RES: " << res << endl);
      
      if ( res < TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE )
	break;
      
      #ifndef UMFPACK
      Matrices.Solve(sol, rhs); 
      #endif
      #ifdef UMFPACK
      DirectSolver(Matrices.sqmatrixA, Matrices.matrixB1T, Matrices.matrixB2T,
		 Matrices.matrixB3T, Matrices.matrixB1, Matrices.matrixB2,
		 Matrices.matrixB3, rhs, sol);
      #endif   
      
//       if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
// 	IntoL20FEFunction3D(p->GetValues(), p->GetLength(), p->GetFESpace3D());
    } // end nonlin loop
    
    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    {
      OutPut("PRJECTING PRESSURE" << endl);
      IntoL20FEFunction3D(p->GetValues(), p->GetLength(), p->GetFESpace3D(),
			  TDatabase::ParamDB->VELOCITY_SPACE, TDatabase::ParamDB->PRESSURE_SPACE);
    }
    
    // move grid
    Daxpy(3*N_Grid, tau, grid_velo, grid_pos);
    grid->DataToGrid();
    if ( TDatabase::ParamDB->USE_ISOPARAMETRIC )
    {
      MoveIsoVertices(u, N_SurfaceJoints, CellNumbers, JointNumbers, tau);
    }
    
    // scale to get grid shift
    Dscal(3*N_Grid, tau, grid_velo);
    
    OutPut(endl);
//     OutPut(Dnorm(3*N_U, rhs) << endl);
//     for (int ii=0;ii<3*N_U;++ii)
//     {
//       OutPut(rhs[ii] << " ");
//     }

    fExactU1->Interpolate(ExactU1);
    fExactU2->Interpolate(ExactU2);
    fExactU3->Interpolate(ExactU3);
    fExactP->Interpolate(ExactP);
    WriteSolution(VTK, img, &Output);
    if ( TDatabase::ParamDB->USE_ISOPARAMETRIC && 1 == (int) TDatabase::ParamDB->P4 )
    {
      CheckIso(Coll, N_SurfaceJoints, CellNumbers, JointNumbers, img);
    }
    ++img;
    PrintMemory();
    
    // scale to get grid velocity
    Dscal(3*N_Grid, 1/tau, grid_velo);
    
    Volumen = GetVolume(Coll);
    OutPut("Volume: " << Volumen << endl);
    
  } // end time loop
  
  OutPut(endl << "ENDE" << endl);
  
  return 0;
}

void MATRICES_NSTYPE::Init(TSquareStructure3D *sqstructureA, TStructure3D *structureB,
			  TStructure3D *structureBT)
{
  NSTYPE = TDatabase::ParamDB->NSTYPE;
  
  sqmatrixM = NULL;
  sqmatrixA = NULL;
  matrixB1 = NULL;
  matrixB2 = NULL;
  matrixB3 = NULL;
  matrixB1T = NULL;
  matrixB2T = NULL;
  matrixB3T = NULL;
  
  sqmatrixA11 = NULL;
  sqmatrixA12 = NULL;
  sqmatrixA13 = NULL;
  sqmatrixA21 = NULL;
  sqmatrixA22 = NULL;
  sqmatrixA23 = NULL;
  sqmatrixA31 = NULL;
  sqmatrixA32 = NULL;
  sqmatrixA33 = NULL;
  sqmatrixM11 = NULL;
  sqmatrixM22 = NULL;
  sqmatrixM33 = NULL;
  
  switch (NSTYPE)
  {
    case 1:
	sqmatrixA = new TSquareMatrix3D (sqstructureA);
	sqmatrixM = new TSquareMatrix3D (sqstructureA);
	matrixB1  = new TMatrix3D (structureB);
	matrixB2  = new TMatrix3D (structureB);
	matrixB3  = new TMatrix3D (structureB);
	
	break;
      
    case 2:
	sqmatrixA = new TSquareMatrix3D (sqstructureA);
	sqmatrixM = new TSquareMatrix3D (sqstructureA);
	matrixB1  = new TMatrix3D (structureB);
	matrixB2  = new TMatrix3D (structureB);
	matrixB3  = new TMatrix3D (structureB);
	matrixB1T  = new TMatrix3D (structureBT);
	matrixB2T  = new TMatrix3D (structureBT);
	matrixB3T  = new TMatrix3D (structureBT);
	
	break;
	
    case 4:
	sqmatrixA11 = new TSquareMatrix3D (sqstructureA);
	sqmatrixA12 = new TSquareMatrix3D (sqstructureA);
	sqmatrixA13 = new TSquareMatrix3D (sqstructureA);
	sqmatrixA21 = new TSquareMatrix3D (sqstructureA);
	sqmatrixA22 = new TSquareMatrix3D (sqstructureA);
	sqmatrixA23 = new TSquareMatrix3D (sqstructureA);
	sqmatrixA31 = new TSquareMatrix3D (sqstructureA);
	sqmatrixA32 = new TSquareMatrix3D (sqstructureA);
	sqmatrixA33 = new TSquareMatrix3D (sqstructureA);
	sqmatrixM11 = new TSquareMatrix3D (sqstructureA);
	sqmatrixM22 = new TSquareMatrix3D (sqstructureA);
	sqmatrixM33 = new TSquareMatrix3D (sqstructureA);
	matrixB1  = new TMatrix3D (structureB);
	matrixB2  = new TMatrix3D (structureB);
	matrixB3  = new TMatrix3D (structureB);
	matrixB1T  = new TMatrix3D (structureBT);
	matrixB2T  = new TMatrix3D (structureBT);
	matrixB3T  = new TMatrix3D (structureBT);
	
	sqmatrixA = sqmatrixA11;
	sqmatrixM = sqmatrixM11;
	break;
	
    default:
	cerr << __FILE__ ":" << __LINE__ << ": NSTYPE " << NSTYPE << " not yet implemented" << endl;
	exit(0);
  }
  
  Solver = new TPardisoSolver();
}

int MATRICES_NSTYPE::PrepareAssemble(int &N_FESP, TFESpace3D **FESP,
				      int &N_SQMAT, TSquareMatrix3D **SQMATRICES,
				      int &N_MAT, TMatrix3D **MATRICES)
{
  TFESpace *ansatzspace, *testspace;
  TStructure3D *structure;
  
  structure = matrixB1->GetStructure();
  ansatzspace = structure->GetAnsatzSpace();
  testspace   = structure->GetTestSpace();
  
  N_FESP = 2;
  FESP[0] = (TFESpace3D*) ansatzspace;
  FESP[1] = (TFESpace3D*) testspace;
  
  switch (NSTYPE)
  {
    case 1:
      N_SQMAT = 2;
      SQMATRICES[0] = sqmatrixA;
      SQMATRICES[1] = sqmatrixM;
      
      N_MAT = 3;
      MATRICES[0] = matrixB1;
      MATRICES[1] = matrixB2;
      MATRICES[2] = matrixB3;
        
      break;
      
    case 2:
      N_SQMAT = 2;
      SQMATRICES[0] = sqmatrixA;
      SQMATRICES[1] = sqmatrixM;
      
      N_MAT = 6;
      MATRICES[0] = matrixB1;
      MATRICES[1] = matrixB2;
      MATRICES[2] = matrixB3;
      MATRICES[3] = matrixB1T;
      MATRICES[4] = matrixB2T;
      MATRICES[5] = matrixB3T;
      
      break;
   
    case 4:
      N_SQMAT = 12;
      
      SQMATRICES[0] = sqmatrixA11;
      SQMATRICES[1] = sqmatrixA12;
      SQMATRICES[2] = sqmatrixA13;
      SQMATRICES[3] = sqmatrixA21;
      SQMATRICES[4] = sqmatrixA22;
      SQMATRICES[5] = sqmatrixA23;
      SQMATRICES[6] = sqmatrixA31;
      SQMATRICES[7] = sqmatrixA32;
      SQMATRICES[8] = sqmatrixA33;
      SQMATRICES[9] = sqmatrixM11;
      SQMATRICES[10] = sqmatrixM22;
      SQMATRICES[11] = sqmatrixM33;
      
      N_MAT = 6;
      MATRICES[0] = matrixB1;
      MATRICES[1] = matrixB2;
      MATRICES[2] = matrixB3;
      MATRICES[3] = matrixB1T;
      MATRICES[4] = matrixB2T;
      MATRICES[5] = matrixB3T;
      
      break;
      
  }
  
  for (int i=0;i<N_SQMAT;++i)
    SQMATRICES[i]->Reset();
  
  for (int i=0;i<N_MAT;++i)
  {
    MATRICES[i]->Reset();
  }
  
  return 2;
}

int MATRICES_NSTYPE::PrepareAssembleNL(int &N_FESP, TFESpace3D **FESP,
					int &N_SQMAT, TSquareMatrix3D **SQMATRICES,
					int &N_MAT, TMatrix3D **MATRICES)
{
  TFESpace3D *fespace;
  
  switch (NSTYPE)
  {
    case 1:
    case 2:
      fespace  = sqmatrixA->GetFESpace();
      
      N_FESP = 1;
      FESP[0] = fespace;
      
      N_SQMAT = 1;
      SQMATRICES[0] = sqmatrixA;
      SQMATRICES[0]->Reset();
      
      N_MAT = 0;
      
      break;
      
    case 4:
      fespace = sqmatrixA11->GetFESpace();
      
      N_FESP = 1;
      FESP[0] = fespace;
      
      N_SQMAT = 3;
      SQMATRICES[0] = sqmatrixA11;
      SQMATRICES[1] = sqmatrixA22;
      SQMATRICES[2] = sqmatrixA33;
      
      SQMATRICES[0]->Reset();
      SQMATRICES[1]->Reset();
      SQMATRICES[2]->Reset();
      
      N_MAT = 0;
      break;
  }
  
  return 1;
}

double MATRICES_NSTYPE::GetDefect(int N_, double *sol, double *rhs)
{
  double norm, *res;
  int N_U, N_P;
  
  res = new double [N_];
  
  switch (NSTYPE)
  {
    case 2:
      N_U = sqmatrixA->GetN_Rows();
      N_P = matrixB1->GetN_Rows();
      
      CoupledDefect(sqmatrixA, matrixB1, matrixB2, matrixB3,
		    matrixB1T, matrixB2T, matrixB3T, sol, rhs, res);
		
      break;
      
    case 4:
      N_U = sqmatrixA11->GetN_Rows();
      N_P = matrixB1->GetN_Rows();
      CoupledDefect(sqmatrixA11, sqmatrixA12, sqmatrixA13,
                   sqmatrixA21, sqmatrixA22, sqmatrixA23,
                   sqmatrixA31, sqmatrixA32, sqmatrixA33,
                   matrixB1, matrixB2, matrixB3,
                   matrixB1T, matrixB2T, matrixB3T,
                   sol, rhs, res);
      break;
      
    default:
      return -1.0;
  }
  
  if ( TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE )
    IntoL20Vector3D(res+3*N_U, N_P, TDatabase::ParamDB->PRESSURE_SPACE);
  
  norm = Dnorm(N_, res);

  delete [] res;
  
  return norm;
}

void MATRICES_NSTYPE::Solve(double *sol, double *rhs)
{  
  switch (NSTYPE)
  {
    case 4:
    case 2:
//       DirectSolver(sqmatrixA, matrixB1T, matrixB2T, matrixB3T,
// 		   matrixB1, matrixB2, matrixB3, rhs, sol);

      BuildSystem();
      Solver->FactorizeSolve(sol, rhs);

      break;
      
    default:
      OutPut(__FILE__ << ":" << __LINE__ << ": NSTYPE not yet implemented!" << endl);
  }
}

void MATRICES_NSTYPE::AddToRhs(double tau, double *sol, double *rhs)
{
  double *rhs_tmp;
  int N_U, N_Active;
  
  switch (NSTYPE)
  {
    case 2:
      N_U = sqmatrixM->GetN_Rows();
      
      rhs_tmp = new double [N_U];
      
      // add mass to rhs and scale rhs (active only)
      for (int i=0;i<3;++i)
      {
	N_Active = sqmatrixM->GetActiveBound();
	
	MatVectActive(sqmatrixM, sol + i*N_U, rhs_tmp);
	
	Dscal(N_Active, tau, rhs + i*N_U);
	Daxpy(N_Active, 1.0, rhs_tmp, rhs + i*N_U);      
      }
      
      break;
      
    case 4:
      N_U = sqmatrixM11->GetN_Rows();
      rhs_tmp = new double [N_U];
      
      N_Active = sqmatrixM11->GetActiveBound();
      
      MatVectActive(sqmatrixM11, sol, rhs_tmp);
      Dscal(N_Active, tau, rhs);
      Daxpy(N_Active, 1.0, rhs_tmp, rhs); 
      
      MatVectActive(sqmatrixM22, sol+N_U, rhs_tmp);
      Dscal(N_Active, tau, rhs+N_U);
      Daxpy(N_Active, 1.0, rhs_tmp, rhs+N_U); 
      
      MatVectActive(sqmatrixM33, sol+2*N_U, rhs_tmp);
      Dscal(N_Active, tau, rhs+2*N_U);
      Daxpy(N_Active, 1.0, rhs_tmp, rhs+2*N_U); 
      
      break;
    default:
      OutPut(__FILE__ << ":" << __LINE__ << ": NSTYPE not yet implemented!" << endl);
  }
  
  delete [] rhs_tmp;
}

void MATRICES_NSTYPE::AddMToA(double tau)
{
  int N_;
  
  switch (NSTYPE)
  {
    case 2:
      MatAdd2(sqmatrixA, sqmatrixM, tau);
      break;
      
    case 4:
      MatAdd2(sqmatrixA11, sqmatrixM11, tau);
      MatAdd2(sqmatrixA22, sqmatrixM22, tau);
      MatAdd2(sqmatrixA33, sqmatrixM33, tau);
      
      // scale non-diagonal blocks
      N_ = sqmatrixA11->GetN_Entries();
      
      Dscal(N_, tau, sqmatrixA12->GetEntries());
      Dscal(N_, tau, sqmatrixA13->GetEntries());
      
      Dscal(N_, tau, sqmatrixA21->GetEntries());
      Dscal(N_, tau, sqmatrixA23->GetEntries());
      
      Dscal(N_, tau, sqmatrixA31->GetEntries());
      Dscal(N_, tau, sqmatrixA32->GetEntries());
      
      
      break;
  }
}

void MATRICES_NSTYPE::AddMToA_diag(double tau)
{
  switch (NSTYPE)
  {
    case 2:
      MatAdd2(sqmatrixA, sqmatrixM, tau);
      break;
      
    case 4:
      MatAdd2(sqmatrixA11, sqmatrixM11, tau);
      MatAdd2(sqmatrixA22, sqmatrixM22, tau);
      MatAdd2(sqmatrixA33, sqmatrixM33, tau);
      
      break;
  }
}


void MATRICES_NSTYPE::BuildSystem()
{
  switch (NSTYPE)
  {
    case 2:
      Solver->SetMatrix(sqmatrixA, matrixB1T, matrixB2T, matrixB3T,
		       matrixB1, matrixB2, matrixB3);
      break;
      
    case 4:
      Solver->SetMatrix(sqmatrixA11, sqmatrixA12, sqmatrixA13, sqmatrixA21, sqmatrixA22,
			sqmatrixA23, sqmatrixA31, sqmatrixA32, sqmatrixA33,
			matrixB1T, matrixB2T, matrixB3T,
		        matrixB1, matrixB2, matrixB3);
      break;
      
    default:
      OutPut(__FILE__ << ":" << __LINE__ << ": NSTYPE not yet implemented!" << endl);
  }
}

void MATRICES_NSTYPE::Scale(double factor)
{  
  switch (NSTYPE)
  {
    case 4:
    case 2:
//       Dscal(matrixB1->GetN_Entries(), factor, matrixB1->GetEntries());
//       Dscal(matrixB2->GetN_Entries(), factor, matrixB2->GetEntries());
//       Dscal(matrixB3->GetN_Entries(), factor, matrixB3->GetEntries());
      Dscal(matrixB1T->GetN_Entries(), factor, matrixB1T->GetEntries());
      Dscal(matrixB2T->GetN_Entries(), factor, matrixB2T->GetEntries());
      Dscal(matrixB3T->GetN_Entries(), factor, matrixB3T->GetEntries()); 
      
      break;
    default:
      OutPut(__FILE__ << ":" << __LINE__ << ": NSTYPE not yet implemented!" << endl);
  }
}

void MATRICES_NSTYPE::Print()
{
  int N_SQMAT, N_MAT, N_, N_Rows;
  int begin, end, pos;
  int *KCol, *RowPtr;
  double *Entries;
  
  char A[] = "A";
  char B1[] = "B1";
  char B2[] = "B2";
  char B3[] = "B3";
  char B1T[] = "B1T";
  char B2T[] = "B2T";
  char B3T[] = "B3T";
  
  char *Names[7] = { A, B1, B2, B3, B1T, B2T, B3T };
  
  TSquareMatrix3D *SqMatrices[4];
  TMatrix3D       *Matrices[6];
  
  switch (NSTYPE)
  {
    case 1:
      N_SQMAT = 1;
      N_MAT   = 3;
      
      SqMatrices[0] = sqmatrixA;
      
      Matrices[0] = matrixB1;
      Matrices[1] = matrixB2;
      Matrices[2] = matrixB3;
      
      break;
      
    case 2:
      N_SQMAT = 1;
      N_MAT   = 6;
      
      SqMatrices[0] = sqmatrixA;
      
      Matrices[0] = matrixB1;
      Matrices[1] = matrixB2;
      Matrices[2] = matrixB3;
      Matrices[3] = matrixB1T;
      Matrices[4] = matrixB2T;
      Matrices[5] = matrixB3T;
      
      break;
  }
  
  N_ = N_SQMAT + N_MAT;
  for (int i=0;i<N_SQMAT;++i)
  {
    KCol = SqMatrices[i]->GetKCol();
    RowPtr = SqMatrices[i]->GetRowPtr();
    Entries = SqMatrices[i]->GetEntries();
    N_Rows = SqMatrices[i]->GetN_Rows();
    
    for (int j=0;j<N_Rows;++j)
    {
      begin = RowPtr[j];
      end   = RowPtr[j+1];
      
      pos = begin;
      for (int k=begin;k<end;++k)
      {
	cout << "\t" << Names[i] << "(" << j << "," << KCol[pos] << ")";
	cout << " = " << Entries[pos] << endl;
	
	++pos;
      }
    }
  }
  
  for (int i=0;i<N_MAT;++i)
  {
    KCol = Matrices[i]->GetKCol();
    RowPtr = Matrices[i]->GetRowPtr();
    Entries = Matrices[i]->GetEntries();
    N_Rows = Matrices[i]->GetN_Rows();
    
    for (int j=0;j<N_Rows;++j)
    {
      begin = RowPtr[j];
      end   = RowPtr[j+1];
      
      pos = begin;
      for (int k=begin;k<end;++k)
      {
	cout << "\t" << Names[N_SQMAT+i] << "(" << j << "," << KCol[pos] << ")";
	cout << " = " << Entries[pos] << endl;
	
	++pos;
      }
    }
  }
  
  cout << endl;
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

void WriteGrid(const char *name, TDomain *Domain)
{
  TCollection *Coll = Domain->GetCollection(It_Finest, 0);
  TOutput3D Output(0,0,0,0,Domain, Coll);
  
  Output.WriteVtk(name);
  
//   PrintMesh(Coll);
  
  delete Coll;
}

void WriteGrid(const char *name, TDomain *Domain, TCollection *Coll)
{
  TOutput3D Output(0,0,0,0,Domain, Coll);
  
  Output.WriteVtk(name);
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

void MATRICES_GRID::Init (TSquareStructure3D *structure)
{
  A11 = new TSquareMatrix3D (structure);
  A12 = new TSquareMatrix3D (structure);
  A13 = new TSquareMatrix3D (structure);
  A21 = new TSquareMatrix3D (structure);
  A22 = new TSquareMatrix3D (structure);
  A23 = new TSquareMatrix3D (structure);
  A31 = new TSquareMatrix3D (structure);
  A32 = new TSquareMatrix3D (structure);
  A33 = new TSquareMatrix3D (structure);
}

void MATRICES_GRID::PrepareAssemble(int &N_FESP, TFESpace3D **FESP,
				      int &N_SQMAT, TSquareMatrix3D **SQMATRICES,
				      int &N_MAT, TMatrix3D **MATRICES)
{
  N_FESP = 1;
  FESP[0] = A11->GetFESpace();
  
  N_SQMAT = 9;
  SQMATRICES[0] = A11;
  SQMATRICES[1] = A12;
  SQMATRICES[2] = A13;
  SQMATRICES[3] = A21;
  SQMATRICES[4] = A22;
  SQMATRICES[5] = A23;
  SQMATRICES[6] = A31;
  SQMATRICES[7] = A32;
  SQMATRICES[8] = A33;
  
  N_MAT = 0;
  
  for (int i=0;i<N_SQMAT;++i)
    SQMATRICES[i]->Reset();
}

void MATRICES_GRID::GetMatrices(TSquareMatrix3D **SQMATRICES)
{
  SQMATRICES[0] = A11;
  SQMATRICES[1] = A12;
  SQMATRICES[2] = A13;
  SQMATRICES[3] = A21;
  SQMATRICES[4] = A22;
  SQMATRICES[5] = A23;
  SQMATRICES[6] = A31;
  SQMATRICES[7] = A32;
  SQMATRICES[8] = A33;
}

void MATRICES_GRID::Solve(double *sol, double *rhs)
{
  TSquareMatrix3D *SQMATRICES[9];
  int N_Grid = A11->GetN_Rows();
  
  SQMATRICES[0] = A11;
  SQMATRICES[1] = A12;
  SQMATRICES[2] = A13;
  SQMATRICES[3] = A21;
  SQMATRICES[4] = A22;
  SQMATRICES[5] = A23;
  SQMATRICES[6] = A31;
  SQMATRICES[7] = A32;
  SQMATRICES[8] = A33;
  
  if ( TDatabase::ParamDB->P1 > 0 )
  {
    OutPut("solve grid equation ... ");
    DirectSolver(SQMATRICES, 3, 3, sol, rhs);
  //   memcpy(sol, rhs, 3*N_Grid*sizeof(double));
    OutPut("done" << endl);
  }
  else if ( TDatabase::ParamDB->P1 == 0 )
  {
    memcpy(sol, rhs, 3*N_Grid*sizeof(double));
  }
  else 
  {
    ;
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
      OutPut("Mesh check failed! Exit" << endl);
      exit(0);
    }
  }
}

void WriteSurface(TDomain *Domain, TCollection *Coll, int N_BoundFaces, int *CellNumbers)
{
  int counter;
  TBaseCell **Cells, *Cell;
  TCollection *SurfColl;
  
  Cells = new TBaseCell* [N_BoundFaces];
  
  for (int i=0;i<N_BoundFaces;++i)
  {
    Cells[i] = NULL;
    Coll->GetCell(CellNumbers[i])->SetClipBoard(-1);
  }
  
  counter = 0;
  for (int i=0;i<N_BoundFaces;++i)
  {
    Cell = Coll->GetCell(CellNumbers[i]);
    
    if ( Cell->GetClipBoard() == -1 )
    {
      Cells[counter] = Cell;
      Cell->SetClipBoard(0);
      ++counter;
    }
  } 
  
  SurfColl = new TCollection (counter, Cells);
  
  WriteGrid("surf.vtk", Domain, SurfColl);
  
  delete SurfColl;
//   delete [] Cells;
}

void PrintRhs(TCollection *Coll, int N_BoundFaces, int *CellNumbers, int *JointNumbers, TFESpace3D *fesp,
	      double *rhs, double *rhs_ref)
{
  int *GlobalNumbers, *BeginIndex, *DOF, N_DOF, N_U;
  int N_JointDOF, *JointDOF, N_Joints;
  TBaseCell *Cell;
  FE3D FeID;
  TFEDesc3D *FEDesc;
  
  GlobalNumbers = fesp->GetGlobalNumbers();
  BeginIndex    = fesp->GetBeginIndex();
  N_U = fesp->GetN_DegreesOfFreedom();
  
  std::ofstream dat("a_rhs_cmp.txt");
  
  for (int i=0;i<N_BoundFaces;++i)
  {
    Cell = Coll->GetCell(CellNumbers[i]);
    DOF = GlobalNumbers + BeginIndex[CellNumbers[i]];
    
    FeID = fesp->GetFE3D(CellNumbers[i], Cell);
    FEDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FeID);
    
    N_DOF = FEDesc->GetN_DOF();
    N_JointDOF = FEDesc->GetN_JointDOF();
    N_Joints = Cell->GetN_Joints();
    
    for (int j=0;j<N_Joints;++j)
    {
      if ( j == JointNumbers[i] )
      {
	JointDOF = FEDesc->GetJointDOF(j);
	
	dat << "Cell: " << CellNumbers[i] << " " << j << endl;
	for (int k=0;k<N_JointDOF;++k)
	{
	  dat << k << ": ";
	  dat << rhs_ref[DOF[JointDOF[k]]] << "  |  " << rhs[DOF[JointDOF[k]]];
	  dat << "  |  " << rhs[DOF[JointDOF[k]]]/rhs_ref[DOF[JointDOF[k]]];
	  dat << "  ||  ";
	  dat << rhs_ref[DOF[JointDOF[k]]+N_U] << "  |  " << rhs[DOF[JointDOF[k]]+N_U];
	  dat << "  |  " << rhs[DOF[JointDOF[k]]+N_U]/rhs_ref[DOF[JointDOF[k]]+N_U];
	  dat << "  ||  ";
	  dat << rhs_ref[DOF[JointDOF[k]]+2*N_U] << "  |  " << rhs[DOF[JointDOF[k]]+2*N_U];
	  dat << "  |  " << rhs[DOF[JointDOF[k]]+2*N_U]/rhs_ref[DOF[JointDOF[k]]+2*N_U];
	  dat << endl;
	  
// 	  if ( k == 0 || k == 2 || k == 5 )
// 	  {
// 	    rhs[DOF[JointDOF[k]]] = rhs_ref[DOF[JointDOF[k]]];
// 	    rhs[DOF[JointDOF[k]]+N_U] = rhs_ref[DOF[JointDOF[k]]+N_U];
// 	    rhs[DOF[JointDOF[k]]+2*N_U] = rhs_ref[DOF[JointDOF[k]]+2*N_U];
// 	  }
	  
	}      
      }
    }
  }
  
  dat.close();
}

void MakeIsoparametric_Sphere(TCollection *Coll, int N_BoundFaces, int *CellNumbers,
			      int *JointNumbers, TFESpace3D *fespace,
			      TVertex **&AllVertices, int *&VertexCell, int &N_AllVertices)
{
  int CellNr, JointNr, N_JointDOF, *JointDOF, N_AddVertex, N_Points;
  double *xi, *eta, *zeta, x, y, z, norm;
  double param1[4], param2[4];
  TBoundComp3D *BdComp;
  TBaseCell *Cell;
  TJoint *Joint, *NewJoint;
  FE3D FeID;
  TFEDesc3D *FEDesc;
  TNodalFunctional3D *NF;
  TVertex *AddVertex[64], **Vertices;
  TTetraAffin *RefTrans;
  
  const int *TmpFV, *TmpLen;
  int MaxLen, counter, len, N_Vertices, hash, a, b, c;
  int N_Edges;
  
  int *EdgeList, **EdgeHash, *BucketCount, *Bucket;
  
  OutPut("Change to isoparametric elements" << endl);
  
//   std::ofstream dat ("scatter.dat");
  
  EdgeList = new int [6*N_BoundFaces];
  
  // reset clipboard
  for (int i=0;i<N_BoundFaces;++i)
  {
    CellNr = CellNumbers[i];
    JointNr = JointNumbers[i];
    Cell = Coll->GetCell(CellNr);
    
    for (int j=0;j<4;++j)
    {
      Cell->GetVertex(j)->SetClipBoard(-1);
    }
  }
  
  // count vertices
  counter = 0;
  for (int i=0;i<N_BoundFaces;++i)
  {
    CellNr = CellNumbers[i];
    JointNr = JointNumbers[i];
    Cell = Coll->GetCell(CellNr);
    
    Cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
    len  = TmpLen[JointNr];
    
    for (int j=0;j<len;++j)
    {
      Cell->GetVertex(TmpFV[JointNr*MaxLen+j])->SetClipBoard(counter);
      counter++;
    }
  }
  
  N_Vertices = counter;
  EdgeHash = new int* [2*N_Vertices];
  BucketCount = new int [2*N_Vertices];
  
  memset(EdgeHash, 0, 2*N_Vertices*sizeof(EdgeHash[0]));
  memset(BucketCount, 0, 2*N_Vertices*sizeof(BucketCount[0]));
  
  // get bucket count
  for (int i=0;i<N_BoundFaces;++i)
  {
    CellNr = CellNumbers[i];
    JointNr = JointNumbers[i];
    Cell = Coll->GetCell(CellNr);
    
    Cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
    len = TmpLen[JointNr];
    
    for (int j=0;j<len;++j)
    {
      hash = 	Cell->GetVertex(TmpFV[JointNr*MaxLen+j])->GetClipBoard()
	      + Cell->GetVertex(TmpFV[JointNr*MaxLen+((j+1)%len)])->GetClipBoard();
	      
      BucketCount[hash]++;
    }
  }
  
  // init hashtable and edge list
  counter = 0;
  for (int i=0;i<N_BoundFaces;++i)
  {
    CellNr = CellNumbers[i];
    JointNr = JointNumbers[i];
    Cell = Coll->GetCell(CellNr);
    
    Cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
    len = TmpLen[JointNr];
    
    for (int j=0;j<len;++j)
    {
      a = Cell->GetVertex(TmpFV[JointNr*MaxLen+j])->GetClipBoard();
      b = Cell->GetVertex(TmpFV[JointNr*MaxLen+((j+1)%len)])->GetClipBoard();
	      
      hash = a + b;
      
      // sort 
      if ( b < a )
      {
	c = a;
	a = b;
	b = c;
      }
      
      if ( EdgeHash[hash] == NULL )
      {
	EdgeHash[hash] = new int [BucketCount[hash]];
	for (int k=0;k<BucketCount[hash];++k)
	  EdgeHash[hash][k] = -1;
	
	EdgeList[2*counter  ] = a;
	EdgeList[2*counter+1] = b;
	EdgeHash[hash][0] = counter;
	counter++;
      }
      else 
      {
	int found = 0, k;
	
	Bucket = EdgeHash[hash];
	for (k=0;Bucket[k] >= 0;++k)
	{
	  int index = Bucket[k];
	  
	  if ( EdgeList[2*index] == a && EdgeList[2*index+1] == b )
	  {
	    found = 1;
	    break;
	  }
	}
	
	if ( found == 0 )
	{
	  EdgeList[counter  ] = a;
	  EdgeList[counter+1] = b;
	  
	  Bucket[k] = counter;
	  counter++;
	}
      }
    }
  } // end init hastable
  
  N_Edges = counter;
  AllVertices = new TVertex* [3*N_BoundFaces];
  VertexCell = new int [N_BoundFaces];
  
  N_AllVertices = 3*N_BoundFaces;
  memset(AllVertices, 0, N_Edges*sizeof(AllVertices[0]));
  
  // set joints
  for (int i=0;i<N_BoundFaces;++i)
  {
    CellNr = CellNumbers[i];
    JointNr = JointNumbers[i];
    Cell = Coll->GetCell(CellNr);
    Joint = Cell->GetJoint(JointNr);
    
    FeID = fespace->GetFE3D(CellNr, Cell);
    
    if ( ! (FeID == C_P2_3D_T_A || FeID == C_B2_3D_T_A) )
    {
      cerr << __FILE__ << ":" << __LINE__ << ": ";
      cerr << "error change to isoparametric: only P2 and P2+ allowed! Exit" << endl;
      exit(0);
    }
    
    FEDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FeID);
    NF     = TFEDatabase3D::GetNodalFunctional3DFromFE3D(FeID);
    RefTrans = (TTetraAffin*) TFEDatabase3D::GetRefTrans3D(TetraAffin);
    N_JointDOF = FEDesc->GetN_JointDOF();
    JointDOF = FEDesc->GetJointDOF(JointNr);
    
    RefTrans->SetCell(Cell);
    
    NF->GetPointsForFace(JointNr, N_Points, xi, eta, zeta);
    
    N_AddVertex = 0;
    for (int j=0;j<N_Points;++j)
    {
      if (j == 0 || j == 2 || j >= 5 )
      {
// 	AddVertex[N_AddVertex] = NULL;
// 	++N_AddVertex;
	continue;
      }
      else
      {
      
	RefTrans->GetOrigFromRef(xi[j], eta[j], zeta[j], x, y, z);
	
// 	norm  = sqrt(x*x + y*y + z*z);
// 	
// 	x /= norm;
// 	y /= norm;
// 	z /= norm;
	
	AddVertex[N_AddVertex] = new TVertex(x,y,z);
	
	++N_AddVertex;
	
// 	dat << x << " " << y << " " << z << endl;
	
// 	OutPut("norm: " << norm << " > " << sqrt(x*x + y*y + z*z) <<  endl);
      }
    }
    
    BdComp = ((TBoundFace*) Joint)->GetBoundComp();
    ((TBoundFace *) Joint)->GetParameters(param1, param2);
    
    NewJoint = new TIsoBoundFace (BdComp, param1, param2);
    
    Vertices = new TVertex* [N_AddVertex];
    for (int j=0;j<N_AddVertex;++j)
    {
      Vertices[j] = AddVertex[j];
      AllVertices[3*N_AddVertex+j] = AddVertex[j];
    }
    VertexCell[i] = CellNr;
    
    ((TIsoBoundFace*) NewJoint)->SetVertices(N_AddVertex, Vertices);
    
    // exchange joints
    Cell->SetJoint(JointNr, NewJoint);
    
    delete Joint;
  }
  
  delete [] EdgeList;
  delete [] EdgeHash;
  delete [] BucketCount;
//   dat.close();
}

void CheckIso(TCollection *Coll, int N_BoundFaces, int *CellNumbers, int *JointNumbers, int img)
{
  int CellNr, JointNr;
  int N_ = 100;
  TBaseCell *Cell;
  RefTrans3D RefTrans;
  TTetraIsoparametric *F_K;
  TTetraAffin *F_K2;
  
  double Coords[4][3] = { { 0.0, 0.0, 0.0 },
		      { 1.0, 0.0, 0.0 },
		      { 0.0, 1.0, 0.0 },
		      { 0.0, 0.0, 1.0 } };
  
  std::ostringstream os;
  
  os << "cell" << ".";
  os.width(5);
  os.fill('0');
  os << img << ".dat" << ends;
  
  std::ofstream dat(os.str().c_str());
  
  for (int i=0;i<N_BoundFaces;++i)
  {
    CellNr = CellNumbers[i];
    JointNr = JointNumbers[i];
    
    Cell = Coll->GetCell(CellNr);
    
    RefTrans = TetraIsoparametric;
    F_K = (TTetraIsoparametric*) TFEDatabase3D::GetRefTrans3D(RefTrans);
    F_K2 = (TTetraAffin*) TFEDatabase3D::GetRefTrans3D(TetraAffin);
    
    F_K->SetApproximationOrder(2);
    F_K->SetQuadFormula(P2Tetra);
    F_K->SetCell(Cell);
    F_K2->SetCell(Cell);
    
    double XI[N_], ETA[N_], ZETA[N_];
    double x, y, z, xi, eta, zeta;
    
    // points
    for (int j=0;j<4;++j)
    {
      Cell->GetVertex(j)->GetCoords(x,y,z);
      
      dat << x << " " << y << " " << z << endl;
    }
    
   
    for(int j=0;j<N_;++j)
    {
      XI[j] = ETA[j] = ZETA[j] = ( 1.0 + j ) / N_;
    }
    
    // face
    for (int j=0;j<N_;++j)
    {
      xi = XI[j];
      
      for (int k=0;k<N_-j;++k)
      {
	eta = ETA[k];
	
	F_K->GetOrigBoundFromRef(JointNr, xi, eta, x, y, z);
	
	dat << x << " " << y << " " << z << endl;
      }
    }
    
    break;
  }
  
  dat.close();
}

void MoveIsoVertices(TFEVectFunct3D *velocity, int N_BoundFaces, int *CellNumbers, int *JointNumbers,
		     double tau)
{
  int CellNr, JointNr, N_Points, N_IsoVertices, Dim, counter, N_U;
  int N_JointDOF, *JointDOF, *GlobalNumbers, *BeginIndex, *DOF;
  double *xi, *eta, *zeta, x, y, z, xvelo, yvelo, zvelo;
  double uref[MaxN_BaseFunctions3D];
  double *Values;
//   double xi[MaxN_BaseFunctions3D], eta[MaxN_BaseFunctions3D], zeta[MaxN_BaseFunctions3D];
  TBaseCell *Cell;
  TJoint *Joint;
  FE3D FeID;
  
  TFESpace3D *fespace;
  TCollection *Coll;
  TVertex **IsoVertices;
  TNodalFunctional3D *nf;
  TBaseFunct3D *bf;
  
  fespace = velocity->GetFESpace3D();
  Coll = fespace->GetCollection();
  N_U = fespace->GetN_DegreesOfFreedom();
  Values = velocity->GetValues();
  
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex    = fespace->GetBeginIndex();
  
  for (int i=0;i<N_BoundFaces;++i)
  {
    CellNr = CellNumbers[i];
    JointNr = JointNumbers[i];
    Cell = Coll->GetCell(CellNr);
    
    FeID = fespace->GetFE3D(CellNr, Cell);
    nf = TFEDatabase3D::GetNodalFunctional3DFromFE3D(FeID);
    bf = TFEDatabase3D::GetBaseFunct3DFromFE3D(FeID);
    Dim = bf->GetDimension();
    
    Joint = Cell->GetJoint(JointNr);
    
    DOF = GlobalNumbers + BeginIndex[CellNr];
    
    if ( Joint->GetType() == IsoBoundFace )
    {
      IsoVertices = ((TIsoBoundFace*) Joint)->GetVertices();
      N_IsoVertices = ((TIsoBoundFace*) Joint)->GetN_Vertices();
      
      nf->GetPointsForFace(JointNr, N_Points, xi, eta, zeta);
      
      counter = 0;
      for (int j=0;j<N_Points;++j)
      {
	if ( j == 0 || j == 2 || j >= 5 ) continue;
	
	bf->GetDerivatives(D000, xi[j], eta[j], zeta[j], uref);
	
	IsoVertices[counter]->GetCoords(x,y,z);
	
	xvelo = yvelo = zvelo = 0;
	for (int k=0;k<Dim;++k)
	{
	  xvelo += uref[k]*Values[DOF[k]      ];
	  yvelo += uref[k]*Values[DOF[k]+  N_U];
	  zvelo += uref[k]*Values[DOF[k]+2*N_U];	  
	}
	
	x += tau*xvelo;
	y += tau*yvelo;
	z += tau*zvelo;
	
	IsoVertices[counter]->SetCoords(x,y,z);
	++counter;
      }
    }
  }
}
