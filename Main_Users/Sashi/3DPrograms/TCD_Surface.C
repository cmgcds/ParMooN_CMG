#include <MooNMD.h>

#include <TCD3D_Surface.h>

#include "../Examples/TNSE_3D/TwoPhase_NavierStokes.h"
#include "../Examples/TNSE_3D/TwoPhase_SurfaceEq.h"

#define Print(x) OutPut(#x<<": "<<x<<endl)

void CheckNeigbourhood(TCollection *SurfColl);

typedef struct _MUTABLE_DATA
{
  TCollection *Coll, *Coll_P0, *Coll_P1, *SurfColl;
  
  int *GlobalCellNo_P0, *GlobalCellNo_P1;
  int *CellNumbers_P0, *CellNumbers_P1;
  int *JointNumbers_P0, *JointNumbers_P1;
  int N_SurfaceJoints_P0, N_SurfaceJoints_P1;
  
  TFEFunction3D *phase_marker;
  
  TFESpace3D *velocity_space, *concen_space;
  TFESpace3D *grid_space, *grid_space_p0, *grid_space_p1;
  
  TFESpace2D *concen_space_surf;
  
  double *sol;
  double *sol_surf, *rhs_surf, *tmp_surf;
  double *grid_velo, *grid_velo_p0, *grid_velo_p1;
  double *grid_rhs_p0, *grid_rhs_p1;
  double *grid_pos;
  
  int N_C_Surf, N_Grid, N_Grid_P0, N_Grid_P1, N_U;
  
  TFEFunction2D *c_surf;
  TFEVectFunct3D *u, *g, *grid;
  TFEFunction3D *ux, *uy, *uz;
  
  TSquareStructure2D *sqstructure_surf;
  TSquareStructure3D *sqstruct_grid_p0, *sqstruct_grid_p1;
  
  TSquareMatrix2D *A_Surf, *M_Surf;
  
  TSquareMatrix3D *G11_P0, *G12_P0, *G13_P0;
  TSquareMatrix3D *G21_P0, *G22_P0, *G23_P0;
  TSquareMatrix3D *G31_P0, *G32_P0, *G33_P0;
  
  TSquareMatrix3D *G11_P1, *G12_P1, *G13_P1;
  TSquareMatrix3D *G21_P1, *G22_P1, *G23_P1;
  TSquareMatrix3D *G31_P1, *G32_P1, *G33_P1;
  
  TAuxParam2D3D *AuxParam_Surf;
  
  TOutput2D *Output_Surf;
  TOutput3D *Output;
  
  void Create(TDomain *Domain);
  
} MUTABLE_DATA;

void WriteSolution(char *basename, int number, TOutput2D *out);

int main(int argc, char **argv)
{
  TDatabase Database;
  TFEDatabase3D FEDatabase3D;
  TFEDatabase2D FEDatabase2D;
  TDomain Domain;
  
  MUTABLE_DATA Data;
  
  int N_FESP, N_SQMAT, N_MAT, img=0;
  double mass, massO, massloss, errors[10], error_L2=0, error_loc;
  TFESpace3D *FESP[10];
  TFESpace2D *FESPSURF[10];
  TSquareMatrix3D *SQMATRICES[10];
  TSquareMatrix2D *SQMATRICESSURF[10];
  
  char Name[] = "name";
  char Desc[] = "desc";
  char Inter[] = "IsoSurface";
  
  if ( argc > 1 )
    Domain.ReadParam(argv[1]);
  else
  {
    cerr << "no readin file given" << endl;
    exit(0);
  }
  
  char *VTK = TDatabase::ParamDB->VTKBASENAME;
//   Print(VTK);exit(0);
  
  BoundCondFunct3D   *GRIDBOUNDCOND[4];
  BoundValueFunct3D  *GRIDBOUNDVALUE[4];
  
  GRIDBOUNDCOND[0] = GridBoundCond;
  GRIDBOUNDCOND[1] = GridBoundCond;
  GRIDBOUNDCOND[2] = GridBoundCond;
  
  GRIDBOUNDVALUE[0] = GridXBoundValues;
  GRIDBOUNDVALUE[1] = GridYBoundValues;
  GRIDBOUNDVALUE[2] = GridZBoundValues;
  
  /// create domain and initialise mutable data
  Domain.Tetgen(TDatabase::ParamDB->SMESHFILE);
  Data.Create(&Domain);
 
//   PRINTVAR(Data.concen_space_surf->GetActiveBound());
// exit(0);  
  /// create discrete forms
  TDiscreteForm3D *DiscreteForm = new TDiscreteForm3D (Name, Desc, TCD3D_Surf::N_Terms,
				      TCD3D_Surf::Derivatives, TCD3D_Surf::FESpaceNumbers,
				      TCD3D_Surf::N_Matrices, TCD3D_Surf::N_Rhs,
				      TCD3D_Surf::RowSpaces, TCD3D_Surf::ColumnSpaces,
				      TCD3D_Surf::RhsSpaces, TCD3D_Surf::MatricesAssemble, 
				      SurfCoeffs, NULL);
				      
  TDiscreteForm3D *DiscreteFormGrid = new TDiscreteForm3D (Name, Desc, ESGridN_Terms, ESGridDerivatives,
					  ESGridSpaceNumbers, ESGridN_Matrices,
					  ESGridN_Rhs, ESGridRowSpace, ESGridColumnSpace,
					  ESGridRhsSpace, ESGridAssemble, GridCoeffs, NULL);
					  
  // dummy aux
  TAuxParam3D aux_dummy (0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL);
				      
  double tau = TDatabase::TimeDB->TIMESTEPLENGTH;
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  
  InterpolateSurface(Data.c_surf, Data.Coll_P0, Data.CellNumbers_P0, Data.JointNumbers_P0,
		     InitialC_Surf);
  
  mass = CalculateSurfaceMass(Data.c_surf, Data.concen_space, Data.Coll_P0, Data.CellNumbers_P0,
			      Data.JointNumbers_P0);
			      
  PRINTVAR(mass);
  massO = mass;
  
  CalcSurfaceErrors ( Data.c_surf, Data.concen_space, Data.CellNumbers_P0, Data.JointNumbers_P0,
			ExactC_Surf, errors);
			
  error_loc = errors[0];
  PRINTVAR(error_loc);
  
  WriteSolution(VTK, img, Data.Output_Surf);
  DumpIsosurface(Inter, img, Data.c_surf,
		  Data.Coll_P0, Data.N_SurfaceJoints_P0, Data.CellNumbers_P0,
		  Data.JointNumbers_P0);
  ++img;
// exit(0);		
  SetNormals(Data.Coll_P1, Data.N_SurfaceJoints_P1, Data.CellNumbers_P1, Data.JointNumbers_P1);
  
  N_FESP = 2;
  FESP[0] = Data.concen_space;
  FESP[1] = Data.velocity_space;
  FESPSURF[0] = Data.concen_space_surf;
  FESPSURF[1] = NULL;
    
  N_SQMAT = 2;
  SQMATRICESSURF[0] = Data.A_Surf;
  SQMATRICESSURF[1] = Data.M_Surf;
  SQMATRICESSURF[0]->Reset();
  SQMATRICESSURF[1]->Reset();
    
  SurfaceAssemble3D(N_FESP, FESPSURF, FESP,
		      Data.CellNumbers_P0, Data.JointNumbers_P0, Data.GlobalCellNo_P0,
		      N_SQMAT, SQMATRICESSURF, DiscreteForm,
		      Data.AuxParam_Surf);
  
  std::ofstream dat("mass.dat");
  dat << TDatabase::TimeDB->CURRENTTIME << " " << mass << endl;
  
  for (int timestep=0;TDatabase::TimeDB->CURRENTTIME<TDatabase::TimeDB->ENDTIME;++timestep)
  {
    TDatabase::TimeDB->CURRENTTIME += tau;
    
    OutPut("=========================================="<<endl);
    OutPut("CURRENTTIME: "<<TDatabase::TimeDB->CURRENTTIME<<endl);
    OutPut("=========================================="<<endl);
    
    memset(Data.grid_rhs_p0, 0, 3*Data.N_Grid_P0*sizeof(Data.grid_rhs_p0[0]));
    memset(Data.grid_rhs_p1, 0, 3*Data.N_Grid_P1*sizeof(Data.grid_rhs_p1[0]));
    
    Data.ux->Interpolate(ExactG1);
    Data.uy->Interpolate(ExactG2);
    Data.uz->Interpolate(ExactG3);
        
    // grid phase 0;
    FESP[0] = Data.grid_space_p0;
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
		
    // grid phase 1
    FESP[0] = Data.grid_space_p1;
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
    
    /// move grid
    SetGridRhs(Data.u, Data.grid_space_p0, Data.N_SurfaceJoints_P0, Data.CellNumbers_P0,
	       Data.JointNumbers_P0, Data.GlobalCellNo_P0, tau, 
	       Data.grid_rhs_p0, Data.grid_rhs_p0+Data.N_Grid_P0, Data.grid_rhs_p0+2*Data.N_Grid_P0);
	       
    SetGridRhs(Data.u, Data.grid_space_p1, Data.N_SurfaceJoints_P1, Data.CellNumbers_P1,
	       Data.JointNumbers_P1, Data.GlobalCellNo_P1, tau, 
	       Data.grid_rhs_p1, Data.grid_rhs_p1+Data.N_Grid_P1, Data.grid_rhs_p1+2*Data.N_Grid_P1);
	       
    SQMATRICES[0] = Data.G11_P0;
    SQMATRICES[1] = Data.G12_P0;
    SQMATRICES[2] = Data.G13_P0;
    SQMATRICES[3] = Data.G21_P0;
    SQMATRICES[4] = Data.G22_P0;
    SQMATRICES[5] = Data.G23_P0;
    SQMATRICES[6] = Data.G31_P0;
    SQMATRICES[7] = Data.G32_P0;
    SQMATRICES[8] = Data.G33_P0;
    
    DirectSolver(SQMATRICES, 3, 3, Data.grid_velo_p0, Data.grid_rhs_p0);
    
    SQMATRICES[0] = Data.G11_P1;
    SQMATRICES[1] = Data.G12_P1;
    SQMATRICES[2] = Data.G13_P1;
    SQMATRICES[3] = Data.G21_P1;
    SQMATRICES[4] = Data.G22_P1;
    SQMATRICES[5] = Data.G23_P1;
    SQMATRICES[6] = Data.G31_P1;
    SQMATRICES[7] = Data.G32_P1;
    SQMATRICES[8] = Data.G33_P1;
    
    DirectSolver(SQMATRICES, 3, 3, Data.grid_velo_p1, Data.grid_rhs_p1);
    
    MapGridVelo(Data.grid_space_p0, Data.g, Data.GlobalCellNo_P0, Data.grid_velo_p0);
    MapGridVelo(Data.grid_space_p1, Data.g, Data.GlobalCellNo_P1, Data.grid_velo_p1);
//     Dscal(3*Data.N_Grid, 1/tau, Data.grid_velo);    
    Daxpy(3*Data.N_Grid, 1, Data.grid_velo, Data.grid_pos); 
    Data.grid->DataToGrid();
    MoveIsoVertices(Data.u, Data.N_SurfaceJoints_P1, Data.CellNumbers_P1, Data.JointNumbers_P1,
		    Data.GlobalCellNo_P1, tau);
    
    SetNormals(Data.Coll_P1, Data.N_SurfaceJoints_P1, Data.CellNumbers_P1, Data.JointNumbers_P1);
    
    ///===================================================================================
    memset(Data.rhs_surf, 0, Data.N_C_Surf*sizeof(Data.rhs_surf[0]));
    MatVect(Data.M_Surf, Data.sol_surf, Data.tmp_surf);
    
    N_FESP = 2;
    FESP[0] = Data.concen_space;
    FESP[1] = Data.velocity_space;
    FESPSURF[0] = Data.concen_space_surf;
    FESPSURF[1] = NULL;
    
    N_SQMAT = 2;
    SQMATRICESSURF[0] = Data.A_Surf;
    SQMATRICESSURF[1] = Data.M_Surf;
    SQMATRICESSURF[0]->Reset();
    SQMATRICESSURF[1]->Reset();
    
    SurfaceAssemble3D(N_FESP, FESPSURF, FESP,
		      Data.CellNumbers_P0, Data.JointNumbers_P0, Data.GlobalCellNo_P0,
		      N_SQMAT, SQMATRICESSURF, DiscreteForm,
		      Data.AuxParam_Surf);
		      
//     MatVect(Data.M_Surf, Data.sol_surf, Data.tmp_surf);
    Daxpy (Data.N_C_Surf, 1/tau, Data.tmp_surf, Data.rhs_surf);
    MatAdd(Data.A_Surf, Data.M_Surf, 1/tau);
    
    DirectSolver(Data.A_Surf, Data.rhs_surf, Data.sol_surf);
    
//     InterpolateSurface(Data.c_surf, Data.Coll_P0, Data.CellNumbers_P0, Data.JointNumbers_P0,
// 		     InitialC_Surf);
    
    WriteSolution(VTK, img, Data.Output_Surf);
    DumpIsosurface(Inter, img, Data.c_surf,
		  Data.Coll_P0, Data.N_SurfaceJoints_P0, Data.CellNumbers_P0,
		  Data.JointNumbers_P0);
    ++img;
    
    mass = CalculateSurfaceMass(Data.c_surf, Data.concen_space, Data.Coll_P0, Data.CellNumbers_P0,
			      Data.JointNumbers_P0);
			      
    CalcSurfaceErrors ( Data.c_surf, Data.concen_space, Data.CellNumbers_P0, Data.JointNumbers_P0,
			ExactC_Surf, errors);
			
    error_loc = errors[0];
    PRINTVAR(error_loc);
    
    if ( error_L2 < error_loc ) error_L2 = error_loc;
			      
    dat << TDatabase::TimeDB->CURRENTTIME << " " << mass << endl;
    PRINTVAR(mass);
  }
  
  massloss = (massO - mass) / massO;
  PRINTVAR(massloss);
  PRINTVAR(error_L2);
  
  dat.close();
  
  return 0;
}

void MUTABLE_DATA::Create(TDomain *Domain)
{
  double hmin, hmax;
  char Name[] = "name";
  char Desc[] = "desc";
  char CSName[] = "c_surf";
  char GName[] = "g";
  char UName[] = "u";
  
  /// create all grid and collections
  SplitGrid_TwoPhase(Domain, Coll, Coll_P0, Coll_P1, phase_marker,
		     GlobalCellNo_P0, GlobalCellNo_P1,
		     N_SurfaceJoints_P0, CellNumbers_P0, JointNumbers_P0,
		     N_SurfaceJoints_P1, CellNumbers_P1, JointNumbers_P1);
		     
  SetInterfaceJoints(Coll, Coll_P1, N_SurfaceJoints_P1, CellNumbers_P1, JointNumbers_P1);
  
  MakeSurfaceGrid(Coll_P0, N_SurfaceJoints_P0, CellNumbers_P0, JointNumbers_P0, SurfColl);
  SurfColl->GetHminHmax(&hmin, &hmax);
//   CheckNeigbourhood(SurfColl);
//   TDatabase::TimeDB->TIMESTEPLENGTH = hmax*hmax;
  
  PRINTVAR(hmax);
  
  /// create fespaces
  velocity_space = new TFESpace3D (Coll, Name, Desc, BoundCondition, 1);
  concen_space   = new TFESpace3D (Coll_P0, Name, Desc, BoundCondition,
				      TDatabase::ParamDB->ANSATZ_ORDER);
  concen_space_surf = new TFESpace2D (SurfColl, Name, Desc, SurfaceBoundCondition,
				      TDatabase::ParamDB->ANSATZ_ORDER, NULL);
				      
  grid_space = new TFESpace3D (Coll, Name, Desc, GridBoundCond, 1);
  grid_space_p0 = new TFESpace3D (Coll_P0, Name, Desc, GridBoundCond, 1);
  grid_space_p1 = new TFESpace3D (Coll_P1, Name, Desc, GridBoundCond, 1);
				      
  /// create arrays
  N_U = velocity_space->GetN_DegreesOfFreedom();
  N_C_Surf = concen_space_surf->GetN_DegreesOfFreedom();
  N_Grid = grid_space->GetN_DegreesOfFreedom();
  N_Grid_P0 = grid_space_p0->GetN_DegreesOfFreedom();
  N_Grid_P1 = grid_space_p1->GetN_DegreesOfFreedom();
  
  sol = new double [3*N_U];
  sol_surf = new double [N_C_Surf];
  rhs_surf = new double [N_C_Surf]; 
  tmp_surf = new double [N_C_Surf];
  grid_velo = new double [3*N_Grid];
  grid_pos = new double [3*N_Grid];
  grid_velo_p0 = new double [3*N_Grid_P0];
  grid_velo_p1 = new double [3*N_Grid_P1];
  grid_rhs_p0 = new double [3*N_Grid_P0];
  grid_rhs_p1 = new double [3*N_Grid_P1];
  
  Print(N_U);
  Print(N_C_Surf);
  Print(N_Grid);
  Print(N_Grid_P0);
  Print(N_Grid_P1);
  
  /// create fe functions
  c_surf = new TFEFunction2D (concen_space_surf, CSName, Desc, sol_surf, N_C_Surf);
  
  u = new TFEVectFunct3D (velocity_space, UName, Desc, sol, N_U, 3);
  ux = u->GetComponent(0);
  uy = u->GetComponent(1);
  uz = u->GetComponent(2);
  
  g = new TFEVectFunct3D (grid_space, GName, Desc, grid_velo, N_Grid, 3);
  grid = new TFEVectFunct3D (grid_space, Name, Desc, grid_pos, N_Grid, 3);
  
  /// create structures
  sqstructure_surf = new TSquareStructure2D (concen_space_surf);
  sqstructure_surf->Sort();
  sqstruct_grid_p0 = new TSquareStructure3D (grid_space_p0);
  sqstruct_grid_p0->Sort();
  sqstruct_grid_p1 = new TSquareStructure3D (grid_space_p1);
  sqstruct_grid_p1->Sort();
  
  /// create matrices;
  A_Surf = new TSquareMatrix2D (sqstructure_surf);
  M_Surf = new TSquareMatrix2D (sqstructure_surf);
  
  G11_P0 = new TSquareMatrix3D (sqstruct_grid_p0);
  G12_P0 = new TSquareMatrix3D (sqstruct_grid_p0);
  G13_P0 = new TSquareMatrix3D (sqstruct_grid_p0);
  G21_P0 = new TSquareMatrix3D (sqstruct_grid_p0);
  G22_P0 = new TSquareMatrix3D (sqstruct_grid_p0);
  G23_P0 = new TSquareMatrix3D (sqstruct_grid_p0);
  G31_P0 = new TSquareMatrix3D (sqstruct_grid_p0);
  G32_P0 = new TSquareMatrix3D (sqstruct_grid_p0);
  G33_P0 = new TSquareMatrix3D (sqstruct_grid_p0);
  
  G11_P1 = new TSquareMatrix3D (sqstruct_grid_p1);
  G12_P1 = new TSquareMatrix3D (sqstruct_grid_p1);
  G13_P1 = new TSquareMatrix3D (sqstruct_grid_p1);
  G21_P1 = new TSquareMatrix3D (sqstruct_grid_p1);
  G22_P1 = new TSquareMatrix3D (sqstruct_grid_p1);
  G23_P1 = new TSquareMatrix3D (sqstruct_grid_p1);
  G31_P1 = new TSquareMatrix3D (sqstruct_grid_p1);
  G32_P1 = new TSquareMatrix3D (sqstruct_grid_p1);
  G33_P1 = new TSquareMatrix3D (sqstruct_grid_p1);
  
  /// create aux objects
  TFEFunction3D *FEFCT[3];
  TFEVectFunct3D *FEVCT[1];
  
  FEFCT[0] = ux;
  FEFCT[1] = uy;
  FEFCT[2] = uz;
  FEVCT[0] = u;
  
  AuxParam_Surf = new TAuxParam2D3D (TCD3D_Surf::N_FEFct, FEFCT, TCD3D_Surf::N_FEVectFct, FEVCT,
				     TCD3D_Surf::N_Parameters, TCD3D_Surf::FctIndex,
				     TCD3D_Surf::ParamDerivatives, TCD3D_Surf::IsGlobal);
				     
  /// create output objects
  Output_Surf = new TOutput2D (1, 1, 0, 0, NULL);
  Output = new TOutput3D (1, 3, 1, 0, Domain);
  
  Output_Surf->AddFEFunction(c_surf);
  Output->AddFEVectFunct(u);
  
  /// init grid
  grid->GridToData();
}

void WriteSolution(char *basename, int number, TOutput2D *out)
{
  std::ostringstream os;
  
  os << basename << ".";
  os.width(5);
  os.fill('0');
  os << number << ".vtk" << ends;
  
  out->WriteVtk(os.str().c_str());
}

void CheckNeigbourhood(TCollection *SurfColl)
{
  int N_Cells = SurfColl->GetN_Cells();
  TBaseCell **Cells, *Cell, *Neib;
  TJoint *Joint;
  TCollection *NeibColl;
  
  for (int i=0;i<N_Cells;++i)
  {
    Cell = SurfColl->GetCell(i);
    std::ostringstream os;
  
    os << "neighbour.";
    os.width(5);
    os.fill('0');
    os << i << ".vtk" << ends;
    
    Cells = new TBaseCell* [4];
    Cells[0] = Cell;
    
    for (int j=0;j<3;++j)
    {
      Joint = Cell->GetJoint(j);
      
      Neib = Joint->GetNeighbour(Cell);
      Cells[j+1] = Neib;
    }
    
    NeibColl = new TCollection(4, Cells);
    TOutput2D out (0,0,0,0, NULL, NeibColl);
    out.WriteVtk(os.str().c_str());
    
    delete NeibColl;
  } // end for i
}

