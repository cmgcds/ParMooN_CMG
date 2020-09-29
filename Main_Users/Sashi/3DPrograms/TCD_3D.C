// =======================================================================
//
// Purpose:     main program
//
// Author:      Volker John
//
// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <FESpace3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <AuxParam3D.h>
#include <QuadAffin.h>
#include <Solver.h>
#include <Assemble3D.h>
#include <Output3D.h>
#include <DiscreteForm3D.h>
#include <LinAlg.h>
#include <Collection.h>
#include <Upwind3D.h>
#include <TCD3D.h>
#include <FEM_TVD_FCT.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <malloc.h>

#include <ItMethod.h>
#include <FixedPointIte.h>
#include <MultiGridScaIte.h>
#include <FgmresIte.h>

#include <MultiGrid3D.h>
#include <MGLevel3D.h>

#include <MainUtilities.h>
#include <TimeUtilities.h>

#define AMG 0
#define GMG 1
double bound = 0;
// =======================================================================
// include current example
// =======================================================================

// #include "../Examples/TCD_3D/test0.h"
// #include "../Examples/TCD_3D/test2.h"
// #include "../Examples/TCD_3D/test3.h"
// #include "../Examples/TCD_3D/Linear.h"
// #include "../Examples/TCD_3D/non_linear_in_t.h"
// #include "../Examples/TCD_3D/bubble.h"
// #include "../Examples/TCD_3D/test2mod.h"
// #include "../Examples/TCD_3D/brenn_cd.h"
// #include "../Examples/TCD_3D/brennStutzenCD.h"
// #include "../Examples/TCD_3D/brennStutzenCD_P.h"
//#include "../Examples/TCD_3D/brennStutzenCD_Real.h"
//#include "../Examples/TCD_3D/bulk_compare.h"
//#include "../Examples/TCD_3D/Bail3D.h"
//#include "../Examples/TCD_3D/TEST_RFB.h"
#include "../Examples/TCD_3D/Sin4.h"



// ======================================================================
// utilities for main program
// ======================================================================

int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase3D *FEDatabase = new TFEDatabase3D();
  TCollection *coll, *mortarcoll = NULL;
  TBaseCell *cell;
  TFESpace3D *concentration_space, *pressure_space, *velocity_space; 
  TFESpace3D **ConcentrationSpaces, **VelocitySpaces;
  TOutput3D *Output;

  double *B, *rhs, *sol, *oldsol, tol, tolerance, *defect, *startsol, *frac_step_sol;
  double *app, *oldrhs, *itmethod_sol, *itmethod_rhs, *current_sol, *current_B;
  double *sol_velo, *oldsol_nlinite, damp_nlinite;
  int i,j,k,l,m,n, N_, Len, low, only_first_time;
  int N_Rows, N_Columns, N_Unknowns, N_Unknowns_Velo, N_P,sold_parameter_type;
  double  *lump_mass, *matrix_D_Entries, *oldrhs_fem_fct0, *tilde_u, *oldrhs_fem_fct1;
  double *l2u1, *l2u2, *h1u1, *h1u2;
  double *l2p, *h1p, *sd, *l_inf, vect[3], exactvect[3];
  int which, *permutation;
  double DiffL2, DiffH1, t;
  char *PRM, *GEO;
  int LEVELS, BASELEVEL, ORDER, order;
  int ret, pde;
  double negPower;
  double x,y,z,max,min,sum;
  double RE_NR;
  double tau1, tau2;
  double errors[7], p1, p2;
  double t1, t2, res, res2, oldres, solver_time, solver_time_curr, residual, oldresidual;
  double total_time, t3, t4, ass_time = 0.0;
  double impuls_residual,limit,linredfac;
  int N_LinIter, N_LinIterCurr, N_LinIterCurrIte, N_SubSteps, N_Active, n_aux;
  double gamma, tau, oldtau;
  int *RowPtr;

  std::ostringstream os;
  char *PsBaseName, *GrapeBaseName, *GmvBaseName, *ReadGrapeBaseName;
  char *VtkBaseName;

  TFEFunction3D *conc, **velo1, **velo2, **velo3, *fefct[3];
  TFEFunction3D **SolArray, **AuxFEFunctArray, *pressure;
  TFEFunction3D *oldconc, **OldSolArray;
  TFEVectFunct3D **velocity, **AuxFEVectFunctArray;
  TFESpace3D *fesp[2], *ferhs[1];
  double delta, end_time;

  TAuxParam3D *aux;
  TSquareStructure3D *sqstructureA;
  TSquareMatrix3D *sqmatrixA, *SQMATRICES[3];
  TSquareMatrix3D *sqmatrixM , *sqmatrixK, *sqmatrixS;
  TSquareMatrix3D **MatricesA, **MatricesM, **MatricesK, **MatricesS;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  MatVecProc *MatVect;
  DefectProc *Defect;
 
  double **RhsArray;

  TMGLevel3D *MGLevel;
  TMultiGrid3D *MG;

  double *RHSs[3];
  int *N_Uarray;

  // discrete forms have to be created
  TDiscreteForm3D *DiscreteForm;

  TDiscreteForm3D *DiscreteFormMatrixMRhs;
  TDiscreteForm3D *DiscreteFormMatrixMRhs_SUPG;

  TDiscreteForm3D *DiscreteFormMatrixARhs;
  TDiscreteForm3D *DiscreteFormMatricesAKRhs_SUPG;
  TDiscreteForm3D *DiscreteFormMatricesAKRhs_SOLD;
  TDiscreteForm3D *DiscreteFormRhs;
  TDiscreteForm3D *DiscreteFormRhs_SUPG;
  TDiscreteForm3D *DiscreteFormMatrixAUpwindRhs;

  int N_SquareMatrices, N_Rhs, N_FESpaces;
  
  BoundCondFunct3D *BoundaryConditions[1];
  BoundValueFunct3D *BoundValues[1];
  CoeffFct3D *Coefficients[1];

  TItMethod *itmethod, *prec;
  int Max_It, FirstSolve;
  double omega, alpha, alpha_fine, divergence;
  int N_Paramters=1, methods, time_discs;
  double Parameters[2], hmin, hmax;

  int N_UConv, level_down, ii, fully_implicit = 0;
  int mg_level,mg_type,CurrentDiscType, last_sq, step_length;
  int very_first_time=0, zerostart, comp_vort, picture;
  int first_matrix_assemble =1 ;
  int *neum_to_diri, N_neum_to_diri = 0;
  double *neum_to_diri_x, *neum_to_diri_y, *neum_to_diri_z;

  char ReadinDat[] = "readin.dat";
  char UString[] = "u";
  char PString[] = "p";
  char NameString[] = "name";
  char MMString[] = "Mass matrix";

#ifdef __FUEL_CELL__
  int N_OutFlowCells, *OutFlowFaces, *OutFlowCells;
  double outflow, area;
#endif // __FUEL_CELL__

#ifdef __STUTZEN__
  TBaseCell **Cells;
  int N_CellsOld, N_CellsNew;
  int N_Vertices;
  double xm, ym, zm;
  double xp, yp, zp;
  double r1, r2;
#endif

  os << " ";

  total_time = GetTime();
//======================================================================
// read parameter file
//======================================================================
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

#ifdef __STUTZEN__
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 1234;
#endif

  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();
//======================================================================
// copy read parameters into local variables
//======================================================================
 /* if( (TDatabase::ParamDB->DISCTYPE==2) )
  {
    OutPut("SDFEM does not work!" << endl);
    Error("SDFEM does not work!" << endl);
    exit(4711);
  }*/
 if( (TDatabase::ParamDB->DISCTYPE==5) )
  {
    OutPut("DISCTYPE 5 NOT IMPLEMENTED!" << endl);
    Error("DISCTYPE 5 NOT IMPLEMENTED!" << endl);
    exit(4711);
  }

   Coefficients[0] = BilinearCoeffs;
  if(TDatabase::ParamDB->DISCTYPE == GALERKIN)
  {
      if (TDatabase::ParamDB->SOLD_TYPE)
      {
	  TDatabase::ParamDB->SOLD_TYPE = 0;
	  OutPut("SOLDTYPE set to 0 !!!" << endl);
      }
  }
  // check parameter consistency, set internal parameters
  if (TDatabase::ParamDB->SOLD_TYPE)
  {
      sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;

    if(TDatabase::ParamDB->DISCTYPE != SDFEM)
    {
      TDatabase::ParamDB->DISCTYPE = SDFEM;
      OutPut("DISCTYPE set to " <<  TDatabase::ParamDB->DISCTYPE  << endl);
    }
    if ((TDatabase::ParamDB->SOLD_PARAMETER_TYPE!=JSW87)&&
	(TDatabase::ParamDB->SOLD_PARAMETER_TYPE!=JSW87_1))
      TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME = 0;
    else
      TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME = 1;
  }
 if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE==FEM_FCT)
      TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME = 0;

  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE==FEM_FCT_LIN)
  {
      TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME = 1;
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE=FEM_FCT;
      sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;
  }

  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;

  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  GrapeBaseName = TDatabase::ParamDB->GRAPEBASENAME;
  GmvBaseName = TDatabase::ParamDB->GMVBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  ReadGrapeBaseName = TDatabase::ParamDB->READGRAPEBASENAME;

  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SADDLE;
  if (mg_type)
    mg_level = 1;
  else
    mg_level = 0;

  LEVELS = TDatabase::ParamDB->LEVELS;
  BASELEVEL = TDatabase::ParamDB->UNIFORM_STEPS;
  l2u1 = new double[LEVELS+1];
  l2u2 = new double[LEVELS+1];
  h1u1 = new double[LEVELS+1];
  h1u2 = new double[LEVELS+1];
  sd = new double[LEVELS+1];
  l_inf = new double[LEVELS+1];

  SolArray = new TFEFunction3D*[LEVELS+1];
  OldSolArray = new TFEFunction3D*[LEVELS+1];
  RhsArray = new double* [LEVELS+1];
  N_Uarray = new int[LEVELS+1];
  velocity = new TFEVectFunct3D*[LEVELS+1];
  velo1 = new TFEFunction3D*[LEVELS+1];
  velo2 = new TFEFunction3D*[LEVELS+1];
  velo3 = new TFEFunction3D*[LEVELS+1];

  ConcentrationSpaces = new TFESpace3D*[LEVELS+1];
  VelocitySpaces = new TFESpace3D*[LEVELS+1];

  MatricesA = new TSquareMatrix3D*[LEVELS+1];
  MatricesM = new TSquareMatrix3D*[LEVELS+1];
  MatricesK = new TSquareMatrix3D*[LEVELS+1];
  // array which points to the stabilization  matrices (sold) on the
  // different levels of the multigrid
  // it is actually used only on the finest level
  MatricesS = new TSquareMatrix3D*[LEVELS+1];


  MatVect = MatVect_Scalar;
  Defect = Defect_Scalar;

//======================================================================
// initialize discrete forms
//======================================================================

// discrete form for assembling mass matrix
   DiscreteFormMatrixMRhs = new TDiscreteForm3D
    (MMString, MMString, N_Terms_MatrixMRhs,
     Derivatives_MatrixMRhs,
     SpacesNumbers_MatrixMRhs, N_Matrices_MatrixMRhs, N_Rhs_MatrixMRhs,
     RowSpace_MatrixMRhs, ColumnSpace_MatrixMRhs, RhsSpace_MatrixMRhs,
     MatrixMRhsAssemble, BilinearCoeffs, NULL);

   DiscreteFormMatrixMRhs_SUPG = new TDiscreteForm3D
    (MMString, MMString, N_Terms_MatrixMRhs_SUPG,
     Derivatives_MatrixMRhs_SUPG,
     SpacesNumbers_MatrixMRhs_SUPG, N_Matrices_MatrixMRhs_SUPG,
     N_Rhs_MatrixMRhs_SUPG,
     RowSpace_MatrixMRhs_SUPG, ColumnSpace_MatrixMRhs_SUPG,
     RhsSpace_MatrixMRhs_SUPG,
     MatrixMRhsAssemble_SUPG, BilinearCoeffs, NULL);

    DiscreteFormMatrixARhs = new TDiscreteForm3D
    (MMString, MMString, N_Terms_MatrixARhs,
     Derivatives_MatrixARhs,
     SpacesNumbers_MatrixARhs, N_Matrices_MatrixARhs,
     N_Rhs_MatrixARhs,
     RowSpace_MatrixARhs, ColumnSpace_MatrixARhs,
     RhsSpace_MatrixARhs,
     MatrixARhsAssemble, BilinearCoeffs, NULL);

    DiscreteFormMatricesAKRhs_SUPG = new TDiscreteForm3D
    (MMString, MMString, N_Terms_MatricesAKRhs_SUPG, Derivatives_MatricesAKRhs_SUPG,
     SpacesNumbers_MatricesAKRhs_SUPG, N_Matrices_MatricesAKRhs_SUPG, N_Rhs_MatricesAKRhs_SUPG,
     RowSpace_MatricesAKRhs_SUPG, ColumnSpace_MatricesAKRhs_SUPG, RhsSpace_MatricesAKRhs_SUPG,
     MatricesAKRhsAssemble_SUPG, BilinearCoeffs, NULL);

    DiscreteFormRhs_SUPG = new TDiscreteForm3D
    (MMString, MMString, N_Terms_Rhs_SUPG, Derivatives_Rhs_SUPG,
     SpacesNumbers_Rhs_SUPG, N_Matrices_Rhs_SUPG, N_Rhs_Rhs_SUPG,
     RowSpace_Rhs_SUPG, ColumnSpace_Rhs_SUPG, RhsSpace_Rhs_SUPG,
     RhsAssemble_SUPG, BilinearCoeffs, NULL);

    DiscreteFormRhs = new TDiscreteForm3D
    (MMString, MMString, N_Terms_Rhs, Derivatives_Rhs,
     SpacesNumbers_Rhs, N_Matrices_Rhs, N_Rhs_Rhs,
     RowSpace_Rhs, ColumnSpace_Rhs, RhsSpace_Rhs,
     RhsAssemble, BilinearCoeffs, NULL);

    DiscreteFormMatrixAUpwindRhs = new TDiscreteForm3D
    (MMString, MMString, N_Terms_MatrixARhs, Derivatives_MatrixARhs,
     SpacesNumbers_MatrixARhs, N_Matrices_MatrixARhs, N_Rhs_MatrixARhs,
     RowSpace_MatrixARhs, ColumnSpace_MatrixARhs, RhsSpace_MatrixARhs,
     MatrixAUpwindRhsAssemble, BilinearCoeffs, NULL);

    
    // discrete form for assembling stiffness matrix, stabilization matrix and rhs (SDFEM)
    DiscreteFormMatricesAKRhs_SOLD = new TDiscreteForm3D
	(MMString, MMString, N_Terms_MatricesAKRhs_SUPG, Derivatives_MatricesAKRhs_SUPG,
	 SpacesNumbers_MatricesAKRhs_SUPG, N_Matrices_MatricesAKRhs_SOLD, N_Rhs_MatricesAKRhs_SUPG,
	 RowSpace_MatricesAKRhs_SOLD, ColumnSpace_MatricesAKRhs_SOLD, RhsSpace_MatricesAKRhs_SUPG,
	 MatricesAKRhsAssemble_SUPG, BilinearCoeffs, NULL);
    

//======================================================================
// read boundary parameterization and initialize coarse grid
//======================================================================
  Domain->Init(PRM, GEO);

  Domain->ConvertQuadToTri(TDatabase::ParamDB->CONVERT_QUAD_TO_TRI);

  for(i=0;i<TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR;i++)
    Domain->RegRefineAll();

  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR;
  omega = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SCALAR;
  alpha = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SCALAR;
  alpha_fine = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_FINE_SCALAR;
  divergence = TDatabase::ParamDB->SC_DIV_FACTOR;
  
  Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
  Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;

  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  SetTimeDiscParameters();

  if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    i=1;
    MG = new TMultiGrid3D(i, N_Paramters, Parameters);
  }
 
  mg_level = LEVELS+mg_level;
  damp_nlinite = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SCALAR;

  t3 = GetTime();
  total_time = t3 - total_time;
//======================================================================
// loop over all levels
// do all refinements
// build all fe spaces
// load velocity on the finest level
//======================================================================
  for(i=0;i<mg_level;i++)
  {
    if (i<LEVELS)
    {
      OutPut("*******************************************************" << endl);
      OutPut("******           GEOMETRY  LEVEL ");
      OutPut(i << "              *******" << endl);
    }
    else
    {
      OutPut("*******************************************************" << endl);
      OutPut("******           GEOMETRY  LEVEL ");
      OutPut(i-1 << "              *******" << endl);
    }
    OutPut("******           MULTIGRID LEVEL ");
    OutPut(i << "              *******" << endl);
    OutPut("*******************************************************" << endl);
    OutPut("memory before: " << setw(10) << GetMemory() << endl);

    // refine domain regularly
    if(i && (i<LEVELS)) Domain->RegRefineAll();

    // create collection of mesh cells
    coll=Domain->GetCollection(It_Finest, 0);
    
    // if multiple discretization multilevel method is used
    // get space for low order disc on finest geo grid
    if ((mg_type==1)&&(i<mg_level-1))
    {
      // nonconforming velocity space
      concentration_space = new TFESpace3D(coll,NameString,UString,BoundCondition,-1);
    }
    // standard multigrid or finest level
    // get fe space of high order disc on finest geo grid
    else
    { 
      ORDER  = TDatabase::ParamDB->ANSATZ_ORDER; 
      concentration_space = new TFESpace3D(coll, NameString, 
                                      UString, BoundCondition, ORDER);
    }
    
    // array of the fe spaces
    ConcentrationSpaces[i] = concentration_space;
    N_Unknowns = concentration_space->GetN_DegreesOfFreedom();
    N_Uarray[i] = N_Unknowns;

    // active dof (i.e. dof without Dirichlet dofs)
    N_Active = concentration_space->GetActiveBound();
    OutPut("degrees of freedom: "<< setw(10) << N_Unknowns << endl);

    // build matrices 
    // first build matrix structure
    sqstructureA = new TSquareStructure3D(concentration_space);
    sqstructureA->Sort();
    
    // two matrices used
    // A contains the non time dependent part of the discretization
    sqmatrixA = new TSquareMatrix3D(sqstructureA);
    MatricesA[i] = sqmatrixA;

    // M is the mass matrix
    // the iterative solver uses M
    sqmatrixM = new TSquareMatrix3D(sqstructureA);
    MatricesM[i] = sqmatrixM;

    if(TDatabase::ParamDB->DISCTYPE == SDFEM)
    {
      // stabilisation matrix K
      sqmatrixK = new TSquareMatrix3D(sqstructureA);
      MatricesK[i] = sqmatrixK;
       if(TDatabase::ParamDB->SOLD_TYPE)
      {
	  // stabilisation matrix S
	  sqmatrixS = new TSquareMatrix3D(sqstructureA);
	  MatricesS[i] = sqmatrixS;
      }
   }

    // allocate rhs 
    rhs = new double[N_Unknowns];
    memset(rhs, 0, N_Unknowns*SizeOfDouble);  
    RhsArray[i] = rhs;

    // allocate array for solution 
    sol = new double[N_Unknowns];
    memset(sol, 0, N_Unknowns*SizeOfDouble);
    current_sol = sol;

    // allocate fe function on the current level
    conc = new TFEFunction3D(concentration_space, UString, UString, sol, N_Unknowns);
    SolArray[i] = conc;

    if (i==mg_level -1)
    {
      oldsol = new double[N_Unknowns];
      oldrhs = new double[N_Unknowns];
      if (!(TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME))
	  oldsol_nlinite = new double[N_Unknowns];	  
    }
    // SOLD schemes
    if ((TDatabase::ParamDB->SOLD_TYPE)&&(i==mg_level-1))
    {
      oldconc =  new TFEFunction3D(concentration_space, UString, UString, oldsol, N_Unknowns);
      OldSolArray[i] = oldconc;
    }

    // allocate array on which the time scheme works
    B = new double [N_Unknowns];
    current_B = B;
    
 // allocate array for lumped mass matrix
    if ((TDatabase::ParamDB->SOLD_PARAMETER_TYPE== FEM_FCT)&&
      (i==mg_level-1))
    { 
        
        
	oldrhs_fem_fct0 = new double[N_Unknowns];
	oldrhs_fem_fct1 = new double[N_Unknowns];
      lump_mass = new double [N_Unknowns];
      matrix_D_Entries = new double[sqmatrixA->GetN_Entries()];
      // matrix K for copy of mass matrix
      sqmatrixK = new TSquareMatrix3D(sqstructureA);
      MatricesK[i] = sqmatrixK;
    }
    // prepare multigrid method
    // if geometric multigrid method
    if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
    {
       // determine number of auxiliary arrays
       if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR)
           || (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR))
	  n_aux=4;
       else
	  n_aux=2;
       // allocate new level and add it
       MGLevel = new TMGLevel3D(i, sqmatrixM, current_B, current_sol,  n_aux, NULL);
       MG->AddLevel(MGLevel);
    }

    // define everything for velocity (convection)
    // array which contains the values
    if ((mg_type==1)&&(i<mg_level-1))
    {
       velocity_space = new TFESpace3D(coll, NameString, UString, BoundCondition,
        Non_USpace, 1);
    }
    else
    {
       ORDER  = TDatabase::ParamDB->ANSATZ_ORDER; 
       velocity_space = new TFESpace3D(coll, NameString, 
                                       UString, BoundCondition, ContP_USpace, ORDER);
    }

    VelocitySpaces[i] = velocity_space; 

    N_Unknowns_Velo = velocity_space->GetN_DegreesOfFreedom();
    sol_velo = new double[3*N_Unknowns_Velo];
    memset(sol_velo, 0, 3*N_Unknowns_Velo*SizeOfDouble);  

    // vector fe function
    velocity[i] = new TFEVectFunct3D(velocity_space, UString, UString, sol_velo, N_Unknowns, 3);
    // individual components of velocity
    velo1[i] = velocity[i]->GetComponent(0);
    velo2[i] = velocity[i]->GetComponent(1);
    velo3[i] = velocity[i]->GetComponent(2);

    // read velocity on finest level from grape file
    if ((i==mg_level-1)&&(TDatabase::ParamDB->READ_GRAPE_FILE))
    {
      // need temporarily the pressure space
      if(TDatabase::ParamDB->ANSATZ_ORDER == -1)
        // piecewise constant for nonconforming velocity
        pressure_space = new TFESpace3D(coll, NameString, PString, 
                                        BoundCondition, 
                                        DiscP_PSpace, 0);
      else
        // implicitely assumed that pressure space is discontinuous linear
        pressure_space = new TFESpace3D(coll, NameString, PString,
                                        BoundCondition, 
                                        DiscP_PSpace, 1);

      N_P = pressure_space->GetN_DegreesOfFreedom();
      pressure = new TFEFunction3D(pressure_space, PString, PString, sol, N_P);

      // do read pressure as scalar fe functions
      AuxFEFunctArray = new TFEFunction3D*[1];
      AuxFEFunctArray[0] = pressure;

      // one vector fe function -> velocity
      AuxFEVectFunctArray = new TFEVectFunct3D*[1];
      AuxFEVectFunctArray[0] = velocity[i];

      // read velocity
      ReadGrapeFile3D(ReadGrapeBaseName, 1, 1, AuxFEFunctArray, AuxFEVectFunctArray);

      // release memory
      delete pressure;
      delete pressure_space;

      // norm of convection
      OutPut("norm of convection " << Dnorm(N_Unknowns,sol_velo) << " ");
      OutPut( Dnorm(N_Unknowns,sol_velo+N_Unknowns) << " " <<  Dnorm(N_Unknowns,sol_velo+2*N_Unknowns ));
      OutPut("  total " <<  Dnorm(3*N_Unknowns,sol_velo) << endl);            
    }
    
    // prepare output, only the concentration will be saved
    if (i==mg_level-1)
    {
       Output = new TOutput3D(1, 1, 0, 1, Domain);
       Output->AddFEFunction(conc);
       os.seekp(std::ios::beg);
       Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());
    }

    // interpolate initial condition
    conc->Interpolate(InitialCondition);
    OutPut(i << " time: " << TDatabase::TimeDB->CURRENTTIME << endl);
  } // endfor i

//======================================================================
// loop over all levels
// restrict velocity
//======================================================================

  if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
  {  
  for(i = mg_level-1 ; i > 0;i--)
  {
     RestrictFunction(VelocitySpaces[i-1], VelocitySpaces[i],
                      velo1[i-1]->GetValues(),
                      velo1[i]->GetValues(),
                      MG->GetLevel(i-1)->GetAuxVector(0));
     RestrictFunction(VelocitySpaces[i-1], VelocitySpaces[i],
                      velo2[i-1]->GetValues(),
                      velo2[i]->GetValues(),
                      MG->GetLevel(i-1)->GetAuxVector(0));
     RestrictFunction(VelocitySpaces[i-1], VelocitySpaces[i],
                      velo3[i-1]->GetValues(),
                      velo3[i]->GetValues(),
                      MG->GetLevel(i-1)->GetAuxVector(0));
  } // endfor i
  }
//======================================================================
// all data are available for assembling matrices
// loop over all levels
//======================================================================
  for(i=0;i<mg_level;i++)
  {
    // set parameters
    N_Rhs = 1;
    N_FESpaces = 2;
    fesp[0] = ConcentrationSpaces[i];
    fesp[1] = VelocitySpaces[i];
    fefct[0] = velo1[i];
    fefct[1] = velo2[i];
    fefct[2] = velo3[i];

    aux =  new TAuxParam3D(TimeCDParamsVeloFieldN_FESpaces,
                           TimeCDParamsVeloFieldN_Fct, 
                           TimeCDParamsVeloFieldN_ParamFct, 
                           TimeCDParamsVeloFieldN_FEValues, 
                           fesp+1, fefct, 
                           TimeCDParamsVeloFieldFct, 
                           TimeCDParamsVeloFieldFEFctIndex,
                           TimeCDParamsVeloFieldFEMultiIndex, 
                           TimeCDParamsVeloFieldN_Params,
                           TimeCDParamsVeloFieldBeginParam); 

    //======================================================================
    // assembling of mass matrix and rhs
    //======================================================================
    // reset matrices
    N_SquareMatrices = 1;
    SQMATRICES[0] = MatricesM[i];
    SQMATRICES[0]->Reset();

    if(TDatabase::ParamDB->DISCTYPE == SDFEM)
      DiscreteForm = DiscreteFormMatrixMRhs_SUPG; 
    else
      DiscreteForm = DiscreteFormMatrixMRhs; 

    BoundaryConditions[0] =  BoundCondition;
    BoundValues[0] = BoundValue;
    if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE==FEM_FCT)
    {
	BoundValues[0] = BoundValue_FEM_FCT;
    }

    memset(RhsArray[i], 0, N_Uarray[i]*SizeOfDouble);  
    RHSs[0] = RhsArray[i];
    ferhs[0] = ConcentrationSpaces[i];

    Assemble3D(N_FESpaces, fesp, 
               N_SquareMatrices, SQMATRICES, 
               0, NULL, 
               N_Rhs, RHSs, ferhs,
               DiscreteForm, 
               BoundaryConditions, 
               BoundValues, 
               aux);
    delete aux;    

    if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE==FEM_FCT)
    {
	BoundValues[0] = BoundValue;
    }
    
    if ((TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_FCT)&&
	(i==mg_level-1))
    {
	lump_mass = new double [N_Unknowns];
	LumpMassMatrixToVector(MatricesM[i], lump_mass);
	tilde_u = new double [N_Unknowns];
	// save mass matrix in matricesK
	memcpy(MatricesK[i]->GetEntries(), MatricesM[i]->GetEntries(),
	       MatricesM[i]->GetN_Entries() * SizeOfDouble);
	// save rhs of previous time step
      memcpy(oldrhs_fem_fct0, rhs, N_Unknowns*SizeOfDouble);
#ifdef __BAIL_3D__
      CheckWrongNeumannNodes(coll, ConcentrationSpaces[i], N_neum_to_diri, neum_to_diri,
        neum_to_diri_x, neum_to_diri_y, neum_to_diri_z);
#endif
#ifdef __SIN4__
      CheckWrongNeumannNodes(coll, ConcentrationSpaces[i], N_neum_to_diri, neum_to_diri,
        neum_to_diri_x, neum_to_diri_y, neum_to_diri_z);
#endif
    
    }


} // endfor i

 
#ifdef __FUEL_CELL__
  // find all cells which have a face on the outflow area

  GetOutFlowCells(coll, N_OutFlowCells, OutFlowCells, OutFlowFaces);
#endif // __FUEL_CELL__


//======================================================================
// end of space cycle, finest grid reached
// everything happens now on the same grid
//======================================================================

 // save solution
  memcpy(oldsol,sol,N_Unknowns*SizeOfDouble);

 
  picture = 0;
  if(TDatabase::ParamDB->WRITE_GRAPE)
  {
    os.seekp(std::ios::beg);
    os << GrapeBaseName << 0 << ".dat" << ends;
    Output->WriteGrape(os.str().c_str());
    picture = 1;
  }


  if (TDatabase::ParamDB->WRITE_GMV)
  {

    os.seekp(std::ios::beg);
    os << GmvBaseName << 0 << ".gmv" << ends;
    Output->WriteGMV(os.str().c_str());
    picture = 1;
 
  }

  if(TDatabase::ParamDB->WRITE_VTK)
  {
      os.seekp(std::ios::beg);
      os << VtkBaseName << 0 << ".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
  }

  // allocate arrays for solver
  defect = new double[N_Unknowns];
  startsol = new double[N_Unknowns];
  frac_step_sol = new double[N_Unknowns];
  oldrhs =  new double[N_Unknowns];

  // number of active d.o.f.
  N_Active = ConcentrationSpaces[mg_level-1]->GetActiveBound();
  N_Unknowns = ConcentrationSpaces[mg_level-1]->GetN_DegreesOfFreedom();

  solver_time = 0.0;
  N_LinIter = 0;
  
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

  // initialize solver
  if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
  {
    switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
    {
      case 11:
        zerostart = 1;
        break;
      case 16:
        zerostart = 0;
        break;
    }
    switch (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR)
    {
      case 5:
        prec = new TMultiGridScaIte(MatVect, Defect, NULL,
                                    0, N_Unknowns, MG, zerostart);
        break;
      default:
        OutPut("Unknown preconditioner !!!" << endl);
       exit(4711);
    }
    switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
    {
       // fixed point iteration 
      case 11:
        itmethod = new TFixedPointIte(MatVect, Defect, prec,
                                      0, N_Unknowns, 1);
        if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
        {
          itmethod_sol = new double[N_Unknowns];
          itmethod_rhs = new double[N_Unknowns];              
          }  
        else
        {
          itmethod_sol = sol;
          itmethod_rhs = rhs;            
        }
        break;
      case 16:
         // FGMRES
        itmethod = new TFgmresIte(MatVect, Defect, prec,
                                  0, N_Unknowns, 1);
        if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
        {
          itmethod_sol = new double[N_Unknowns];
          itmethod_rhs = new double[N_Unknowns];              
        } 
        else
        {
          itmethod_sol = sol;
          itmethod_rhs = rhs;            
        }
        break;
      default:
        OutPut("Unknown solver !!!" << endl);
        exit(4711);
    }
  }

  coll->GetHminHmax(&hmin,&hmax);
  OutPut("h_min : " << hmin << " h_max : " << hmax << endl);

//======================================================================
// start of time cycle
// everything happens on the same grid
//======================================================================
  while(TDatabase::TimeDB->CURRENTTIME< end_time)
  {                              // time cycle
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
    
    for (methods=0;methods<time_discs;methods++)
    {
      if (time_discs==2)
      {
        if (methods==0) // fractional-step-theta-scheme
        {
          TDatabase::TimeDB->TIME_DISC = 3;
          memcpy(startsol,sol,N_Unknowns*SizeOfDouble); // save start sol
          memcpy(oldrhs,rhs,N_Unknowns*SizeOfDouble); // save rhs
        }
        else           // crank nicolson scheme
        {              // take solution of first scheme as initial iterate
          TDatabase::TimeDB->TIME_DISC = 2;
          TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->INTERNAL_STARTTIME;
          memcpy(frac_step_sol,sol,N_Unknowns*SizeOfDouble); // save solution of fract.step
          memcpy(sol,startsol,N_Unknowns*SizeOfDouble); // get former startsol
          memcpy(rhs,oldrhs,N_Unknowns*SizeOfDouble); // get old rhs
        }
        N_SubSteps = GetN_SubSteps();
      }

      for(l=0;l<N_SubSteps;l++)      // sub steps of fractional step theta
      {
        if (!very_first_time)
        {
          SetTimeDiscParameters();
        }
        if (m==1)
        {
          OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
          OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
          OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
          OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
        }
        
        tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
        if (!very_first_time)
          TDatabase::TimeDB->CURRENTTIME += tau;
         j=0; 
        OutPut(endl << "CURRENT TIME: ");
        OutPut(TDatabase::TimeDB->CURRENTTIME << endl);

	solver_time_curr = 0;
   
       N_LinIterCurr = 0;
       // start iteration for solving nonlinear problems
      while(1)
      {
        // working array for rhs is B, initialize B
        memset(B, 0, N_Unknowns*SizeOfDouble);
        // old rhs multiplied with current subtime step and theta3 on B
        Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, rhs, B);
       
        
	OutPut("z rhs1 " << rhs[29792] << endl);
 
        // assembling of A and rhs
        // for SDFEM: in addition stabilisation matrix K
        if ((TDatabase::ParamDB->INTERNAL_STEADY_STATE_MATRICES_OR_RHS==0) ||
            (first_matrix_assemble))
        {
          for(i=0;i<mg_level;i++)
          {  
            if ((TDatabase::ParamDB->SOLD_PARAMETER_TYPE==FEM_FCT) && (j>0))
		  break;
	    t1 = GetTime();
	    
            // set parameters
            N_Rhs = 1;
	    N_FESpaces = 1;
	    fesp[0] = ConcentrationSpaces[i];
	    if (TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
            {
		N_FESpaces = 2;
		fesp[1] = VelocitySpaces[i];
		fefct[0] = velo1[i];
		fefct[1] = velo2[i];
		fefct[2] = velo3[i];
          
		aux =  new TAuxParam3D(TimeCDParamsVeloFieldN_FESpaces,
				       TimeCDParamsVeloFieldN_Fct, 
				       TimeCDParamsVeloFieldN_ParamFct, 
				       TimeCDParamsVeloFieldN_FEValues, 
				       fesp+1, fefct, 
				       TimeCDParamsVeloFieldFct, 
				       TimeCDParamsVeloFieldFEFctIndex,
				       TimeCDParamsVeloFieldFEMultiIndex, 
				       TimeCDParamsVeloFieldN_Params,
				       TimeCDParamsVeloFieldBeginParam); 
            }
            // SOLD methods, only on finest grid
            // on coarser grids do SDFEM
            else if ((TDatabase::ParamDB->DISCTYPE == SDFEM) && 
		     (TDatabase::ParamDB->SOLD_TYPE)&&(i==mg_level-1))
            {
              TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE = 1;
              // current solution
              fefct[0] = SolArray[i];
              // old solution
              fefct[1] = OldSolArray[i];

              aux =  new TAuxParam3D(TimeCDParamsSOLDN_FESpaces,
                TimeCDParamsSOLDN_Fct,
                TimeCDParamsSOLDN_ParamFct,
                TimeCDParamsSOLDN_FEValues,
                fesp, fefct,
                TimeCDParamsSOLDFct,
                TimeCDParamsSOLDFEFctIndex,
                TimeCDParamsSOLDFEMultiIndex,
                TimeCDParamsSOLDN_Params,
                TimeCDParamsSOLDBeginParam);
            }
            else
            {
              aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
            }

            //======================================================================
            // assembling of mass matrix and rhs
            //======================================================================
            if ((mg_type==0) ||(i == mg_level-1))
            {
              if(TDatabase::ParamDB->DISCTYPE == SDFEM)
              {
                N_SquareMatrices = 2;
                SQMATRICES[0] = MatricesA[i];
                SQMATRICES[0]->Reset();
                SQMATRICES[1] = MatricesK[i];
                SQMATRICES[1]->Reset();
                DiscreteForm = DiscreteFormMatricesAKRhs_SUPG; 
		if (TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE)
		{
		    N_SquareMatrices = 3;
		    SQMATRICES[2] = MatricesS[i];
		    SQMATRICES[2]->Reset();		  
		    DiscreteForm = DiscreteFormMatricesAKRhs_SOLD;
		}
              }
              else
              {
                N_SquareMatrices = 1;
                SQMATRICES[0] = MatricesA[i];
                SQMATRICES[0]->Reset();
                if(TDatabase::ParamDB->DISCTYPE == UPWIND)
                  DiscreteForm = DiscreteFormMatrixAUpwindRhs;
                else
                  DiscreteForm = DiscreteFormMatrixARhs; 
              }
            }
            else
            {
              N_SquareMatrices = 1;
              SQMATRICES[0] = MatricesA[i];
              SQMATRICES[0]->Reset();
              DiscreteForm = DiscreteFormMatrixAUpwindRhs;              
            }
            
            BoundaryConditions[0] =  BoundCondition;
            BoundValues[0] = BoundValue;
 	    if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE==FEM_FCT)
	    {
		BoundValues[0] = BoundValue_FEM_FCT;
	    }
            memset(RhsArray[i], 0, N_Uarray[i]*SizeOfDouble);  
            RHSs[0] = RhsArray[i];
            ferhs[0] = ConcentrationSpaces[i];
            
            Assemble3D(N_FESpaces, fesp, 
                       N_SquareMatrices, SQMATRICES, 
                       0, NULL, 
                       N_Rhs, RHSs, ferhs,
                       DiscreteForm, 
                       BoundaryConditions, 
                       BoundValues, 
                       aux);
              // save rhs without Dirichlet values
	    if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE==FEM_FCT)
	    {
		BoundValues[0] = BoundValue;
		memcpy(oldrhs_fem_fct1, rhs, N_Unknowns*SizeOfDouble);
 	    }

            if (DiscreteForm == DiscreteFormMatrixAUpwindRhs)
            {
              UpwindForNavierStokes3D(SQMATRICES[0],velo1[i],velo2[i],  velo3[i]);  
            }
            if ((TDatabase::ParamDB->SOLD_TYPE)&&(i==mg_level-1))
            {
		// set parameter back to zero
		TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE = 0;
	    } // end SOLD type
            delete aux;   
	    t2 = GetTime();
	    ass_time += t2-t1;
	  } // endfor i
          first_matrix_assemble = 0;
	  OutPut("time for assembling " << t2-t1 << "s " << ass_time  << "s "<<  endl);
        }        
	OutPut("z rhs1 " << rhs[29792] << endl);
        //======================================================================
        // only rhs needs to be assembles  
        // only on finest level  
        //======================================================================
        /*if (TDatabase::ParamDB->INTERNAL_STEADY_STATE_MATRICES_OR_RHS==1) 
        {
          // set parameters
          N_Rhs = 1;
          N_SquareMatrices = 0;

          if(TDatabase::ParamDB->DISCTYPE == SDFEM)
          {
            N_FESpaces = 2;
            fesp[0] = ConcentrationSpaces[mg_level-1];
            fesp[1] = VelocitySpaces[mg_level-1];
            fefct[0] = velo1[mg_level-1];
            fefct[1] = velo2[mg_level-1];
            fefct[2] = velo3[mg_level-1];
            
            aux =  new TAuxParam3D(TimeCDParamsVeloFieldN_FESpaces,
                                   TimeCDParamsVeloFieldN_Fct, 
                                   TimeCDParamsVeloFieldN_ParamFct, 
                                   TimeCDParamsVeloFieldN_FEValues, 
                                   fesp+1, fefct, 
                                   TimeCDParamsVeloFieldFct, 
                                   TimeCDParamsVeloFieldFEFctIndex,
                                   TimeCDParamsVeloFieldFEMultiIndex, 
                                   TimeCDParamsVeloFieldN_Params,
                                   TimeCDParamsVeloFieldBeginParam); 
            DiscreteForm = DiscreteFormRhs_SUPG; 
          }
          else
          {
            N_FESpaces = 1;
            fesp[0] = ConcentrationSpaces[mg_level-1];
            
            aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
            DiscreteForm = DiscreteFormRhs; 
          }

          //======================================================================
          // assembling of rhs
          //======================================================================
           
          BoundaryConditions[0] =  BoundCondition;
          BoundValues[0] = BoundValue;
            
          memset(RhsArray[mg_level-1], 0, N_Uarray[mg_level-1]*SizeOfDouble);  
          RHSs[0] = RhsArray[mg_level-1];
          ferhs[0] = ConcentrationSpaces[mg_level-1];
            
          Assemble3D(N_FESpaces, fesp, 
                     0, NULL, 
                     0, NULL, 
                     N_Rhs, RHSs, ferhs,
                     DiscreteForm, 
                     BoundaryConditions, 
                     BoundValues, 
                     aux);
          delete aux;    
	  } */
   
        if (very_first_time==1)
        {
           very_first_time=0;
           l--;
           continue;
        }
	OutPut("z rhs1 " << rhs[29792] << endl);
        // add rhs from current sub time step to rhs array B
        Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, rhs, B);
          
        //======================================================================
        // manipulation of matrices due to current time discretization
        // the stiffness matrix is stored on M 
        //======================================================================
        oldtau = tau;
        if (!(TDatabase::ParamDB->SOLD_PARAMETER_TYPE== FEM_FCT))
	{
	    // update rhs by Laplacian and convective term from previous
          // time step
          // scaled by current sub time step length and theta2
          // currently : M := M + gamma A
          // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2)
          for(i=0;i<mg_level;i++)
          {
            MatAdd(MatricesM[i], MatricesA[i],
              -gamma - tau*TDatabase::TimeDB->THETA2);
          }
          // set current factor of steady state matrix
          gamma = -tau*TDatabase::TimeDB->THETA2;

          if(TDatabase::ParamDB->DISCTYPE == SDFEM)
          {
            for(i=0;i<mg_level;i++)
              MatAdd(MatricesM[i], MatricesK[i], 1);       
          }
          // defect = M * sol
          // B:= B + defect
          MatVectActive(MatricesM[mg_level-1], sol, defect);
          Daxpy(N_Active, 1, defect, B);

	  // right hand side needs to be computed only at the first loop 
	  if (j==0)
	      memcpy(oldrhs, B, N_Unknowns*SizeOfDouble);
	  else
	      memcpy(B, oldrhs, N_Unknowns*SizeOfDouble);

          // set Dirichlet values
          // RHSs[0] still available from assembling
          memcpy(B+N_Active, RHSs[0]+N_Active, (N_Unknowns-N_Active)*SizeOfDouble);
          // copy Dirichlet values from rhs into sol
          memcpy(sol+N_Active, RHSs[0]+N_Active, (N_Unknowns-N_Active)*SizeOfDouble);

          // M = M + (-gamma + tau*TDatabase::TimeDB->THETA1) A
          for(i=0;i<mg_level;i++)
          {
            MatAdd(MatricesM[i], MatricesA[i],
              -gamma + tau*TDatabase::TimeDB->THETA1);
          }
          // set current factor of steady state matrix
          gamma = tau*TDatabase::TimeDB->THETA1;
	  // SOLD
          if ((TDatabase::ParamDB->DISCTYPE == SDFEM)&&(TDatabase::ParamDB->SOLD_TYPE))
	      MatAdd(MatricesM[mg_level-1], MatricesS[mg_level-1],1); 
        }
        else
        {
          // FEM-FCT methods
          only_first_time = 1;
          if(j>0) 
	      only_first_time = 0;

	  // Here, already the new matrix A is used. This is slightly incorrect
	  // in the case of non-constant diffusion or convection. In this case,
	  // using the old matrix is correct. The new matrix should be used
	  // later in FEM_FCT_SystemMatrix(...)
          FEM_FCT_ForConvDiff(MatricesK[mg_level-1], MatricesA[mg_level-1],
			      N_Unknowns, N_Active, 
			      lump_mass, matrix_D_Entries, 
			      sol, oldsol, 
			      B, RhsArray[mg_level-1], oldrhs_fem_fct0, tilde_u,
			      N_neum_to_diri, neum_to_diri, 
			      neum_to_diri_x, neum_to_diri_y, neum_to_diri_z,
			      only_first_time,
			      BoundValue,NULL);
	
          SQMATRICES[0] = MatricesM[mg_level-1];
	  // only first iteration
	  if (j==0)
	  {
	      MatricesM[mg_level-1]->Reset();
	      // system matrix for FEM-FCT   M_lump + theta1*tau*A
	      // A = Galerkin + D
	      FEM_FCT_SystemMatrix(MatricesM[mg_level-1], MatricesA[mg_level-1],
					  lump_mass, N_Unknowns);
	  }

	  OutPut("entries " << Ddot(MatricesM[mg_level-1]->GetN_Entries(),MatricesM[mg_level-1]->GetEntries(),MatricesM[mg_level-1]->GetEntries()) << 
		 " lump " << Ddot(N_Unknowns,lump_mass,lump_mass) << endl);
 	  // set Diriclet nodes
	  if (N_neum_to_diri)
	      SetDirichletNodesFromNeumannNodes(SQMATRICES, B, sol, 
						N_neum_to_diri, neum_to_diri,
						neum_to_diri_x, neum_to_diri_y, neum_to_diri_z,
						BoundValue);
	  


         }

        //======================================================================
        // solution of linear system
        //======================================================================
        
        memset(defect, 0, N_Unknowns*SizeOfDouble);
        SQMATRICES[0] = MatricesM[mg_level-1];
           
        // compute defect
        Defect(sqmatrices,NULL,sol,B,defect);

        residual =  Ddot(N_Unknowns, defect, defect);
        residual = sqrt(residual);
        
        OutPut("initial residual ");
        OutPut(setw(14) << residual << endl);
        
        OutPut("nonlinear step " << j << " residual ");
        OutPut(setw(14) << residual << endl);
        if ((((residual<= limit) || (j >= Max_It-1)))
          && (j>=1))
        {
          if (j==Max_It-1)
            j++;
          OutPut("ITE : " << setw(3) << j);
          OutPut(" (" << setw(3) << N_LinIterCurr << "/");
          OutPut(setw(3) << N_LinIter << " LINITE)");
          OutPut("  TIME FOR SOLVER : " << solver_time_curr << "/" << solver_time << "s");
          OutPut("  RES : " <<  residual << endl);
          // count total running time
          t4 =  GetTime();
          total_time += t4 - t3;
          t3 = t4;
          OutPut(TDatabase::TimeDB->CURRENTTIME << " total elapsed running time " <<    total_time << endl);
	  // save oldrhs 
          if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE==FEM_FCT)
	  {
	      memcpy(oldrhs_fem_fct0, oldrhs_fem_fct1, N_Unknowns*SizeOfDouble);
	  }
	  if (!(TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_FCT))
	  {
	      if ((TDatabase::ParamDB->DISCTYPE == SDFEM)&&(TDatabase::ParamDB->SOLD_TYPE))
		  MatAdd(MatricesM[mg_level-1], MatricesS[mg_level-1],-1); 
	      
	      // restore mass matrices by subtracting the A-matrices
	      for(i=0;i<mg_level;i++)
	      {
		  MatAdd(MatricesM[i], MatricesA[i], -gamma);
	      }
	      // set current factor of steady state matrix
	      gamma = 0;
	      
	      if(TDatabase::ParamDB->DISCTYPE == SDFEM)
	      {
		  for(i=0;i<mg_level;i++)
		      MatAdd(MatricesM[i], MatricesK[i], -1);
	      }
	  }
         
          break;
        }
        j++;
	// save current solution for damping
	if (!(TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME))
	    memcpy(oldsol_nlinite, sol, N_Unknowns*SizeOfDouble);
        //======================================================================
        // solve linear system
        //======================================================================
        switch(TDatabase::ParamDB->SOLVER_TYPE)
        {
          case AMG:
	    TDatabase::ParamDB->SC_VERBOSE=0;     
            t1 = GetTime();
            Solver(SQMATRICES[0], B, sol);
            t2 = GetTime();
            solver_time_curr = t2-t1;
            solver_time += t2-t1;
          break;
             
          case GMG:
            t1 = GetTime();
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
            {
               memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
               memcpy(itmethod_rhs, B, N_Unknowns*SizeOfDouble);
            }
           
            N_LinIterCurrIte = itmethod->Iterate(sqmatrices,NULL,itmethod_sol,itmethod_rhs);
	    N_LinIterCurr += N_LinIterCurrIte;
	    N_LinIter += N_LinIterCurrIte ;
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
            {
               memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
               memcpy(B, itmethod_rhs, N_Unknowns*SizeOfDouble);
            }
            t2 = GetTime();
            solver_time_curr += t2-t1;
            solver_time += t2-t1;
          break; 
        } // endswitch SOLVER_TYPE
/*	OutPut("solution " << sqrt(Ddot(N_Unknowns,sol,sol)) << endl);
	int j1,j2,j3;
	double xx,xy,xz,valxx[4];
	for (j1=0;j1<33;j1++)
	{
	    xx = j1/32.0;
	    for (j2=0;j2<33;j2++)
	    {
		xy = j2/32.0;
		for (j3=0;j3<33;j3++)
		{
		    xz = j3/32.0;
		    SolArray[mg_level-1]->FindGradient(xx,xy,xz,valxx);
		    OutPut(xx << " "<< xy << " " << xz << " :ex  " << 
			   sin(TDatabase::TimeDB->CURRENTTIME)*(sin(2*Pi*xx)*sin(2*Pi*xy)*sin(2*Pi*xz)+1)
			   << " :sol " << valxx[0] << endl);
		}
	    }
	}*/

        //======================================================================
        // end solve linear system
        //======================================================================
	// reset matrices for SUPG and SOLD
        if (!(TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_FCT))
        {
	  if ((TDatabase::ParamDB->DISCTYPE == SDFEM)&&(TDatabase::ParamDB->SOLD_TYPE))
	      MatAdd(MatricesM[mg_level-1], MatricesS[mg_level-1],-1); 

          // restore mass matrices by subtracting the A-matrices
          for(i=0;i<mg_level;i++)
          {
            MatAdd(MatricesM[i], MatricesA[i], -gamma);
          }
          // set current factor of steady state matrix
          gamma = 0;

          if(TDatabase::ParamDB->DISCTYPE == SDFEM)
          {
            for(i=0;i<mg_level;i++)
              MatAdd(MatricesM[i], MatricesK[i], -1);
          }
	}// end if !FEM_FCT
       if ((TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME))
	{
	  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE==FEM_FCT)
	  {
	      memcpy(oldrhs_fem_fct0, oldrhs_fem_fct1, N_Unknowns*SizeOfDouble);
	  }
          break;
	}
      }   //end iteration for solving nonlinear problem
  
      // damping
      if (fabs(damp_nlinite - 1.0)>1e-6)
      {
	  for (ii=0;ii<N_Unknowns;ii++)
	      sol[ii] = damp_nlinite * sol[ii] + (1.0-damp_nlinite)*oldsol_nlinite[ii];
      }

    // save solution
     memcpy(oldsol,sol,N_Unknowns*SizeOfDouble);
    }   // endfor l (sub steps of fractional step theta)                                      
     
      // measure errors to known solution
      if(TDatabase::ParamDB->MEASURE_ERRORS)
      {
        aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
        conc->GetErrors(Exact, 4, AllDerivatives, 2, L2H1Errors, 
                        BilinearCoeffs, aux, 1, fesp, errors);
        delete aux;
  
        OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
        OutPut(" L2: " << errors[0]);
        OutPut(" H1-semi: " << errors[1] << endl);
      } // endif MEASURE_ERRORS
  
    } // endfor two time discs of adaptive time step control


#ifdef  __BAIL_3D__
    ComputeExtremalValues(N_Unknowns,sol,errors);
    OutPut(TDatabase::TimeDB->CURRENTTIME <<  " min: " << errors[0] << " max: " << errors[1] 
	   << " " << errors[1] - errors[0] <<endl );
#endif // 
  


#ifdef __FUEL_CELL__
    max = -100;
    min =  100;
    for(i=0;i<N_Unknowns;i++)
    {
      if(sol[i] < min) min = sol[i];
      if(sol[i] > max) max = sol[i];
    }

    OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
    x = 17; y = 120; z = 0;
    BoundValue(x, y, z, t1);
    OutPut(" Inflow: " << t1);
    OutPut(" Minimum: " << min << " Maximum: " << max << endl);

/*
    for(i=0;i<N_Unknowns;i++)
    {
      if(sol[i] < 0) sol[i] = 0;
      if(sol[i] > t1) sol[i] = t1;
    }
*/

    GetOutFlow(N_OutFlowCells, OutFlowCells, OutFlowFaces, conc,
               outflow, area);
    OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
    //OutPut("outflow: " << outflow << " area: " << area);
    OutPut(" mean value: " << outflow/area << endl);
#endif // __FUEL_CELL__

    if ((TDatabase::ParamDB->WRITE_GRAPE)
	||(TDatabase::ParamDB->WRITE_GMV)
	||(TDatabase::ParamDB->WRITE_VTK))
    {
      picture = 0;
      if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      {
         if (TDatabase::ParamDB->WRITE_GRAPE)
         {
            os.seekp(std::ios::beg);
            os << GrapeBaseName << m << ".dat" << ends;
            Output->WriteGrape(os.str().c_str());
            picture++;
         }
         if (TDatabase::ParamDB->WRITE_GMV)
         {
            os.seekp(std::ios::beg);
            os << GmvBaseName << m << ".gmv" << ends;
            Output->WriteGMV(os.str().c_str());
            picture++;
         }

	if(TDatabase::ParamDB->WRITE_VTK)
	{
	    os.seekp(std::ios::beg);
	    os << VtkBaseName << m << ".vtk" << ends;
	    Output->WriteVtk(os.str().c_str());
	}
      }
    }
  } // while
//======================================================================
// end of time cycle
//======================================================================

  if(TDatabase::ParamDB->WRITE_GRAPE)
  {
    os.seekp(std::ios::beg);
    os << GrapeBaseName <<  "end." << m << ".dat" << ends;
    Output->WriteGrape(os.str().c_str());
  }
  
  if(TDatabase::ParamDB->WRITE_GMV)
  {
    os.seekp(std::ios::beg);
    os << GmvBaseName <<  "end." << m << ".dat" << ends;
    Output->WriteGMV(os.str().c_str());
  }
 
  if(TDatabase::ParamDB->WRITE_VTK)
  {
      os.seekp(std::ios::beg);
      os << VtkBaseName << "end." << m << ".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
  }

  t4 =  GetTime();
  total_time += t4 - t3;
  OutPut("total running time: " << total_time << endl);
  CloseFiles();
  return 0;
}

