// =======================================================================
//
// Purpose:     main program
//
// Author:      Volker Behns  22.07.97
//
// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <AuxParam2D.h>
#include <QuadAffin.h>
#include <Solver.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <DiscreteForm2D.h>
#include <LinAlg.h>
#include <Collection.h>
#include <TCD2D.h>
//#include <ConvDiff2D_Routines.h>
#include <LocalProjection.h>
#include <TimeConvDiff2D.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <malloc.h>

double bound = 0;

#include <Upwind.h>
#include <FEM_TVD_FCT.h>
#include <ItMethod.h>
#include <FixedPointIte.h>
#include <MultiGridScaIte.h>
#include <FgmresIte.h>

#include <MultiGrid2D.h>
#include <MGLevel2D.h>

#include <MainUtilities.h>
#include <TimeUtilities.h>

#define AMG 0
#define GMG 1
#define DIRECT 2

// =======================================================================
// include current example
// =======================================================================

//#include "../Examples/TCD_2D/ansatz1.h"
//#include "../Examples/TCD_2D/SinCos1.h"
//#include "../Examples/TCD_2D/SinCos2.h"
//#include "../Examples/TCD_2D/SinCos4.h"
//#include "../Examples/TCD_2D/Robin_indian.01.h"
//#include "../Examples/TCD_2D/Robin_indian.00.h"
//#include "../Examples/TCD_2D/Robin_indian.02.h"
//#include "../Examples/TCD_2D/Robin_indian.03.h"
//#include "../Examples/TCD_2D/Robin_indian.04.h"
//#include "../Examples/TCD_2D/Robin_scaled.00.h"
#include "../Examples/TCD_2D/RotatingBodies.h"
//#include "../Examples/TCD_2D/Bulk_Academic_Test.h"
//#include "../Examples/TCD_2D/JohnMaubachTobiska1997inst.h"
//#include "../Examples/TCD_2D/JohnMaubachTobiska1997inst2.h"
//#include "../Examples/ConvDiff2D/Hemker1996.h"
//#include "../Examples/TCD_2D/Sin3.h"
//#include "../Examples/TCD_2D/Sin3_0_pi.h"
//#include "../Examples/TCD_2D/exp_0.h"
//#include "../Examples/TCD_2D/exp_sin_pi.h"
//#include "../Examples/TCD_2D/Sin3_0_pi_stabil.h"
//#include "../Examples/TCD_2D/Production.h"
//#include "../Examples/TCD_2D/Impulse.h"

// ======================================================================
// utilities for main program
// ======================================================================

int main(int argc, char* argv[])
{
  // ======================================================================
  // variable declaration
  // ======================================================================

  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *coll, *mortarcoll = NULL;
  TBaseCell *cell;
  TFESpace2D *concentration_space, *pressure_space, *velocity_space;
  TFESpace2D **ConcentrationSpaces, **VelocitySpaces;
  TOutput2D *Output;

  double *B, *rhs, *sol, *oldsol, tol, tolerance, *defect, *startsol, *frac_step_sol;
  double *oldrhs, *itmethod_sol, *itmethod_rhs, *current_sol, *current_B;
  double *sol_velo, *lump_mass, *matrix_D_Entries, *oldrhs_fem_fct0, *tilde_u;
  double  *oldrhs_fem_fct1, *rhs_edge, *oldsol_nlinite, damp_nlinite;
  double *entries_mass_matrix, tau_small = -1;
  int i,j,k,l,m,n, N_, Len, low, only_first_time;
  int N_Rows, N_Columns, N_Unknowns, N_Unknowns_Velo, N_P, sold_parameter_type;
  double *l2u1, *l2u2, *h1u1, *h1u2;
  double *l2p, *h1p, *sd, *l_inf, vect[3], exactvect[3];
  int which, *permutation;
  double DiffL2, DiffH1, t;
  char *PRM, *GEO;
  int LEVELS, ORDER, order;
  int ret, pde;
  double negPower;
  double x,y,max,min,sum;
  double tau1, tau2, tau_array[4];
  double errors[9], olderror, olderror1, olderror2, p1, p2;
  double errors_smooth[9], olderror_smooth[3];
  double t1, t2, res, res2, oldres, solver_time, solver_time_curr, residual, oldresidual;
  double total_time, t3, t4;
  double impuls_residual,linredfac;
  int N_LinIter, N_LinIterCurr, N_LinIterCurrIte, N_SubSteps, N_Active, n_aux;
  double gamma, tau, oldtau;
  int *RowPtr;

  std::ostringstream os;
  char *PsBaseName, *GnuBaseName;
  char *VtkBaseName, *MatlabBaseName, *GmvBaseName;

  TFEFunction2D *conc, **velo1, **velo2, *fefct[3];
  TFEFunction2D **SolArray, **AuxFEFunctArray, *pressure;
  TFEFunction2D *oldconc, **OldSolArray;
  TFEVectFunct2D **velocity, **AuxFEVectFunctArray;
  TFESpace2D *fesp[2], *ferhs[1];
  double delta, end_time;

  TAuxParam2D *aux;

  TSquareStructure2D *sqstructureA;
  TSquareMatrix2D *sqmatrixA, *SQMATRICES[3];
  TSquareMatrix2D *sqmatrixM , *sqmatrixK, *sqmatrixS;
  TSquareMatrix2D **MatricesA, **MatricesM, **MatricesK, **MatricesS;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  MatVecProc *MatVect;
  DefectProc *Defect;

  double **RhsArray;

  TMGLevel2D *MGLevel;
  TMultiGrid2D *MG;

  double *RHSs[3];
  int *N_Uarray;

  // discrete forms have to be created
  TDiscreteForm2D *DiscreteForm;

  TDiscreteForm2D *DiscreteFormMatrixMRhs;
  TDiscreteForm2D *DiscreteFormMatrixMRhs_SUPG;

  TDiscreteForm2D *DiscreteFormMatrixARhs;
  TDiscreteForm2D *DiscreteFormMatricesAKRhs_SUPG;
  TDiscreteForm2D *DiscreteFormMatricesAKRhs_SOLD;

  int N_SquareMatrices, N_Rhs, N_FESpaces;

  BoundCondFunct2D *BoundaryConditions[1];
  BoundValueFunct2D *BoundValues[1];
  CoeffFct2D *Coefficients[1];

  TItMethod *itmethod, *prec;
  int FirstSolve, Max_It;
  int N_Paramters=1, time_discs;
  double Parameters[2], hmin, hmax, limit;

  int N_UConv, level_down, ii, fully_implicit = 0;
  int mg_level,mg_type,CurrentDiscType, last_sq, step_length;
  int very_first_time=0, zerostart, comp_vort, picture;
  int *neum_to_diri, N_neum_to_diri = 0, *neum_to_diri_bdry;
  double *neum_to_diri_param;

  char UString[] = "u";
  char PString[] = "p";
  char Readin[] = "readin.dat";
  char MassMatrix[] = "Mass matrix";
  char Mass[] = "Mass";
  char Name[] = "name";

  // ======================================================================
  // end of variable declaration
  // ======================================================================

  os << " ";

  //======================================================================
  // read parameter file
  //======================================================================
  total_time = GetTime();
  if(argc>=2)
    ret=Domain->ReadParam(argv[1]);
  else
    ret=Domain->ReadParam(Readin);

  if(ret==-1)
  {
    exit(-1);
  }

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  if( (TDatabase::ParamDB->DISCTYPE>=5)&&(TDatabase::ParamDB->DISCTYPE!=LOCAL_PROJECTION) )
  {
    OutPut("DISCTYPE " << TDatabase::ParamDB->DISCTYPE << " NOT IMPLEMENTED!" << endl);
    exit(4711);
  }
  // set some parameters
  if(TDatabase::ParamDB->DISCTYPE != LOCAL_PROJECTION)
  {
    // switch off all local projection terms
    TDatabase::ParamDB->LP_FULL_GRADIENT = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT = 1;

    TDatabase::ParamDB->LP_STREAMLINE = 0;
    TDatabase::ParamDB->LP_STREAMLINE_COEFF = 0;
    TDatabase::ParamDB->LP_STREAMLINE_EXPONENT = 1;
  }
  if(TDatabase::ParamDB->DISCTYPE != LOCAL_PROJECTION)
  {
    // switch off all local projection terms
    TDatabase::ParamDB->LP_FULL_GRADIENT = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT = 1;

    TDatabase::ParamDB->LP_STREAMLINE = 0;
    TDatabase::ParamDB->LP_STREAMLINE_COEFF = 0;
    TDatabase::ParamDB->LP_STREAMLINE_EXPONENT = 1;
  }

  if(TDatabase::ParamDB->LP_FULL_GRADIENT)
  {
    if(TDatabase::ParamDB->LP_STREAMLINE)
    {
      TDatabase::ParamDB->LP_STREAMLINE = 0;
      TDatabase::ParamDB->LP_STREAMLINE_COEFF = 0;
      TDatabase::ParamDB->LP_STREAMLINE_EXPONENT = 1;
      OutPut("local projection stabilisation in streamline direction ");
      OutPut("is switched off due to stabilisation of full gradient." << endl);
    }
  }

  if(TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE == -123)
    TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;

  if(TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE == -123)
    TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;

  // write parameters into outfile
  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();

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
      TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME = TDatabase::ParamDB->FEM_FCT_LINEAR_TYPE;
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE=FEM_FCT;
      sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;
  }

  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;

  // assign names for output files
  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  GnuBaseName = TDatabase::ParamDB->GNUBASENAME;
  GmvBaseName = TDatabase::ParamDB->GMVBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  MatlabBaseName = TDatabase::ParamDB->MATLABBASENAME;

  // set type of multilevel
  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SCALAR;
  if (mg_type)
    mg_level = 1;
  else
    mg_level = 0;

  // number of levels
  LEVELS = TDatabase::ParamDB->LEVELS;

  // allocate arrays for error computation
  l2u1 = new double[LEVELS+1];
  l2u2 = new double[LEVELS+1];
  h1u1 = new double[LEVELS+1];
  h1u2 = new double[LEVELS+1];
  sd = new double[LEVELS+1];
  l_inf = new double[LEVELS+1];

  // array for pointers to the solutio on the
  // different levels of the multigrid
  SolArray = new TFEFunction2D*[LEVELS+1];
  OldSolArray = new TFEFunction2D*[LEVELS+1];

  // array for pointers to right hand sides on the
  // different levels of the multigrid
  RhsArray = new double* [LEVELS+1];
  N_Uarray = new int[LEVELS+1];

  // allocate pointers to arrays for velocity which
  // may play the role of a convection
  velocity = new TFEVectFunct2D*[LEVELS+1];
  velo1 = new TFEFunction2D*[LEVELS+1];
  velo2 = new TFEFunction2D*[LEVELS+1];

  // array which points to the finite element spaces on the
  // different levels of the multigrid
  ConcentrationSpaces = new TFESpace2D*[LEVELS+1];

  // array which points to the finite element spaces of the
  // different levels of the multigrid
  VelocitySpaces = new TFESpace2D*[LEVELS+1];

  // array which points to the system matrices on  the
  // different levels of the multigrid
  MatricesA = new TSquareMatrix2D*[LEVELS+1];

  // array which points to the mass matrices on  the
  // different levels of the multigrid
  MatricesM = new TSquareMatrix2D*[LEVELS+1];

  // array which points to the stabilization  matrices (sdfem) on the
  // different levels of the multigrid
  MatricesK = new TSquareMatrix2D*[LEVELS+1];

  // array which points to the stabilization  matrices (sold) on the
  // different levels of the multigrid
  // it is actually used only on the finest level
  MatricesS = new TSquareMatrix2D*[LEVELS+1];

  // pointers to the routines which compute matrix-vector
  // products and the defect
  MatVect = MatVect_Scalar;
  Defect = Defect_Scalar;

  //======================================================================
  // initialize discrete forms
  //======================================================================

  // discrete form for assembling mass matrix and rhs (Galerkin FEM)
  DiscreteFormMatrixMRhs = new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_MatrixMRhs, Derivatives_MatrixMRhs,
    SpacesNumbers_MatrixMRhs, N_Matrices_MatrixMRhs, N_Rhs_MatrixMRhs,
    RowSpace_MatrixMRhs, ColumnSpace_MatrixMRhs, RhsSpace_MatrixMRhs,
    MatrixMRhsAssemble, BilinearCoeffs, NULL);

  // discrete form for assembling mass matrix and rhs (SDFEM)
  DiscreteFormMatrixMRhs_SUPG = new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_MatrixMRhs_SUPG, Derivatives_MatrixMRhs_SUPG,
    SpacesNumbers_MatrixMRhs_SUPG, N_Matrices_MatrixMRhs_SUPG, N_Rhs_MatrixMRhs_SUPG,
    RowSpace_MatrixMRhs_SUPG, ColumnSpace_MatrixMRhs_SUPG, RhsSpace_MatrixMRhs_SUPG,
    MatrixMRhsAssemble_SUPG, BilinearCoeffs, NULL);

  // discrete form for assembling stiffness matrix and rhs (Galerkin FEM)
  DiscreteFormMatrixARhs = new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_MatrixARhs, Derivatives_MatrixARhs,
    SpacesNumbers_MatrixARhs, N_Matrices_MatrixARhs, N_Rhs_MatrixARhs,
    RowSpace_MatrixARhs, ColumnSpace_MatrixARhs, RhsSpace_MatrixARhs,
    MatrixARhsAssemble, BilinearCoeffs, NULL);

  // discrete form for assembling stiffness matrix, stabilization matrix and rhs (SDFEM)
  DiscreteFormMatricesAKRhs_SUPG = new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_MatricesAKRhs_SUPG, Derivatives_MatricesAKRhs_SUPG,
    SpacesNumbers_MatricesAKRhs_SUPG, N_Matrices_MatricesAKRhs_SUPG, N_Rhs_MatricesAKRhs_SUPG,
    RowSpace_MatricesAKRhs_SUPG, ColumnSpace_MatricesAKRhs_SUPG, RhsSpace_MatricesAKRhs_SUPG,
    MatricesAKRhsAssemble_SUPG, BilinearCoeffs, NULL);

  // discrete form for assembling stiffness matrix, stabilization matrix and rhs (SDFEM)
  DiscreteFormMatricesAKRhs_SOLD = new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_MatricesAKRhs_SUPG, Derivatives_MatricesAKRhs_SUPG,
    SpacesNumbers_MatricesAKRhs_SUPG, N_Matrices_MatricesAKRhs_SOLD, N_Rhs_MatricesAKRhs_SUPG,
    RowSpace_MatricesAKRhs_SOLD, ColumnSpace_MatricesAKRhs_SOLD, RhsSpace_MatricesAKRhs_SUPG,
    MatricesAKRhsAssemble_SUPG, BilinearCoeffs, NULL);

  //======================================================================
  // read boundary parameterization and initialize coarse grid
  //======================================================================
  Domain->Init(PRM, GEO);
  Domain->ConvertQuadToTri(TDatabase::ParamDB->CONVERT_QUAD_TO_TRI);

  // refine grid up to the coarsest level
  for(i=0;i<TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR;i++)
    Domain->RegRefineAll();

  // initialize multigrid
  if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
    Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
    i=1;
    MG = new TMultiGrid2D(i, N_Paramters, Parameters);
  }
  mg_level = LEVELS+mg_level;
  damp_nlinite = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SCALAR;

  // initializ time
  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = TDatabase::TimeDB->TIMESTEPLENGTH;
  SetTimeDiscParameters();
  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR;

  t3 = GetTime();
  total_time = t3 - total_time;
  SetPolynomialDegree();
  //======================================================================
  // loop over all levels
  // do all refinements
  // build all fe spaces
  // load velocity on the finest level (if required)
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
      concentration_space = new TFESpace2D(coll,Name,UString,BoundCondition,-1,NULL);
    }
    // standard multigrid or finest level
    // get fe space of high order disc on finest geo grid
    else
    {
      ORDER  = TDatabase::ParamDB->ANSATZ_ORDER;
      concentration_space = new TFESpace2D(coll, Name,
        UString, BoundCondition, ORDER, NULL);
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
    sqstructureA = new TSquareStructure2D(concentration_space);
    sqstructureA->Sort();

    // two matrices used
    // A contains the non time dependent part of the discretization
    sqmatrixA = new TSquareMatrix2D(sqstructureA);
    MatricesA[i] = sqmatrixA;

    // M is the mass matrix
    // the iterative solver uses M
    sqmatrixM = new TSquareMatrix2D(sqstructureA);
    MatricesM[i] = sqmatrixM;

    if(TDatabase::ParamDB->DISCTYPE == SDFEM)
    {
      // stabilisation matrix K
      sqmatrixK = new TSquareMatrix2D(sqstructureA);
      MatricesK[i] = sqmatrixK;
      if(TDatabase::ParamDB->SOLD_TYPE)
      {
	  // stabilisation matrix S
	  sqmatrixS = new TSquareMatrix2D(sqstructureA);
	  MatricesS[i] = sqmatrixS;
      }
      if ((sold_parameter_type >= BH04)&& (sold_parameter_type <= BE05_2))
          rhs_edge = new double[N_Unknowns];
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
    conc = new TFEFunction2D(concentration_space, UString, UString, sol, N_Unknowns);
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
      oldconc =  new TFEFunction2D(concentration_space, UString, UString, oldsol, N_Unknowns);
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
      sqmatrixK = new TSquareMatrix2D(sqstructureA);
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
      MGLevel = new TMGLevel2D(i, sqmatrixM, current_B, current_sol,  n_aux, NULL);
      MG->AddLevel(MGLevel);
    }

    // define everything for velocity (convection)
    // array which contains the values
    if (TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
    {
      ORDER  = TDatabase::ParamDB->ANSATZ_ORDER;
      velocity_space = new TFESpace2D(coll, Name,
        UString, BoundConditionNSE, ORDER, NULL);
      VelocitySpaces[i] = velocity_space;

      N_Unknowns_Velo = velocity_space->GetN_DegreesOfFreedom();
      sol_velo = new double[2*N_Unknowns_Velo];
      memset(sol_velo, 0, 2*N_Unknowns_Velo*SizeOfDouble);

      // vector fe function
      velocity[i] = new TFEVectFunct2D(velocity_space, UString, UString, sol_velo, N_Unknowns, 2);
      // individual components of velocity
      velo1[i] = velocity[i]->GetComponent(0);
      velo2[i] = velocity[i]->GetComponent(1);
    }

    // prepare output, only the concentration will be saved
    if (i==mg_level-1)
    {
      Output = new TOutput2D(1, 1, 0, 1, Domain);
      Output->AddFEFunction(conc);
      os.seekp(std::ios::beg);
      Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());
    }

    // interpolate initial condition
    conc->Interpolate(InitialCondition);
  }                                               // endfor i

  //======================================================================
  // loop over all levels
  // restrict velocity
  //======================================================================

  if (TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
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
  }                                               // endfor i

  //======================================================================
  // all data are available for assembling matrices
  // loop over all levels
  //======================================================================
  for(i=0;i<mg_level;i++)
  {
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

      aux =  new TAuxParam2D(TimeCDParamsVeloFieldN_FESpaces,
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
    else
    {
      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
    }
    //======================================================================
    // assembling of mass matrix and rhs
    //======================================================================
    // reset matrices
    N_SquareMatrices = 1;
    SQMATRICES[0] = MatricesM[i];
    SQMATRICES[0]->Reset();

    // SUPG term of rhs is assembled
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

    Assemble2D(N_FESpaces, fesp,
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
#ifdef __ROTATING_BODIES__
      CheckWrongNeumannNodes(coll, ConcentrationSpaces[i], N_neum_to_diri, neum_to_diri,
        neum_to_diri_bdry, neum_to_diri_param);
#endif
#ifdef __SINCOS1__
      CheckWrongNeumannNodes(coll, ConcentrationSpaces[i], N_neum_to_diri, neum_to_diri,
        neum_to_diri_bdry, neum_to_diri_param);
#endif
#ifdef __JOHNMAUBACHTOBISKA__
      CheckWrongNeumannNodes(coll, ConcentrationSpaces[i], N_neum_to_diri, neum_to_diri,
        neum_to_diri_bdry, neum_to_diri_param);
#endif
#ifdef __HEMKER1996__
      CheckWrongNeumannNodes(Domain, i, ConcentrationSpaces[i], N_neum_to_diri, neum_to_diri,
        neum_to_diri_bdry, neum_to_diri_param);
#endif

#ifdef __BULK_ACAD_TEST__
      CheckWrongNeumannNodes(coll, ConcentrationSpaces[i], N_neum_to_diri, neum_to_diri,
        neum_to_diri_bdry, neum_to_diri_param);
#endif
#ifdef __SIN3__
      CheckWrongNeumannNodes(coll, ConcentrationSpaces[i], N_neum_to_diri, neum_to_diri,
        neum_to_diri_bdry, neum_to_diri_param);
#endif
#ifdef __TCD_2D_PRODUCTION__
      CheckWrongNeumannNodes(coll, ConcentrationSpaces[i], N_neum_to_diri, neum_to_diri,
        neum_to_diri_bdry, neum_to_diri_param);
#endif
#ifdef __IMPULSE__
      CheckWrongNeumannNodes(coll, ConcentrationSpaces[i], N_neum_to_diri, neum_to_diri,
        neum_to_diri_bdry, neum_to_diri_param);
#endif

    }
#ifdef __ROTATING_BODIES__
    else
    {
	// needed for computing l1 and l2 error
      lump_mass = new double [N_Unknowns];
      LumpMassMatrixToVector(MatricesM[i], lump_mass);
    }
#endif
  }                                               // endfor i
  //======================================================================
  // end of space cycle, finest grid reached
  // everything happens now on the same grid
  //======================================================================
  // save solution
  memcpy(oldsol,sol,N_Unknowns*SizeOfDouble);

  // initialize
  for (i=0;i<9;i++)
      errors[i] = errors_smooth[i] = 0;
  olderror=0; 
  olderror1=0;
  olderror2=0;
  olderror_smooth[0] = olderror_smooth[1] = olderror_smooth[2] = 0;
  picture = 0;

  if (TDatabase::ParamDB->WRITE_GNU)
  {
    os.seekp(std::ios::beg);
    os << GnuBaseName << 0 << ".gnu" << ends;
    Output->WriteGnuplot(os.str().c_str());
    picture = 1;
  }

  if(TDatabase::ParamDB->WRITE_GMV)
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
  
  if(TDatabase::ParamDB->WRITE_MATLAB)
  {
      os.seekp(std::ios::beg);
      os << MatlabBaseName << 0 << ".m" << ends;
      Output->WriteMatlab(os.str().c_str());
  }
  
  // allocate arrays for solver
  defect = new double[N_Unknowns];
  startsol = new double[N_Unknowns];
  frac_step_sol = new double[N_Unknowns];

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
          0, N_Unknowns, 0);
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

  // set some scalings for paper with Julia Novo
  if (TDatabase::ParamDB->P9 == 191)
  {
      hmin/=sqrt(2);
      TDatabase::TimeDB->TIMESTEPLENGTH = pow(hmin,2.0*(ORDER+1)/3.0);
      //TDatabase::TimeDB->TIMESTEPLENGTH = hmin*hmin;
      //TDatabase::TimeDB->TIMESTEPLENGTH = pow(hmin,0.8*(ORDER+0.5));
       TDatabase::TimeDB->TIMESTEPLENGTH = pow(hmin,1.5);
     OutPut("CHANGED TDatabase::TimeDB->TIMESTEPLENGTH TO " << TDatabase::TimeDB->TIMESTEPLENGTH 
	     << endl);
  }  
  /* TDatabase::TimeDB->CURRENTTIME = 1e-4;
  Exact(0.55,0.4,errors);
  OutPut(errors[0] << " "<< errors[1] << " "<< errors[2] <<  " " << errors[3] << endl);
  double X[1], Y[1], *coeffss[1];

  X[0] = 0.55;
  Y[0] = 0.4;
  coeffss[0] = new double[7];
  BilinearCoeffs(1, X, Y, NULL, coeffss); 
  exit(4711);
  */
  if (TDatabase::TimeDB->CURRENTTIMESTEPLENGTH < tau_small)
  {
      entries_mass_matrix = new double[MatricesM[mg_level-1]->GetN_Entries()];
      memcpy(MatricesM[mg_level-1]->GetEntries(), entries_mass_matrix,
	     MatricesM[mg_level-1]->GetN_Entries() * SizeOfDouble);
  }
  
  //======================================================================
  // start of time cycle
  // everything happens on the same grid
  //======================================================================
  while(TDatabase::TimeDB->CURRENTTIME< end_time)
  {                                               // time cycle
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

    for(l=0;l<N_SubSteps;l++)                     // sub steps of fractional step theta
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
      OutPut("memory: " << setw(10) << GetMemory() << endl);

      // edge stabilization SOLD methods
      if ((sold_parameter_type >= BH04)&&(sold_parameter_type <= BE05_2))
	  memset(rhs_edge, 0,  N_Unknowns*SizeOfDouble);

      // start iteration for solving nonlinear problems
      while(1)
      {
	  // working array for rhs is B, initialize B
	  memset(B, 0, N_Unknowns*SizeOfDouble);
	  // compute terms with data from previous time step
	  // old rhs multiplied with current subtime step and theta3 on B
	  Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, rhs, B);
  
        // assembling of A and rhs
        // for SDFEM: in addition stabilisation matrix K
        if (TDatabase::ParamDB->INTERNAL_STEADY_STATE_MATRICES_OR_RHS<=1)
        {
          for(i=0;i<mg_level;i++)
          {
	      if ((TDatabase::ParamDB->SOLD_PARAMETER_TYPE==FEM_FCT) && (j>0))
		  break;
            // set parameters
            N_Rhs = 1;
            N_FESpaces = 1;
            fesp[0] = ConcentrationSpaces[i];
            // convection is given as finite element velocity field
            if (TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
            {
              N_FESpaces = 2;
              fesp[1] = VelocitySpaces[i];
              fefct[0] = velo1[i];
              fefct[1] = velo2[i];

              aux =  new TAuxParam2D(TimeCDParamsVeloFieldN_FESpaces,
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

              aux =  new TAuxParam2D(TimeCDParamsSOLDN_FESpaces,
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
              aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
            }

            //======================================================================
            // assembling of convection-diffusion matrix and rhs
            //======================================================================
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
              DiscreteForm = DiscreteFormMatrixARhs;
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

            Assemble2D(N_FESpaces, fesp,
              N_SquareMatrices, SQMATRICES,
              0, NULL,
              N_Rhs, RHSs, ferhs,
              DiscreteForm,
              BoundaryConditions,
              BoundValues,
              aux);
            delete aux;
	   
	    // save rhs without Dirichlet values
	    if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE==FEM_FCT)
	    {
		BoundValues[0] = BoundValue;
		memcpy(oldrhs_fem_fct1, rhs, N_Unknowns*SizeOfDouble);
	    }

            if ((TDatabase::ParamDB->SOLD_TYPE)&&(i==mg_level-1))
            {
		// set parameter back to zero
              TDatabase::ParamDB->INTERNAL_SOLD_ACTIVE = 0;
	    } // end SOLD type
	    // LPS 
	    if(TDatabase::ParamDB->DISCTYPE == LOCAL_PROJECTION)
	    {
		//if(TDatabase::ParamDB->LP_FULL_GRADIENT)
		//  UltraLocalProjection(MatricesA[i], FALSE, Coefficients[0]);
		
		if(TDatabase::ParamDB->LP_STREAMLINE)
		{
		    OutPut("local projection stabilisation in streamline direction ");
		    OutPut("is currently not available." << endl);
		    exit(4711);
		}
	    }
	  }                                       // endfor i
          if (TDatabase::ParamDB->INTERNAL_STEADY_STATE_MATRICES_OR_RHS==1)
            TDatabase::ParamDB->INTERNAL_STEADY_STATE_MATRICES_OR_RHS = 4711;
        }
        if (very_first_time==1)
        {
          very_first_time=0;
          l--;
          continue;
        }
        // add rhs from current sub time step to rhs array B
        Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, rhs, B);

	if (TDatabase::ParamDB->SOLD_TYPE)
	{
	    // prepare rhs for edge stabilizations, only on the finest grid
	    if ((sold_parameter_type >= BH04)&&(sold_parameter_type <=  BE05_2))
	    {
		// restore rhs
		//Daxpy(N_Active, 1, rhs_edge, rhs);
		memset(rhs_edge, 0,  N_Unknowns*SizeOfDouble);
		tau_array[0] = tau*TDatabase::TimeDB->THETA1;
		tau_array[1] = tau*TDatabase::TimeDB->THETA2;
		tau_array[2] = tau*TDatabase::TimeDB->THETA3;
		tau_array[3] = tau*TDatabase::TimeDB->THETA4;
		//EdgeStabilization(ConcentrationSpaces[mg_level-1], SolArray[mg_level-1],
		//	  Coefficients[0], rhs_edge, 1, 
		//	  tau_array, oldconc);
		// subtract rhs_edge from rhs
		OutPut("EdgeStabilization removed !!!" << endl);
		exit(4711);
		Daxpy(N_Active, -1, rhs_edge, B);
	      }
	} // end SOLD type
      	
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

          FEM_FCT_ForConvDiff(MatricesK[mg_level-1], MatricesA[mg_level-1],
			      N_Unknowns, N_Active, 
			      lump_mass, matrix_D_Entries, 
			      sol, oldsol, 
			      B, RhsArray[mg_level-1], oldrhs_fem_fct0, tilde_u,
			      N_neum_to_diri, neum_to_diri, 
			      neum_to_diri_bdry, neum_to_diri_param,
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
	  // set Diriclet nodes
	  if (N_neum_to_diri)
	      SetDirichletNodesFromNeumannNodes(SQMATRICES, B, sol, 
						N_neum_to_diri, neum_to_diri,
						neum_to_diri_bdry, neum_to_diri_param,
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

        OutPut("nonlinear step " << j << " residual ");
        OutPut(setw(14) << residual << endl);

        if ((((residual<= limit) || (j >= Max_It-1)))
          && (j>=1))
        {
          if (j==Max_It-1)
            j++;
          OutPut("ITE : " << setw(3) << j);
          //OutPut(" (" << setw(3) << N_LinIterCurr << "/");
          //OutPut(setw(3) << N_LinIter << " LINITE)");
          //OutPut("  TIME FOR SOLVER : " << solver_time_curr << "/" << solver_time << "s");
          OutPut("  RES : " <<  residual << endl);
          // count total running time
          //t4 =  GetTime();
          //total_time += t4 - t3;
          //t3 = t4;
          //OutPut(TDatabase::TimeDB->CURRENTTIME << " total elapsed running time " <<    total_time << endl);
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
          case DIRECT:
            t1 = GetTime();
            DirectSolver(SQMATRICES[0], B, sol);
            t2 = GetTime();
            OutPut( "time for direct solving: " << t2-t1 << endl);
            OutPut("solution " << sqrt(Ddot(N_Unknowns,sol,sol)) << endl);
            break;

          case AMG:
            TDatabase::ParamDB->SC_VERBOSE=0;
            t1 = GetTime();
            Solver(SQMATRICES[0], B, sol);
            t2 = GetTime();
            solver_time_curr = t2-t1;
            solver_time += solver_time_curr;
            break;

          case GMG:
            t1 = GetTime();
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
            {
              memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
              memcpy(itmethod_rhs, B, N_Unknowns*SizeOfDouble);
            }
            // for (j=0;j<N_Unknowns;j++)
            //  OutPut(j << " " << itmethod_rhs[j] << " " << rhs[j] << " " << itmethod_sol[j] << endl);

            N_LinIterCurrIte = itmethod->Iterate(sqmatrices,NULL,itmethod_sol,itmethod_rhs);
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
            {
              memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
              memcpy(B, itmethod_rhs, N_Unknowns*SizeOfDouble);
            }
            t2 = GetTime();
            solver_time_curr += t2-t1;
            solver_time += solver_time_curr;
            break;
        }                                         // endswitch SOLVER_TYPE

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
	  // for very small time steps, set original entries
	  if (TDatabase::TimeDB->CURRENTTIMESTEPLENGTH < tau_small)
	  {
	      entries_mass_matrix = new double[MatricesM[mg_level-1]->GetN_Entries()];
	      memcpy(entries_mass_matrix, MatricesM[mg_level-1]->GetEntries(),
		     MatricesM[mg_level-1]->GetN_Entries() * SizeOfDouble);
	  }
        }
        // only one iteration of linear schmes
        if ((TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME))
	{
	  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE==FEM_FCT)
	  {
	      memcpy(oldrhs_fem_fct0, oldrhs_fem_fct1, N_Unknowns*SizeOfDouble);
	  }
          break;
	}
	// damping
	if (fabs(damp_nlinite - 1.0)>1e-6)
	{
	    for (ii=0;ii<N_Unknowns;ii++)
		sol[ii] = damp_nlinite * sol[ii] + (1.0-damp_nlinite)*oldsol_nlinite[ii];
	}
      }                                           //end iteration for solving nonlinear problem

      // save solution
      memcpy(oldsol,sol,N_Unknowns*SizeOfDouble);
      /*	if (TDatabase::ParamDB->SOLD_TYPE)
      {
      for(i=0;i<mg_level-1;i++)
      memcpy(OldSolArray[i],SolArray[i],N_Uarray*SizeOfDouble);
      }*/
    }                                             // endfor l (sub steps of fractional step theta)

    /*****************************************************************************/
    /* postprocessing of current time step                                       */
    /*****************************************************************************/

    // measure errors to known solution
    if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      conc->GetErrors(Exact, 3, AllDerivatives, 3, SDFEMErrors,
        BilinearCoeffs, aux, 1, fesp, errors);
      delete aux;

      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
      OutPut(" L2: " << errors[0]);
      OutPut(" H1-semi: " << errors[1]);
      OutPut(" SD: " << errors[2] << endl);

      errors[4] += (errors[0]*errors[0] +olderror * olderror)*
	  TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      olderror = errors[0];
      OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(errors[4]) << " ");

      errors[5] += (errors[1]*errors[1] +olderror1 * olderror1)*
	  TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      OutPut( "L2(0,T;H1) " << sqrt(errors[5]));
      olderror1 = errors[1];

      errors[6] += tau*(errors[2]*errors[2] +olderror2 * olderror2);
	  //*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      OutPut( " SD(0,T) " << sqrt(errors[6]) << endl);
      olderror2 = errors[2];

#ifdef __ROTATING_BODIES__
  if(TDatabase::ParamDB->DISCTYPE != LOCAL_PROJECTION)
  {
      ComputeDiscreteErrors(ConcentrationSpaces[mg_level-1], SolArray[mg_level-1],
			    sol, errors, lump_mass);
      OutPut(TDatabase::TimeDB->CURRENTTIME <<  " l1: " << errors[0] << " l2: " << errors[1] <<endl );
      ComputeExtremalValues(N_Unknowns,sol,errors);
      OutPut(TDatabase::TimeDB->CURRENTTIME <<  " min: " << errors[0] << " max: " << errors[1] 
	     << " " << errors[1] - errors[0] <<endl );
      //exit(1);
  }
  else
  {
      //UltraLocalError_TCD(SolArray[mg_level-1], errors+5, lump_mass);
      OutPut(TDatabase::TimeDB->CURRENTTIME <<  " min: " << errors[7] << " max: " << errors[8] 
	     << " " << errors[8] - errors[7] <<endl );
  }
#endif
#ifdef __JOHNMAUBACHTOBISKA__
  if(TDatabase::ParamDB->DISCTYPE != LOCAL_PROJECTION)
  {
      ComputeExtremalValues(N_Unknowns,sol,errors);
      OutPut(TDatabase::TimeDB->CURRENTTIME <<  " min: " << errors[0] << " max: " << errors[1] 
	     << " " << errors[1] - errors[0] <<endl );
  }
  else
  {
      //UltraLocalError_TCD(SolArray[mg_level-1], errors+5, lump_mass);
      OutPut(TDatabase::TimeDB->CURRENTTIME <<  " min: " << errors[7] << " max: " << errors[8] 
	     << " " << errors[8] - errors[7] <<endl );
  }
 
  // errors in smooth region
  aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
  conc->GetErrors(Exact, 3, AllDerivatives, 3, SDFEMErrorsSmooth_JohnMaubachTobiska1997,
		  BilinearCoeffs, aux, 1, fesp, errors_smooth);
  delete aux;
  
  OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
  OutPut(" L2_sm: " << errors_smooth[0]);
  OutPut(" H1-semi_sm: " << errors_smooth[1]);
  OutPut(" SD_sm: " << errors_smooth[2] << endl);
  
  errors_smooth[4] += (errors_smooth[0]*errors_smooth[0] +olderror_smooth[0] * olderror_smooth[0])*
      TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
  olderror_smooth[0] = errors_smooth[0];
  OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2)_sm " << sqrt(errors_smooth[4]) << " ");

  errors_smooth[5] += (errors_smooth[1]*errors_smooth[1] +olderror_smooth[1] * olderror_smooth[1])*
      TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
  OutPut( "L2(0,T;H1)_sm " << sqrt(errors_smooth[5]));
  olderror_smooth[1] = errors_smooth[1];
  
  errors_smooth[6] += tau*(errors_smooth[2]*errors_smooth[2] +olderror_smooth[2] * olderror_smooth[2]);
  //*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
  OutPut( " SD(0,T)_sm " << sqrt(errors_smooth[6]) << endl);
  olderror_smooth[2] = errors_smooth[2];  
#endif
#ifdef   __HEMKER1996__
ComputeExtremalValues(N_Unknowns, sol,errors);
OutPut(setprecision(4) << "ext min " << errors[0] << "& max " << errors[1]<<endl);
#endif
    }                                             // endif MEASURE_ERRORS

#ifdef __IMPULSE__
 SolArray[mg_level-1]->FindGradient(1-1e-10,0.5,errors);
 OutPut(TDatabase::TimeDB->CURRENTTIME << " value at outlet center "<< errors[0] << endl); 
#endif

    if ((TDatabase::ParamDB->WRITE_GNU)
	||(TDatabase::ParamDB->WRITE_GMV)||(TDatabase::ParamDB->WRITE_MATLAB)
	||(TDatabase::ParamDB->WRITE_VTK))
    {
      picture = 0;
      if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      {
        if (TDatabase::ParamDB->WRITE_GNU)
        {
          os.seekp(std::ios::beg);
          os << GnuBaseName << m << ".gnu" << ends;
          Output->WriteGnuplot(os.str().c_str());
          picture++;
        }
        if(TDatabase::ParamDB->WRITE_GMV)
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
	
	if(TDatabase::ParamDB->WRITE_MATLAB)
	{
	    os.seekp(std::ios::beg);
	    os << MatlabBaseName << m << ".m" << ends;
	    Output->WriteMatlab(os.str().c_str());
	}

      }
    }
  }                                               // while
  //======================================================================
  // end of time cycle
  //======================================================================

  if (TDatabase::ParamDB->WRITE_GNU)
  {
    os.seekp(std::ios::beg);
    os << GnuBaseName <<  "end." << m << ".gnu" << ends;
    Output->WriteGnuplot(os.str().c_str());
  }

  if(TDatabase::ParamDB->WRITE_GMV)
  {
    os.seekp(std::ios::beg);
    os << GmvBaseName << "end." << m << ".gmv" << ends;
    Output->WriteGMV(os.str().c_str());
  }

  if(TDatabase::ParamDB->WRITE_VTK)
  {
      os.seekp(std::ios::beg);
      os << VtkBaseName << "end." << m << ".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
  }

  if(TDatabase::ParamDB->WRITE_MATLAB)
  {
    os.seekp(std::ios::beg);
    os << MatlabBaseName << "end." << m << ".m" << ends;
    Output->WriteMatlab(os.str().c_str());
  }

  t4 =  GetTime();
  total_time += t4 - t3;
  OutPut("total running time: " << total_time << endl);
  CloseFiles();
  return 0;
}
