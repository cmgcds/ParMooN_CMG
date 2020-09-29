// =======================================================================
//
// Purpose:     main program
//
// Author:      Volker Behns  22.07.97
//              Volker John   Jan. 2000
//
// =======================================================================

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

double bound = 0;

// ======================================================================
// utilities for main program
// ======================================================================
#include <MainUtilities.h>
#include <Upwind.h>

#include <NSE_MultiGrid.h>
#include <NSE_MGLevel1.h>
#include <NSE_MGLevel2.h>
#include <NSE_MGLevel3.h>
#include <NSE_MGLevel4.h>
#include <NSE_MGLevel14.h>
#include <ItMethod.h>
#include <FixedPointIte.h>
#include <MultiGridIte.h>
#include <FgmresIte.h>
#include <MultiGrid2D.h>
#include <MGLevel2D.h>
#include <MultiGridScaIte.h>
#include <RFB.h>

#define AMG 0
#define GMG 1
#define DIRECT 2

// =======================================================================
// include current example
// =======================================================================
//#include "Examples/NSE_2D/Linear.h"
// #include "Examples/NSE_2D/Const.h"
//#include "Examples/NSE_2D/FSHabil.h"
//#include "Examples/NSE_2D/FSHabil_slip.h"
// #include "Examples/NSE_2D/DC2.h"
//#include "Examples/NSE_2D/DrivenCavity.h"
#include "../Examples/NSE_2D/Benchmark.h"
//#include "Examples/NSE_2D/Benchmark_Neum.h"
// #include "Examples/NSE_2D/Frequence.h"
//#include "Examples/NSE_2D/Poiseuille.h"
// #include "Examples/NSE_2D/Poiseuille_Neum.h"
// #include "Examples/NSE_2D/Poiseuille2.h"
// #include "Examples/NSE_2D/Einfach.h"
//#include "Examples/NSE_2D/SinCos.h"
//#include "Examples/NSE_2D/BraessSarazin.h"
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
//#include "Examples/NSE_2D/SinCos.h"
//#include "Examples/NSE_2D/Calotte_test.h"
//#include "Examples/NSE_2D/LinkeNoFlow.h"
// #include "Examples/NSE_2D/ExpLayer.h"

int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *coll, *mortarcoll = NULL;
  TBaseCell *cell;
  TFESpace2D *velocity_space, *pressure_space, *streamfunction_space;
  TFESpace2D *velocity_space_low, *pressure_space_low;
  TFESpace2D *old_u_space, *old_p_space, *projection_space;
  TFESpace2D *pressure_separation_space;
  TFESpace2D **USpaces, **PSpaces, **PsiSpaces, **ProjectionSpaces;
  TOutput2D *Output;

  double *rhs, *sol, *oldsol, tol, tolerance, *psi, *defect, *fesol, *soldiff;
  double *rhs_low, *sol_low, *old_sol, *itmethod_sol, *itmethod_rhs;
  double *rhs_high, *nosep_p;
  int i,j,k,l, N_, Len, low, ii;
  int N_Rows, N_Columns, N_U, N_P, N_Unknowns, N_V;
  int N_Active, N_NonActive, N_U_low,  N_P_low, N_Unknowns_low;
  double *l2u1, *l2u2, *h1u1, *h1u2;
  double *l2p, *h1p, *sd, *l_inf, vect[3], exactvect[3];
  int which;
  double DiffL2, DiffH1;
  char *PRM, *GEO, *MAP;
  int LEVELS, BASELEVEL;
  int ret, pde;
  double negPower;
  double x,y,max,min,sum;
  double RE_NR;
  double tau1, tau2;
  double errors[4], p1, p2, errors_mg[4],velo_l2;
  double t1, t2, res, res2, oldres, solver_time,residual;
  double impuls_residual,limit, total_time1, total_time2;
  int N_LinIter;

  std::ostringstream os;
  char *PsBaseName, *GrapeBaseName, *GnuBaseName, *ReadGrapeBaseName;
  char *VtkBaseName, *MatlabBaseName, *GmvBaseName;

  double *val;
  TFEFunction2D *u1, *u2, *p, *fefct[3], *StreamFct;
  TFEFunction2D *u1_low, *u2_low, *p_low, *old_p, *soldiff_fe1,*soldiff_fe2;
  TFEFunction2D **U1Array, **U2Array, **AuxFEFunctArray;
  TFEFunction2D **PArray, *AuxPArray, *separated_pressure_fe_funct;
  TFEFunction2D *separated_pressure_rhs_fe_funct;
  TFEVectFunct2D *u, **UArray, *u_low, *old_u, **AuxFEVectFunctArray;
  TFESpace2D *fesp[2], *ferhs[3];

  TAuxParam2D *aux;

  TSquareStructure2D *sqstructureA, *sqstructurePressSep, *sqstructureC;
  TStructure2D *structureB, *structureBT;
  TSquareMatrix2D *sqmatrixA, *SQMATRICES[5];
  TSquareMatrix2D *sqmatrixA11, *sqmatrixA12;
  TSquareMatrix2D *sqmatrixA21, *sqmatrixA22;
  TSquareMatrix2D *sqmatrixC, **MatricesC;
  TSquareMatrix2D **MatricesA;
  TSquareMatrix2D **MatricesA11, **MatricesA12;
  TSquareMatrix2D **MatricesA21, **MatricesA22;
  TMatrix2D *matrixB1, *matrixB2, *MATRICES[6];
  TMatrix2D *matrixB1T, *matrixB2T;
  TMatrix2D **MatricesB1, **MatricesB2, **MatricesB1T, **MatricesB2T;
  TMatrix **matrices = (TMatrix **)MATRICES;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  MatVecProc *MatVect;
  DefectProc *Defect;
  TSquareMatrix2D *sqmatrixPressSep;

  TSquareStructure2D *sqstructureA_low, *sqstructureC_low;
  TStructure2D *structureB_low, *structureBT_low;
  TSquareMatrix2D *sqmatrixA_low;
  TSquareMatrix2D *sqmatrixA11_low, *sqmatrixA12_low;
  TSquareMatrix2D *sqmatrixA21_low, *sqmatrixA22_low;
  TSquareMatrix2D *sqmatrixC_low;
  TSquareMatrix2D **MatricesA_low;
  TSquareMatrix2D **MatricesA11_low, **MatricesA12_low;
  TSquareMatrix2D **MatricesA21_low, **MatricesA22_low;
  TMatrix2D *matrixB1_low, *matrixB2_low;
  TMatrix2D *matrixB1T_low, *matrixB2T_low;
  TMatrix2D **MatricesB1_low, **MatricesB2_low, **MatricesB1T_low, **MatricesB2T_low;

  TSquareStructure2D *sqstructureL;
  TStructure2D *structure_tilde_G, *structure_G;
  TSquareMatrix2D *sqmatrixL, **MatricesL;
  TMatrix2D *matrix_tilde_G11, *matrix_tilde_G22;
  TMatrix2D *matrix_G11, *matrix_G22;
  TMatrix2D **Matrices_tilde_G11, **Matrices_tilde_G22;
  TMatrix2D **Matrices_G11, **Matrices_G22;
  int N_L;

  int N_SquareMatrices, N_RectMatrices;
  int N_Rhs, N_FESpaces;

  double **RhsArray;

  TNSE_MGLevel *MGLevel, *MGLevel_low;
  TNSE_MultiGrid *MG;

  double *RHSs[3];
  int *N_Uarray, *N_Parray;

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

  BoundCondFunct2D *BoundaryConditions[3], *BoundaryConditionsPressureSeparation[1];
  BoundValueFunct2D *BoundValues[3], *BoundaryValuesPressureSeparation[1];
  CoeffFct2D *Coefficients[1];
  double average, hmin, hmax;

  TItMethod *itmethod, *prec, *Auxprec, *Auxitmethod;
  int Max_It, FirstSolve;
  double omega, alpha[2], alpha_fine[2], cd,cl,dp;
  int N_Parameters=2,n_aux, **downwind;
  double Parameters[4],delta0,delta1, reatt_pt,  reatt_point[3];
  double convergence_speed, residuals[10], firsterror,lasterror;
  double firsterrorl2,lasterrorl2,p3,p4;
  int last_digit_ite,slow_conv,last_sq,mg_level,mg_type, pre_calculation=1;
  int calculations=1,ll, zerostart;
  int velocity_space_code, pressure_space_code;
  int pressure_separation = 0, N_P_sep;
  double *separated_pressure_array, *separated_pressure_aux_array;
  double *pressure_aux_array;
  double *rhsPressSep;

  TMultiGrid2D *AuxMG;
  TMGLevel2D *AuxMGLevel;
  double *Auxitmethod_sol, *Auxitmethod_rhs;
  int number_pressep;

  char Readin[] = "readin.dat";
  char Name[] = "name";
  char UString[] = "u";
  char PString[] = "p";
  char PsiString[] = "psi";
  char DString[] = "d";
  char PsepString[] = "psep";

  #ifdef __BENCH__
  double Cd, Cl, dP1[3], dP2[3];
  #endif

  os << " ";

  //======================================================================
  // read parameter file
  //======================================================================
  if (argc >= 2)
    ret = Domain->ReadParam(argv[1]);
  else
    ret = Domain->ReadParam(Readin);

  if (ret == -1)
  {
    exit(-1);
  }

  RE_NR = TDatabase::ParamDB->RE_NR;
  delta0 = TDatabase::ParamDB->DELTA0;
  delta1 = TDatabase::ParamDB->DELTA1;
  convergence_speed = TDatabase::ParamDB->SC_DIV_FACTOR;

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->CheckParameterConsistencyNSE();
  Database->WriteParamDB(argv[0]);
  ExampleFile();

  //======================================================================
  // copy read parameters into local variables
  //======================================================================

  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
  MAP = TDatabase::ParamDB->MAPFILE;

  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  GrapeBaseName = TDatabase::ParamDB->GRAPEBASENAME;
  GnuBaseName = TDatabase::ParamDB->GNUBASENAME;
  ReadGrapeBaseName = TDatabase::ParamDB->READGRAPEBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  MatlabBaseName = TDatabase::ParamDB->MATLABBASENAME;
  GmvBaseName = TDatabase::ParamDB->GMVBASENAME;

  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SADDLE;
  if (mg_type)
    mg_level = 0;
  else
    mg_level = -1;
  LEVELS = TDatabase::ParamDB->LEVELS;
  BASELEVEL = TDatabase::ParamDB->UNIFORM_STEPS;
  l2u1 = new double[LEVELS+1];
  l2u2 = new double[LEVELS+1];
  l2p = new double[LEVELS+1];
  h1u1 = new double[LEVELS+1];
  h1u2 = new double[LEVELS+1];
  h1p = new double[LEVELS+1];
  sd = new double[LEVELS+1];
  l_inf = new double[LEVELS+1];

  U1Array = new TFEFunction2D*[LEVELS+1];
  U2Array = new TFEFunction2D*[LEVELS+1];
  PArray = new TFEFunction2D*[LEVELS+1];
  UArray = new TFEVectFunct2D*[LEVELS+1];

  RhsArray = new double* [LEVELS+1];
  N_Uarray = new int[LEVELS+1];
  N_Parray = new int[LEVELS+1];

  USpaces = new TFESpace2D*[LEVELS+1];
  PSpaces = new TFESpace2D*[LEVELS+1];
  PsiSpaces = new TFESpace2D*[LEVELS+1];

  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      MatricesA = new TSquareMatrix2D*[LEVELS+1];

      MatricesB1 = new TMatrix2D*[LEVELS+1];
      MatricesB2 = new TMatrix2D*[LEVELS+1];
      MatVect = MatVect_NSE1;
      Defect = Defect_NSE1;
      break;

    case 2:
      MatricesA = new TSquareMatrix2D*[LEVELS+1];

      MatricesB1 = new TMatrix2D*[LEVELS+1];
      MatricesB2 = new TMatrix2D*[LEVELS+1];
      MatricesB1T = new TMatrix2D*[LEVELS+1];
      MatricesB2T = new TMatrix2D*[LEVELS+1];
      MatVect = MatVect_NSE2;
      Defect = Defect_NSE2;
      break;

    case 3:
      MatricesA11 = new TSquareMatrix2D*[LEVELS+1];
      MatricesA12 = new TSquareMatrix2D*[LEVELS+1];
      MatricesA21 = new TSquareMatrix2D*[LEVELS+1];
      MatricesA22 = new TSquareMatrix2D*[LEVELS+1];

      MatricesB1 = new TMatrix2D*[LEVELS+1];
      MatricesB2 = new TMatrix2D*[LEVELS+1];
      MatVect = MatVect_NSE3;
      Defect = Defect_NSE3;
      break;

    case 4:
      MatricesA11 = new TSquareMatrix2D*[LEVELS+1];
      MatricesA12 = new TSquareMatrix2D*[LEVELS+1];
      MatricesA21 = new TSquareMatrix2D*[LEVELS+1];
      MatricesA22 = new TSquareMatrix2D*[LEVELS+1];

      MatricesB1 = new TMatrix2D*[LEVELS+1];
      MatricesB2 = new TMatrix2D*[LEVELS+1];
      MatricesB1T = new TMatrix2D*[LEVELS+1];
      MatricesB2T = new TMatrix2D*[LEVELS+1];
      MatVect = MatVect_NSE4;
      Defect = Defect_NSE4;
      break;

    case 14:
      MatricesA11 = new TSquareMatrix2D*[LEVELS+1];
      MatricesA12 = new TSquareMatrix2D*[LEVELS+1];
      MatricesA21 = new TSquareMatrix2D*[LEVELS+1];
      MatricesA22 = new TSquareMatrix2D*[LEVELS+1];
      MatricesC = new TSquareMatrix2D*[LEVELS+1];

      MatricesB1 = new TMatrix2D*[LEVELS+1];
      MatricesB2 = new TMatrix2D*[LEVELS+1];
      MatricesB1T = new TMatrix2D*[LEVELS+1];
      MatricesB2T = new TMatrix2D*[LEVELS+1];
      MatVect = MatVect_EquOrd_NSE4;
      Defect = Defect_EquOrd_NSE4;
      break;
  }                              // endswitch

  // matrices for VMS_PROJECTION
  if (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION)
  {
    MatricesL = new TSquareMatrix2D*[LEVELS+1];
    Matrices_tilde_G11 = new TMatrix2D*[LEVELS+1];
    Matrices_tilde_G22 = new TMatrix2D*[LEVELS+1];
    Matrices_G11 = new TMatrix2D*[LEVELS+1];
    Matrices_G22 = new TMatrix2D*[LEVELS+1];
    ProjectionSpaces = new TFESpace2D*[LEVELS+1];
  }
  downwind = new int*[LEVELS+1];

  //======================================================================
  // read boundary parameterization and initialize coarse grid
  //======================================================================
  Domain->Init(PRM, GEO);

  //======================================================================
  // do some special mortar stuff
  //======================================================================
  #ifdef __MORTAR__
  int N_Mortar;
  TVelocity_SpaceD *fespace_mortar;
  TStructure2D *struct_mortar;
  TMatrix2D *matrix_mortar;

  if (!strcmp(MAP, "NO_MAP_FILE"))
  {
    OutPut("switch off 'MORTAR = -D__MORTAR__' in the makefile !!!" << endl);
    OutPut("set 'MORTAR = ' in the makefile !!!" << endl);
    exit(1);
  }
  else
  {
    Domain->ReadMapFile(MAP, Database);
    //Domain->TestMortar();

    if (TDatabase::ParamDB->CONVERT_QUAD_TO_TRI > 0.5)
    {
      Domain->ConvertQuadToTri();
      OutPut("ConvertQuadToTri" << endl);
    }

    Domain->RegRefineSub(0);
    Domain->PS("Grid1.ps",It_Finest,0);
    //Domain->PS("Grid2.ps",It_Finest,0);
    //Domain->RegRefineAll();
  }
  #endif                         // __MORTAR__

  //======================================================================
  // initialize all discrete forms
  //======================================================================

  /*if( (TDatabase::ParamDB->DISCTYPE == 2) &&
    ((TDatabase::ParamDB->NSTYPE == 1) ||
    (TDatabase::ParamDB->NSTYPE == 3) ))
  {
    OutPut("DISCTYPE=2 works only with NSTYPE=2 or NSTYPE=4!" << endl);
    Error("DISCTYPE=2 works only with NSTYPE=2 or NSTYPE=4!" << endl);
    return -1;
  }
  */
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

  BoundaryConditions[0] = BoundCondition;
  BoundaryConditions[1] = BoundCondition;
  BoundaryConditions[2] = BoundConditionNoBoundCondition;

  BoundValues[0] = U1BoundValue;
  BoundValues[1] = U2BoundValue;
  BoundValues[2] = BoundaryValueNoBoundaryValue;

  BoundaryConditionsPressureSeparation[0] = BoundaryConditionPressSep;
  BoundaryValuesPressureSeparation[0] = BoundaryValuePressSep;

  Coefficients[0] = LinCoeffs;
  // refine up to user defined coarsest level

  for(i=0;i<TDatabase::ParamDB->SC_COARSEST_LEVEL_SADDLE;i++)
    Domain->RegRefineAll();

  // initialize solver parameters

  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
  omega = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SADDLE;
  velocity_space_code =   TDatabase::ParamDB->VELOCITY_SPACE;

  Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
  Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
  // Parameters[2] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;

  if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    i=1;
    MG = new TNSE_MultiGrid(i, N_Parameters, Parameters);
  }
  Parameters[2] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
  Parameters[3] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
  if ((TDatabase::ParamDB->PRESSURE_SEPARATION==3)||
    (TDatabase::ParamDB->PRESSURE_SEPARATION==4)||
    (TDatabase::ParamDB->PRESSURE_SEPARATION==5))
  {
    AuxMG = new TMultiGrid2D(i, N_Parameters, Parameters+2);
  }

  FirstSolve = TDatabase::ParamDB->SC_FIRST_SOLUTION_LEVEL_SADDLE;

  
  //======================================================================
  // loop over all levels
  //======================================================================

  for(i=0;i<LEVELS;i++)
  {
    mg_level++;
    total_time1 = GetTime();
    OutPut("*******************************************************" << endl);
    OutPut("******           GEOMETRY  LEVEL ");
    OutPut(i << "              *******" << endl);
    OutPut("******           MULTIGRID LEVEL ");
    OutPut(mg_level << "              *******" << endl);
    OutPut("*******************************************************" << endl);
    solver_time = 0.0;
    N_LinIter = 0;
    OutPut("memory before: " << setw(10) << GetMemory() << endl);
    for (j=0;j<10;j++)
      residuals[j]=1e10;
    slow_conv = 0;

    // refine grid if level is greater than 0
    if (i)
      Domain->RegRefineAll();

    coll=Domain->GetCollection(It_Finest, 0);
    Output = new TOutput2D(2, 2, 1, 1,Domain);
    cout << endl << endl;

    if(TDatabase::ParamDB->WRITE_PS)
    {
      // write grid into an Postscript file
      os.seekp(std::ios::beg);
      os << PsBaseName << i << ".ps" << ends;
      Domain->PS(os.str().c_str(),It_Finest,0);
    }

    #ifdef __MORTAR__
    mortarcoll = Domain->GetMortarColl(It_Mortar1, MAX_ItLevel);
    Domain->InitMortarJoints(It_Mortar1, MAX_ItLevel, mortarcoll);
    #endif

    // get spaces for low order disc on finest geo grid
    if (mg_type==1)
    {
      velocity_space_low = new TFESpace2D(coll,Name,UString,BoundCondition,
        Non_USpace,1, mortarcoll);
      pressure_space_low = new TFESpace2D(coll,Name,PString,BoundCondition,
        DiscP_PSpace,0, mortarcoll);
    }
    // get spaces of high order disc on finest geo grid
    if ((i>=FirstSolve)||(mg_type==0))
      GetVelocityAndPressureSpace(coll,BoundCondition,
        mortarcoll, velocity_space,
        pressure_space, &pressure_space_code,
        TDatabase::ParamDB->VELOCITY_SPACE,
        TDatabase::ParamDB->PRESSURE_SPACE);

    TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;

    // fe space for stream function
    streamfunction_space = new TFESpace2D(coll,Name,PsiString,BoundCondition,
      1, mortarcoll);

    #ifdef __MORTAR__
    fespace_mortar = new TVelocity_SpaceD(mortarcoll, "mortar space",
      "mortar space", velocity_space);

    N_Mortar =  fespace_mortar->GetN_DegreesOfFreedom();

    struct_mortar = new TStructure2D(fespace_mortar, velocity_space);
    matrix_mortar = new TMatrix2D(struct_mortar);
    #endif

    // build fespace hierarchy
    // set values and pointers for low order fe space
    if (mg_type==1)
    {
      USpaces[i] = velocity_space_low;
      PSpaces[i] = pressure_space_low;
      N_U_low = velocity_space_low->GetN_DegreesOfFreedom();
      N_U = N_U_low;
      N_P_low = pressure_space_low->GetN_DegreesOfFreedom();
      N_P = N_P_low;
      N_Uarray[i] = velocity_space_low->GetN_DegreesOfFreedom();
      N_Parray[i] = pressure_space_low->GetN_DegreesOfFreedom();
    }
    // set values and pointers for high order fe space
    if ((i>=FirstSolve)||(mg_type==0))
    {
      USpaces[mg_level] = velocity_space;
      PSpaces[mg_level] = pressure_space;
      N_U = velocity_space->GetN_DegreesOfFreedom();
      N_P = pressure_space->GetN_DegreesOfFreedom();
      N_Uarray[mg_level] = N_U;
      N_Parray[mg_level] = N_P;
      N_Active = velocity_space->GetActiveBound();
      N_NonActive = N_U - N_Active;
      PsiSpaces[i] = streamfunction_space;
      N_V = streamfunction_space->GetN_DegreesOfFreedom();
    }

    if (TDatabase::ParamDB->PRESSURE_SEPARATION>=1)
      //  ((i>=FirstSolve)||(mg_type==0)))
    {
      if (i<FirstSolve)
        pressure_space_code=0;
      OutPut("pressure_space_code " << pressure_space_code << endl);
      // allocate finite element space for separated pressure
      switch (pressure_space_code)
      {
        case 0:
          pressure_separation_space = new TFESpace2D(coll,Name,PsiString,BoundaryConditionPressSep,
            1, mortarcoll);
          break;
        case 1:
        case -11:
          pressure_separation_space = new TFESpace2D(coll,Name,PsiString,BoundConditionVMM,
            2, mortarcoll);
          break;
        default:
          OutPut("case for pressure_space_code not implemented" << endl);
          exit(4711);
      }
      N_P_sep = pressure_separation_space->GetN_DegreesOfFreedom();
      separated_pressure_array = new double[N_P_sep];
      memset(separated_pressure_array,0, N_P_sep*SizeOfDouble);
      separated_pressure_aux_array = new double[N_P_sep];
      memset(separated_pressure_aux_array,0, N_P_sep*SizeOfDouble);
      separated_pressure_fe_funct = new TFEFunction2D(pressure_separation_space,
        PsepString, PsepString,separated_pressure_array,
        N_P_sep);
      if (i<FirstSolve)
        N_P = N_P_sep;
      nosep_p = new double[N_P];
      if (TDatabase::ParamDB->PRESSURE_SEPARATION==1)
        pressure_separation = 0;
      else
        pressure_separation = 1;
      if ((TDatabase::ParamDB->PRESSURE_SEPARATION==3)||
        (TDatabase::ParamDB->PRESSURE_SEPARATION==4)||
        (TDatabase::ParamDB->PRESSURE_SEPARATION==5))
      {
        // allocate matrices
        sqstructurePressSep = new TSquareStructure2D(pressure_separation_space);
        sqmatrixPressSep = new TSquareMatrix2D(sqstructurePressSep);
        // rhs for auxiliary problem
        rhsPressSep = new double[N_P_sep];
        memset(rhsPressSep, 0, N_P_sep*SizeOfDouble);
        separated_pressure_rhs_fe_funct = new TFEFunction2D(pressure_separation_space,
          PsepString, PsepString,rhsPressSep,
          N_P_sep);
        // auxiliary array for prolongation
        pressure_aux_array = new double[N_P];
        memset(pressure_aux_array,0, N_P*SizeOfDouble);
        //OutPut("NP " << N_P << " NPSEP " << N_P_sep << endl);
        // allocate multigrid level
        n_aux = 4;
        AuxMGLevel = new TMGLevel2D(i, sqmatrixPressSep,
          rhsPressSep, separated_pressure_array,
          n_aux, NULL);
        AuxMG->AddLevel(AuxMGLevel);
        if (i <= FirstSolve)
        {
          // the matrix for the auxiliary problem must be build here
          // it is needed on the lower levels of the multigrid method
          OutPut("assemble pressure separation problem"<<endl);
          N_SquareMatrices = 1;
          N_RectMatrices = 0;
          N_Rhs = 1;
          SQMATRICES[0] = sqmatrixPressSep;
          SQMATRICES[0]->Reset();
          fesp[0] = pressure_separation_space;
          ferhs[0] = pressure_separation_space;
          // form of rhs on lower levels not important
          // use for for PRESSURE_SEPARATION 3
          number_pressep = TDatabase::ParamDB->PRESSURE_SEPARATION;
          TDatabase::ParamDB->PRESSURE_SEPARATION=3;
          N_FESpaces = 1;
          aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
          RHSs[0] = rhsPressSep;

          // assemble
          Assemble2D(N_FESpaces, fesp,
            N_SquareMatrices, SQMATRICES,
            0, NULL,
            N_Rhs, RHSs, ferhs,
            DiscreteFormAuxProbPressSep,
            BoundaryConditionsPressureSeparation,
            BoundaryValuesPressureSeparation,
            aux);
          delete aux;
          TDatabase::ParamDB->PRESSURE_SEPARATION =  number_pressep;
        }
      }
    }

    // build matrices for low order disc
    if (mg_type==1)
    {
	AllocateMatricesNSE_2D(i, velocity_space_low, pressure_space_low,
			       sqstructureA_low, sqstructureC_low, structureB_low,  structureBT_low,
			       sqmatrixA_low, 
			       sqmatrixA11_low, sqmatrixA12_low, sqmatrixA21_low, sqmatrixA22_low, 
			       sqmatrixC_low,
			       matrixB1_low, matrixB2_low, matrixB1T_low, matrixB2T_low, 
			       MatricesA, 
			       MatricesA11, MatricesA12, MatricesA21, MatricesA22,
			       MatricesC,
			       MatricesB1,  MatricesB2,  MatricesB1T,  MatricesB2T);
    }                            // end if (mg_type==1)
   
    // build matrices for high order disc
    if ((i>=FirstSolve)||(mg_type==0))
    {
	AllocateMatricesNSE_2D(mg_level, velocity_space, pressure_space,
			       sqstructureA, sqstructureC, structureB,  structureBT,
			       sqmatrixA, 
			       sqmatrixA11, sqmatrixA12, sqmatrixA21, sqmatrixA22,
			       sqmatrixC,
			       matrixB1, matrixB2, matrixB1T, matrixB2T, 
			       MatricesA, 
			       MatricesA11, MatricesA12, MatricesA21, MatricesA22,
			       MatricesC,
			       MatricesB1,  MatricesB2,  MatricesB1T,  MatricesB2T);
    }
   
    #ifdef __MORTAR__
    N_Unknowns = 2*N_U + N_P + 2* N_Mortar;
    OutPut("dof mortar   : "<< setw(10) << 2*N_Mortar << endl);
    #else
    N_Unknowns = 2*N_U + N_P;
    if (mg_type==1)
      N_Unknowns_low = 2*N_U_low + N_P_low;
    #endif

    coll->GetHminHmax(&hmin,&hmax);
    OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
    OutPut("dof velocity : "<< setw(10) << 2* N_U << endl);
    OutPut("dof pressure : "<< setw(10) << N_P << endl);
    OutPut("dof all      : "<<  setw(10) << N_Unknowns  << endl);

    // matrices for VMS_PROJECTION
    if ((i>=FirstSolve)||(mg_type==0))
    {
    if (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION)
    {
	switch(TDatabase::ParamDB->VMS_LARGE_VELOCITY_SPACE)
	{
	    case -1:
		projection_space = new TFESpace2D(coll, Name, UString, BoundCondition,
					      DiscP_PSpace,102, mortarcoll);
		break;
	    case 0:
		projection_space = new TFESpace2D(coll, Name, UString, BoundCondition,
					      DiscP_PSpace,0, mortarcoll);
		break;
	    case 1:
		projection_space = new TFESpace2D(coll, Name, UString, BoundCondition,
						  DiscP_PSpace,1, mortarcoll);
		break;
	    default:
		OutPut("VMS_LARGE_VELOCITY_SPACE: " <<
		       TDatabase::ParamDB->VMS_LARGE_VELOCITY_SPACE <<
		       " not available !!!" << endl);
		exit(4711);
	}

      ProjectionSpaces[mg_level] = projection_space;
      sqstructureL = new TSquareStructure2D(projection_space);
      sqstructureL->Sort();
      structure_tilde_G = new TStructure2D(velocity_space, projection_space);
      structure_G = new TStructure2D(projection_space, velocity_space);
      sqmatrixL = new TSquareMatrix2D(sqstructureL);
      MatricesL[mg_level] = sqmatrixL;
      matrix_tilde_G11 = new TMatrix2D(structure_tilde_G);
      Matrices_tilde_G11[mg_level] = matrix_tilde_G11;
      matrix_tilde_G22 = new TMatrix2D(structure_tilde_G);
      Matrices_tilde_G22[mg_level] = matrix_tilde_G22;
      matrix_G11 = new TMatrix2D(structure_G);
      Matrices_G11[mg_level] = matrix_G11;
      matrix_G22 = new TMatrix2D(structure_G);
      Matrices_G22[mg_level] = matrix_G22;
      N_L = projection_space->GetN_DegreesOfFreedom();
      OutPut("dof projection : " << setw(10) << N_L << endl);
    }
    }

    if (TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALE_SADDLE)
    {
      limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE*sqrt(1.0*N_Unknowns);
      OutPut("stopping tolerance for nonlinear iteration " << limit << endl);
    }
    if (mg_type==1)
    {
      OutPut("dof low order disc     : "<<  setw(10) << N_Unknowns_low  << endl);

      // initialize solver
      // low order disc
      rhs_low = new double[N_Unknowns_low];
      memset(rhs_low, 0, N_Unknowns_low*SizeOfDouble);
      RhsArray[i] = rhs_low;
      sol_low = new double[N_Unknowns_low];
      memset(sol_low, 0, N_Unknowns_low*SizeOfDouble);
    }

    // high order disc
    if ((i>=FirstSolve)||(mg_type==0))
    {
      rhs_high = new double[N_Unknowns];
      memset(rhs_high, 0, N_Unknowns*SizeOfDouble);
      RhsArray[mg_level] = rhs_high;
      sol = new double[N_Unknowns];
      oldsol = new double[N_Unknowns];
      memset(sol, 0, N_Unknowns*SizeOfDouble);
      memset(oldsol, 0, N_Unknowns*SizeOfDouble);
      /* for (k=0;k<N_U;k++)
        sol[k] = 5;
      for (k=0;k<N_U;k++)
      sol[N_U+k] = -4;*/
    }

    // build multigrid level(s)
    // ( A B' )
    // ( B 0  )
    switch(TDatabase::ParamDB->SOLVER_TYPE)
    {
      case AMG:
      case DIRECT:
        low = mg_level;
        break;

      case GMG:
        // coarsest grid number
        low = 0;
        // determine number of auxiliary arrays
        if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE)
          || (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE))
          n_aux=4;
        else
          n_aux=2;

        if (i==0)
        {
          alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
          alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
        }
        else
        {
          alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
          alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
        }
        if (mg_type==1)
        {
          alpha_fine[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
          alpha_fine[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_FINE_SADDLE;
        }
        else
        {
          alpha_fine[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
          alpha_fine[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
        }

        downwind[i] = new int[coll->GetN_Cells()];
        for (j=0;j<coll->GetN_Cells();j++)
          downwind[i][j] = j;
      #ifdef __DOWNWIND__
        DownwindNumberingCells(coll, downwind[i]);
      #endif
   
        // build fe multigrid levels
        switch(TDatabase::ParamDB->NSTYPE)
        {
          case 1:
            // low order disc
            if (mg_type==1)
            {
              MGLevel_low = new TNSE_MGLevel1(i, sqmatrixA_low,
                matrixB1_low, matrixB2_low,
                structureBT_low,
                rhs_low,  sol_low,
                n_aux, alpha, -1, 0, NULL,downwind[i] );
              if (i==0)
                MG->AddLevel(MGLevel_low);
              else
                MG->ReplaceLevel(i,MGLevel_low);
            }
            // high order disc
            if ((i>=FirstSolve)||(mg_type==0))
            {
              MGLevel = new TNSE_MGLevel1(mg_level, sqmatrixA,
                matrixB1, matrixB2,
                structureBT,
                rhs_high,  sol,  n_aux, alpha_fine,
                velocity_space_code, pressure_space_code, NULL,downwind[i]);
              MG->AddLevel(MGLevel);
            }
            break;

          case 2:
            if (mg_type==1)
            {
              MGLevel_low = new TNSE_MGLevel2(i, sqmatrixA_low,
                matrixB1_low, matrixB2_low,
                matrixB1T_low, matrixB2T_low,
                rhs_low, sol_low,
                n_aux, alpha, -1, 0, NULL,downwind[i]);
              if (i==0)
                MG->AddLevel(MGLevel_low);
              else
                MG->ReplaceLevel(i,MGLevel_low);
            }
            // high order disc
            if ((i>=FirstSolve)||(mg_type==0))
            {
              MGLevel = new TNSE_MGLevel2(mg_level, sqmatrixA,
                matrixB1, matrixB2,
                matrixB1T, matrixB2T,
                rhs_high, sol, n_aux, alpha_fine,
                velocity_space_code, pressure_space_code, NULL,downwind[i]);
              MG->AddLevel(MGLevel);
            }
            break;

          case 3:
            if (mg_type==1)
            {
              MGLevel_low = new TNSE_MGLevel3(i,
                sqmatrixA11_low, sqmatrixA12_low,
                sqmatrixA21_low, sqmatrixA22_low,
                matrixB1_low, matrixB2_low,
                structureBT_low,
                rhs_low, sol_low,
                n_aux, alpha, -1, 0, NULL,downwind[i]);
              if (i==0)
                MG->AddLevel(MGLevel_low);
              else
                MG->ReplaceLevel(i,MGLevel_low);
            }
            // high order disc
            if ((i>=FirstSolve)||(mg_type==0))
            {
              MGLevel = new TNSE_MGLevel3(mg_level,
                sqmatrixA11, sqmatrixA12,
                sqmatrixA21, sqmatrixA22, matrixB1, matrixB2,
                structureBT,
                rhs_high, sol, n_aux, alpha_fine,
                velocity_space_code, pressure_space_code, NULL,downwind[i]);
              MG->AddLevel(MGLevel);
            }
            break;

          case 4:
	      // low order disc
            if (mg_type==1)
            {
              MGLevel_low = new TNSE_MGLevel4(i,
                sqmatrixA11_low,  sqmatrixA12_low,
                sqmatrixA21_low, sqmatrixA22_low,
                matrixB1_low, matrixB2_low,
                matrixB1T_low, matrixB2T_low,
                rhs_low, sol_low,
                n_aux, alpha, -1, 0, NULL,downwind[i]);
              if (i==0)
                MG->AddLevel(MGLevel_low);
              else
                MG->ReplaceLevel(i,MGLevel_low);
	    }
            // high order disc
            if ((i>=FirstSolve)||(mg_type==0))
            {
              MGLevel = new TNSE_MGLevel4(mg_level, sqmatrixA11,  sqmatrixA12,
                sqmatrixA21, sqmatrixA22, matrixB1, matrixB2,
                matrixB1T, matrixB2T,
                rhs_high, sol, n_aux, alpha_fine,
                velocity_space_code, pressure_space_code, NULL,downwind[i]);
              MG->AddLevel(MGLevel);
            }
            break;
	    case 14:
	      // low order disc
            if (mg_type==1)
            {
              MGLevel_low = new TNSE_MGLevel14(i,
                sqmatrixA11_low,  sqmatrixA12_low,
                sqmatrixA21_low, sqmatrixA22_low,
		sqmatrixC_low,
                matrixB1_low, matrixB2_low,
                matrixB1T_low, matrixB2T_low,
                rhs_low, sol_low,
                n_aux, alpha, -1, 0, NULL,downwind[i]);
              if (i==0)
                MG->AddLevel(MGLevel_low);
              else
                MG->ReplaceLevel(i,MGLevel_low);
	    }
            // high order disc
            if ((i>=FirstSolve)||(mg_type==0))
            {
              MGLevel = new TNSE_MGLevel14(mg_level, sqmatrixA11,  sqmatrixA12,
                sqmatrixA21, sqmatrixA22, sqmatrixC, matrixB1, matrixB2,
                matrixB1T, matrixB2T,
                rhs_high, sol, n_aux, alpha_fine,
                velocity_space_code, pressure_space_code, NULL,downwind[i]);
              MG->AddLevel(MGLevel);
            }
            break;
       }                        // end switch(NSTYPE)
        break;
    }
   
    #ifdef __MORTAR__
    Assemble(matrix_mortar);
    #endif

    // build new fe functions
    // high order fe space
    if ((i>=FirstSolve)||(mg_type==0))
    {
      u = new TFEVectFunct2D(velocity_space, UString, UString, sol, N_U, 2);
      u1 = u->GetComponent(0);
      u2 = u->GetComponent(1);
      p = new TFEFunction2D(pressure_space, PString, PString, sol+2*N_U, N_P);

      U1Array[mg_level] = u1;
      U2Array[mg_level] = u2;
      PArray[mg_level] = p;
      UArray[mg_level] = u;
      
      /*  if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == OSEEN_PROBLEM)
      {
	  u1->Interpolate(InitialU1);
	  u2->Interpolate(InitialU2);
	  p->Interpolate(InitialP);
	  }*/
    }

    // low order fe space
    if (mg_type==1)
    {
      u_low = new TFEVectFunct2D(velocity_space_low, UString, UString, sol_low, N_U_low, 2);
      u1_low = u_low->GetComponent(0);
      u2_low = u_low->GetComponent(1);
      p_low = new TFEFunction2D(pressure_space_low, PString, PString, sol_low+2*N_U_low, N_P_low);
      U1Array[i] = u1_low;
      U2Array[i] = u2_low;
      PArray[i] = p_low;
      UArray[i] = u_low;
    }

                                 // CHECK THIS !!!
    if ((i>=FirstSolve)||(mg_type==0))
    {
      Output->AddFEVectFunct(u);
      Output->AddFEFunction(p);
    }

    // prolongation, to get a good starting iterate
    if(i && i>FirstSolve)
    {
      Prolongate(old_u_space, USpaces[mg_level],
        old_u->GetComponent(0)->GetValues(), U1Array[mg_level]->GetValues(),
        oldsol);
      Prolongate(old_u_space, USpaces[mg_level],
        old_u->GetComponent(1)->GetValues(), U2Array[mg_level]->GetValues(),
        oldsol);
      Prolongate(old_p_space, PSpaces[mg_level],
        old_p->GetValues(), PArray[mg_level]->GetValues(),
        oldsol);

      if ((TDatabase::ParamDB->SOLVER_TYPE==AMG)||(TDatabase::ParamDB->SOLVER_TYPE==DIRECT))
      {
        delete USpaces[i-1];
        delete PSpaces[i-1];
        delete U1Array[i-1]->GetValues();
      }
      // copy current solution for assembling the nonlinear terms
      memcpy(oldsol,sol, N_Unknowns*SizeOfDouble);

      if (mg_type==1)
      {
        delete old_sol;
        delete old_u;
        delete old_p;
        delete old_u_space;
        delete old_p_space;
      }
    }                            // end of prolongate

    if (TDatabase::ParamDB->P9==123456789)
    {
      for (k=0;k<2*N_U;k++)
        sol[k] = 1.0;
      for (k=0;k<N_P;k++)
        sol[2*N_U+k] = 0.0;
      fesol = new double[N_Unknowns];
      memset(fesol, 0, N_Unknowns*SizeOfDouble);
      soldiff = new double[N_Unknowns];
      memset(soldiff, 0, N_Unknowns*SizeOfDouble);
      soldiff_fe1 = new TFEFunction2D(velocity_space, DString, DString, soldiff,N_U);
      soldiff_fe2 = new TFEFunction2D(velocity_space, DString, DString, soldiff+N_U,N_U);
      pre_calculation = 1;
      calculations = 2;
      TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SADDLE = 1e-13;
    }

    // restrict solution to all grids
    if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
      MG->RestrictToAllGrids();

    // if no solution on this grid, continue
    if(FirstSolve>i)
      continue;

    if(TDatabase::ParamDB->READ_GRAPE_FILE)
    {
      AuxFEFunctArray = new TFEFunction2D*[2];
      AuxPArray =  new TFEFunction2D(velocity_space, PString, PString, oldsol, N_V);
      AuxFEFunctArray[0] = PArray[mg_level];
      AuxFEFunctArray[1] = AuxPArray;
      AuxFEVectFunctArray = new TFEVectFunct2D*[1];
      AuxFEVectFunctArray[0] = UArray[mg_level];
      ReadGrapeFile(ReadGrapeBaseName, 2, 1, AuxFEFunctArray,AuxFEVectFunctArray);
      delete AuxPArray;
      TDatabase::ParamDB->READ_GRAPE_FILE = 0;
      OutPut("u " << Ddot(2*N_U,sol,sol)<< endl);
      // for (ii=0;ii<N_Unknowns;ii++)
      //  OutPut(ii << " " << sol[ii] << endl);
      OutPut("p " << Ddot(N_P,sol+2*N_U,sol+2*N_U)<< endl);
      memcpy(oldsol,sol, N_Unknowns*SizeOfDouble);
      if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
        MG->RestrictToAllGrids();
    }

    // build the discretizations
    for(k=low;k<=mg_level;k++)
    {
      rhs = RhsArray[k];
      N_U = N_Uarray[k];
      N_P = N_Parray[k];
      N_Active = USpaces[k]->GetActiveBound();
      N_NonActive = N_U - N_Active;

      RHSs[0] = rhs;
      RHSs[1] = rhs + N_U;
      RHSs[2] = rhs + 2*N_U;
      memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);

      fesp[0] = USpaces[k];
      fesp[1] = PSpaces[k];

      fefct[0] = U1Array[k];
      fefct[1] = U2Array[k];
      ferhs[0] = USpaces[k];
      ferhs[1] = USpaces[k];
      ferhs[2] = PSpaces[k];

      // find discrete form
      if ((mg_type==1) && (k<i+1))
        DiscreteForm = DiscreteFormUpwind;
      else
        switch(TDatabase::ParamDB->DISCTYPE)
        {
          case GALERKIN:
	  case NSE_RFB:
            DiscreteForm = DiscreteFormGalerkin;
            break;

          case SDFEM:
            DiscreteForm = DiscreteFormSDFEM;
	    break;

          case UPWIND:
            DiscreteForm = DiscreteFormUpwind;
            break;

          case SMAGORINSKY:
            DiscreteForm = DiscreteFormSmagorinsky;
            break;

          case VMS_PROJECTION:
	      DiscreteForm = DiscreteFormVMSProjection;
	      if (TDatabase::ParamDB->NSTYPE != 1)
	    {
		OutPut("VMS only for NSTYPE 1 implemented !!!"<<endl);
		exit(4711);
	    }
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
          SQMATRICES[0] = MatricesA[k];
          MATRICES[0] = MatricesB1[k];
          MATRICES[1] = MatricesB2[k];

          SQMATRICES[0]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();

          N_SquareMatrices = 1;
          N_RectMatrices = 2;

          N_Rhs = 2;
          N_FESpaces = 2;
	  if (DiscreteForm == DiscreteFormVMSProjection)
	  {
	      SQMATRICES[1] = MatricesL[k];
	      MATRICES[2] = Matrices_tilde_G11[k];
	      MATRICES[3] = Matrices_tilde_G22[k];
	      MATRICES[4] = Matrices_G11[k];
	      MATRICES[5] = Matrices_G22[k];
	      
	      SQMATRICES[1]->Reset();
	      MATRICES[2]->Reset();
	      MATRICES[3]->Reset();
	      MATRICES[4]->Reset();
	      MATRICES[5]->Reset();
	      
	      N_SquareMatrices = 2;
	      N_RectMatrices = 6;
	      
	      N_FESpaces = 3;
	      fesp[2] = ProjectionSpaces[k];
	  }
          break;

        case 2:
          SQMATRICES[0] = MatricesA[k];
          MATRICES[0] = MatricesB1[k];
          MATRICES[1] = MatricesB2[k];
          MATRICES[2] = MatricesB1T[k];
          MATRICES[3] = MatricesB2T[k];

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
          SQMATRICES[0] = MatricesA11[k];
          SQMATRICES[1] = MatricesA12[k];
          SQMATRICES[2] = MatricesA21[k];
          SQMATRICES[3] = MatricesA22[k];
          MATRICES[0] = MatricesB1[k];
          MATRICES[1] = MatricesB2[k];

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
          SQMATRICES[0] = MatricesA11[k];
          SQMATRICES[1] = MatricesA12[k];
          SQMATRICES[2] = MatricesA21[k];
          SQMATRICES[3] = MatricesA22[k];
          MATRICES[0] = MatricesB1[k];
          MATRICES[1] = MatricesB2[k];
          MATRICES[2] = MatricesB1T[k];
          MATRICES[3] = MatricesB2T[k];

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
        case 14:
          SQMATRICES[0] = MatricesA11[k];
          SQMATRICES[1] = MatricesA12[k];
          SQMATRICES[2] = MatricesA21[k];
          SQMATRICES[3] = MatricesA22[k];
          SQMATRICES[4] = MatricesC[k];
          MATRICES[0] = MatricesB1[k];
          MATRICES[1] = MatricesB2[k];
          MATRICES[2] = MatricesB1T[k];
          MATRICES[3] = MatricesB2T[k];

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          SQMATRICES[4]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();

          N_SquareMatrices = 5;
          N_RectMatrices = 4;

          N_Rhs = 3;
          N_FESpaces = 2;
          break;
      }
      // get auxiliary values
      // fixed point iteration
      if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
	  if ((DiscreteForm == DiscreteFormSmagorinsky)||
	      (DiscreteForm == DiscreteFormVMSProjection))
	{
	            aux =  new TAuxParam2D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo,
			       NSN_ParamFctVelo_GradVelo,
			       NSN_FEValuesVelo_GradVelo,
			       fesp, fefct,
			       NSFctVelo_GradVelo,
			       NSFEFctIndexVelo_GradVelo, NSFEMultiIndexVelo_GradVelo,
			       NSN_ParamsVelo_GradVelo, NSBeginParamVelo_GradVelo);
	}
	else
	{
		aux =  new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo, NSN_ParamFctVelo,
				   NSN_FEValuesVelo,
				   fesp, fefct,
				   NSFctVelo,
				   NSFEFctIndexVelo, NSFEMultiIndexVelo,
				   NSN_ParamsVelo, NSBeginParamVelo);
	}
      else                       // Newton method
      {
	  aux =  new TAuxParam2D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo,
				 NSN_ParamFctVelo_GradVelo,
				 NSN_FEValuesVelo_GradVelo,
				 fesp, fefct,
				 NSFctVelo_GradVelo,
				 NSFEFctIndexVelo_GradVelo, NSFEMultiIndexVelo_GradVelo,
				 NSN_ParamsVelo_GradVelo, NSBeginParamVelo_GradVelo);
      }
  

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
            UpwindForNavierStokes(Coefficients[0], SQMATRICES[0], U1Array[k], U2Array[k]);
            cout << "UPWINDING DONE : level " << k << endl;
            break;

          case 3:
          case 4:
            // do upwinding with two matrices
            cout << "UPWINDING DONE : level " << k << endl;
            UpwindForNavierStokes(Coefficients[0], SQMATRICES[0], U1Array[k], U2Array[k]);
            UpwindForNavierStokes(Coefficients[0],SQMATRICES[3], U1Array[k], U2Array[k]);
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

        SQMATRICES[0] = MatricesA11[k];
        SQMATRICES[1] = MatricesA22[k];
        SQMATRICES[2] = MatricesA12[k];
        SQMATRICES[3] = MatricesA21[k];

        MATRICES[0] = MatricesB1T[k];
        MATRICES[1] = MatricesB2T[k];

        fesp[0] = USpaces[k];
        ferhs[0] = USpaces[k];
        ferhs[1] = USpaces[k];

        RHSs[0] = RhsArray[k];
        RHSs[1] = RhsArray[k]+N_Uarray[k];

        Assemble2DSlipBC(N_FESpaces, fesp,
          N_SquareMatrices, SQMATRICES,
          N_RectMatrices, MATRICES,
          N_Rhs, RHSs, ferhs,
          DiscreteForm0,
          BoundaryConditions,
          BoundValues,
          aux,
          U1Array[k],U2Array[k]);

        // reset MATRICES for solver
        SQMATRICES[0] = MatricesA11[k];
        SQMATRICES[1] = MatricesA12[k];
        SQMATRICES[2] = MatricesA21[k];
        SQMATRICES[3] = MatricesA22[k];
        MATRICES[0] = MatricesB1[k];
        MATRICES[1] = MatricesB2[k];
        MATRICES[2] = MatricesB1T[k];
        MATRICES[3] = MatricesB2T[k];
	}
      delete aux;

      // pressure separation with Neumann problem
      if (((TDatabase::ParamDB->PRESSURE_SEPARATION==3)||
        (TDatabase::ParamDB->PRESSURE_SEPARATION==4))
        &&(k==mg_level))
      {
        OutPut("assemble pressure separation problem"<<endl);
        N_SquareMatrices = 1;
        N_RectMatrices = 0;
        N_Rhs = 1;
        SQMATRICES[0] = sqmatrixPressSep;
        SQMATRICES[0]->Reset();
        fesp[0] = pressure_separation_space;
        ferhs[0] = pressure_separation_space;
        if (TDatabase::ParamDB->PRESSURE_SEPARATION==3)
        {
          N_FESpaces = 1;
          aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
        }
        else
        {
          N_FESpaces = 2;
          fesp[1] = USpaces[mg_level];
          fefct[0] = U1Array[mg_level];
          fefct[1] = U2Array[mg_level];
          aux =  new TAuxParam2D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo,
            NSN_ParamFctVelo_GradVelo,
            NSN_FEValuesVelo_GradVelo,
            fesp+1, fefct,
            NSFctVelo_GradVelo,
            NSFEFctIndexVelo_GradVelo, NSFEMultiIndexVelo_GradVelo,
            NSN_ParamsVelo_GradVelo, NSBeginParamVelo_GradVelo);
        }
        RHSs[0] = rhsPressSep;
        memset(rhsPressSep, 0, N_P_sep*SizeOfDouble);
        // assemble
        Assemble2D(N_FESpaces, fesp,
          N_SquareMatrices, SQMATRICES,
          0, NULL,
          N_Rhs, RHSs, ferhs,
          DiscreteFormAuxProbPressSep,
          BoundaryConditionsPressureSeparation,
          BoundaryValuesPressureSeparation,
          aux);

        // solve linear system
        //        Solver(sqmatrixPressSep, RHSs[0],separated_pressure_array,1);
        Auxprec = new TMultiGridScaIte(MatVect_Scalar, Defect_Scalar, NULL,
          0, N_P_sep, AuxMG, 0);
        Auxitmethod = new TFgmresIte(MatVect_Scalar, Defect_Scalar, Auxprec,
          0, N_P_sep, 1);
        Auxitmethod_sol = new double[N_P_sep];
        Auxitmethod_rhs = new double[N_P_sep];
        memcpy(Auxitmethod_sol, separated_pressure_array, N_P_sep*SizeOfDouble);
        memcpy(Auxitmethod_rhs, rhsPressSep, N_P_sep*SizeOfDouble);
        SQMATRICES[0] = sqmatrixPressSep;
        Auxitmethod->Iterate(sqmatrices,NULL,Auxitmethod_sol,Auxitmethod_rhs);
        memcpy(separated_pressure_array, Auxitmethod_sol, N_P_sep*SizeOfDouble);
        memcpy(rhsPressSep, Auxitmethod_rhs, N_P_sep*SizeOfDouble);
        delete Auxitmethod_sol;
        delete Auxitmethod_rhs;
        delete Auxprec;
        delete Auxitmethod;

        // store separated pressure in the original pressure space
        Prolongate(pressure_separation_space, pressure_space,
          separated_pressure_array, nosep_p,
          pressure_aux_array);

        // assemble rhs for NSE
        // the gradient of the separated pressure is needed for assembling
        // this has to be said to the assembling routine by an aux object
        fesp[0] = USpaces[mg_level];
        fesp[1] = pressure_separation_space;

        fefct[0] = separated_pressure_fe_funct;

        aux =  new TAuxParam2D(NSN_FESpacesPressSep, NSN_FctPressSep, NSN_ParamFctPressSep,
          NSN_FEValuesPressSep,
          fesp+1, fefct,
          NSFctPressSep,
          NSFEFctIndexPressSep, NSFEMultiIndexPressSep,
          NSN_ParamsPressSep, NSBeginParamPressSep);

        // assemble the right hand side
        N_FESpaces = 2;
        N_SquareMatrices = 0;
        N_RectMatrices = 0;
        N_Rhs = 2;
        RHSs[0] = rhs;
        RHSs[1] = rhs + N_U;
        memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);
        ferhs[0] = USpaces[mg_level];
        ferhs[1] = USpaces[mg_level];
        DiscreteForm0 = DiscreteFormPressSep;

        Assemble2D(N_FESpaces, fesp,
          N_SquareMatrices, SQMATRICES,
          N_RectMatrices, MATRICES,
          N_Rhs, RHSs, ferhs,
          DiscreteForm0,
          BoundaryConditions,
          BoundValues,
          aux);

        // initialize solution array for separated pressure
        memset(sol+2*N_U,0, N_P*SizeOfDouble);
      }
	
      if (DiscreteForm == DiscreteFormVMSProjection)
      {
	  SQMATRICES[0] = MatricesA[k];
	  SQMATRICES[1] = MatricesL[k];
	  MATRICES[2] = Matrices_tilde_G11[k];
	  MATRICES[3] = Matrices_tilde_G22[k];
	  MATRICES[4] = Matrices_G11[k];
	  MATRICES[5] = Matrices_G22[k];

	  LumpMassMatrixToDiag(MatricesL[k]);
	  VMSProjectionUpdateMatrices(N_Uarray[k],USpaces[k]->GetActiveBound(),
				      ProjectionSpaces[k]->GetN_DegreesOfFreedom(),
				      SQMATRICES,MATRICES);
      }
    }                            // endfor, assembling done

    // modify rhs for RFB stabilization
    if (TDatabase::ParamDB->DISCTYPE==NSE_RFB)
    {
	//ApproximateRFBSolutionQuad_Q2_NSE2D(coll, U1Array[mg_level], U2Array[mg_level],
	ApproximateRFBSolutionQuadNSE2D(coll, U1Array[mg_level], U2Array[mg_level],
					Coefficients[0], rhs);
	OutPut("RFB DONE"<<endl);
    }

    // set Dirichlet nodes
    memcpy(sol+N_Active, rhs_high+N_Active, N_NonActive*SizeOfDouble);
    memcpy(sol+N_U+N_Active, rhs_high+N_U+N_Active, N_NonActive*SizeOfDouble);

    // compute defect
    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
      IntoL20FEFunction(sol+2*N_U, N_P,PSpaces[mg_level],
        velocity_space_code, pressure_space_code);
    defect = new double[N_Unknowns];
    memset(defect,0,N_Unknowns*SizeOfDouble);
    Defect(sqmatrices,matrices,sol,rhs_high,defect);

    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
      IntoL20Vector2D(defect+2*N_U, N_P,pressure_space_code);
    residual =  Ddot(N_Unknowns,defect,defect);
    impuls_residual = Ddot(2*N_U,defect,defect);
    OutPut("nonlinear iteration step   0");
    OutPut(setw(14) << impuls_residual);
    OutPut(setw(14) << residual-impuls_residual);
    OutPut(setw(14) << sqrt(residual) << endl);

    // solve system
    switch(TDatabase::ParamDB->SOLVER_TYPE)
    {
      case AMG:
        TDatabase::ParamDB->SC_VERBOSE=1;
        TDatabase::ParamDB->CC_VERBOSE=1;
        t1 = GetTime();
        switch(TDatabase::ParamDB->NSTYPE)
        {
          case 1: 
            Solver(sqmatrixA, matrixB1, matrixB2, rhs_high, sol);
            break;

          case 2:
        #ifdef __MORTAR__
            Solver(sqmatrixA, matrixB1T, matrixB2T,
              matrixB1, matrixB2, matrix_mortar, rhs_high, sol);
        #else
            Solver(sqmatrixA, matrixB1T, matrixB2T,
              matrixB1, matrixB2, rhs_high, sol);
        #endif
            break;

          case 3:
            Solver(sqmatrixA11, sqmatrixA12, sqmatrixA21,
              sqmatrixA22, matrixB1, matrixB2, rhs_high, sol);
            break;

          case 4:
            Solver(sqmatrixA11, sqmatrixA12, sqmatrixA21,
              sqmatrixA22, matrixB1T, matrixB2T,
              matrixB1, matrixB2, rhs_high, sol);
            break;
	    default:
		OutPut("AMG solver for NSTYPE " << TDatabase::ParamDB->NSTYPE 
		       << " not implemented !!!" << endl);
		exit(4711);

        }
        t2 = GetTime();
        solver_time += (t2-t1);
        break;
      case DIRECT:
        t1 = GetTime();
        switch(TDatabase::ParamDB->NSTYPE)
        {
	    case 1: 
		DirectSolver(sqmatrixA, matrixB1, matrixB2, rhs_high, sol);
            break;

	    case 2:
		DirectSolver(sqmatrixA, matrixB1T, matrixB2T,
			     matrixB1, matrixB2, rhs_high, sol);
		break;
	    case 4:
		DirectSolver(sqmatrixA11, sqmatrixA12, sqmatrixA21,
			     sqmatrixA22, matrixB1T, matrixB2T,
			     matrixB1, matrixB2, rhs_high, sol);
		break;

	    case 14:
		DirectSolver(sqmatrixA11, sqmatrixA12, sqmatrixA21,
			     sqmatrixA22, sqmatrixC, matrixB1T, matrixB2T,
			     matrixB1, matrixB2, rhs_high, sol);
		break;
	    default:
		OutPut("Direct solver for NSTYPE " << TDatabase::ParamDB->NSTYPE
		       << " not implemented !!!" << endl);
		exit(4711);
        }
        t2 = GetTime();
        solver_time += (t2-t1);
	break;

      case GMG:
        t1 = GetTime();
        switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
        {
          case 11:
            zerostart = 1;
            break;
          case 16:
            zerostart = 0;
            break;
        }
        // build preconditioner
        switch (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE)
        {
          case 5:
            prec = new TMultiGridIte(MatVect, Defect, NULL,
              0, N_Unknowns, MG, zerostart);
            break;
          default:
            OutPut("Unknown preconditioner !!!" << endl);
            exit(4711);
        }
        switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
        {
          case 11:
            itmethod = new TFixedPointIte(MatVect, Defect, prec,
              0, N_Unknowns, 0);
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
            {
              itmethod_sol = new double[N_Unknowns];
              itmethod_rhs = new double[N_Unknowns];
              memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
              memcpy(itmethod_rhs, rhs_high, N_Unknowns*SizeOfDouble);
            }
            else
            {
              itmethod_sol = sol;
              itmethod_rhs = rhs_high;
            }
            break;
          case 16:
            itmethod = new TFgmresIte(MatVect, Defect, prec,
              0, N_Unknowns, 0);
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
            {
              itmethod_sol = new double[N_Unknowns];
              itmethod_rhs = new double[N_Unknowns];
              memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
              memcpy(itmethod_rhs, rhs_high, N_Unknowns*SizeOfDouble);
            }
            else
            {
              itmethod_sol = sol;
              itmethod_rhs = rhs_high;
            }
            break;
          default:
            OutPut("Unknown solver !!!" << endl);
            exit(4711);
        }
        t2 = GetTime();
        solver_time += (t2-t1);
        for (ll=0;ll<calculations;ll++)
        {
          if(TDatabase::ParamDB->P9 == 123456789)
          {
            fesp[0] = USpaces[mg_level];
            fefct[0] = U1Array[mg_level];
            fefct[1] = U2Array[mg_level];
            aux =  new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo,
              NSN_ParamFctVelo,
              NSN_FEValuesVelo,
              fesp, fefct,
              NSFctVelo,
              NSFEFctIndexVelo, NSFEMultiIndexVelo,
              NSN_ParamsVelo, NSBeginParamVelo);

            // errors in first velocity component
            U1Array[mg_level]->GetErrors(ExactU1, 3, NSAllDerivatives, 2,
              L2H1Errors,
              NULL, aux, 1, USpaces+mg_level, errors_mg);
            p1 = errors_mg[0] * errors_mg[0];

            // errors in second velocity component
            U2Array[mg_level]->GetErrors(ExactU2, 3, NSAllDerivatives, 2,
              L2H1Errors,
              NULL, aux, 1, USpaces+mg_level, errors_mg);
            p1 += errors_mg[0] * errors_mg[0];
            // errors in first velocity component
            //          soldiff_fe1->GetErrors(ExactNull, 3, NSAllDerivatives, 2,
            //                               L2H1Errors,
            //                     NULL, aux, 1, USpaces+mg_level, errors_mg);
            p1 = errors_mg[0] * errors_mg[0];

            // errors in second velocity component
            //soldiff_fe2->GetErrors(ExactNull, 3, NSAllDerivatives, 2,
            //                     L2H1Errors,
            //                     NULL, aux, 1, USpaces+mg_level, errors_mg);
            p1 += errors_mg[0] * errors_mg[0];
            firsterror =  p1 = sqrt(p1);
            delete aux;

            for(l=0;l<N_Unknowns;l++)
              soldiff[l] = sol[l]-fesol[l];

            p3 = sqrt(Ddot(2*N_U,soldiff,soldiff));
            firsterrorl2 = p3;

            OutPut("iteration -1  L2(u): " <<  p1 << " l2(u) " << p3 << endl);
            p2 = p1;
            p4 = p3;
          }                      // endif MEASURE_ERRORS

          t1 = GetTime();
          // solve linear system
          N_LinIter+=itmethod->Iterate(sqmatrices,matrices,itmethod_sol,itmethod_rhs);
          switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
          {
            case 11:
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
              {
                memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
                memcpy(rhs_high, itmethod_rhs, N_Unknowns*SizeOfDouble);
              }
              break;
            case 16:
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
              {
                memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
                memcpy(rhs_high, itmethod_rhs, N_Unknowns*SizeOfDouble);
              }
              break;
          }
          t2 = GetTime();
          solver_time += (t2-t1);

          if(TDatabase::ParamDB->P9 == 123456789)
          {
            fesp[0] = USpaces[mg_level];
            fefct[0] = U1Array[mg_level];
            fefct[1] = U2Array[mg_level];
            aux =  new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo,
              NSN_ParamFctVelo,
              NSN_FEValuesVelo,
              fesp, fefct,
              NSFctVelo,
              NSFEFctIndexVelo, NSFEMultiIndexVelo,
              NSN_ParamsVelo, NSBeginParamVelo);

            // errors in first velocity component
            U1Array[mg_level]->GetErrors(ExactU1, 3, NSAllDerivatives, 2,
              L2H1Errors,
              NULL, aux, 1, USpaces+mg_level, errors_mg);
            p1 = errors_mg[0] * errors_mg[0];

            // errors in second velocity component
            U2Array[mg_level]->GetErrors(ExactU2, 3, NSAllDerivatives, 2,
              L2H1Errors,
              NULL, aux, 1, USpaces+mg_level, errors_mg);
            // errors in first velocity component
            /* soldiff_fe1->GetErrors(ExactNull, 3, NSAllDerivatives, 2,
               L2H1Errors,
               NULL, aux, 1, USpaces+mg_level, errors_mg);
               p1 = errors_mg[0] * errors_mg[0];

               // errors in second velocity component
               soldiff_fe2->GetErrors(ExactNull, 3, NSAllDerivatives, 2,
               L2H1Errors,
               NULL, aux, 1, USpaces+mg_level, errors_mg);
            */
            p1 += errors_mg[0] * errors_mg[0];
            p1 = sqrt(p1);
            delete aux;

            for(l=0;l<N_Unknowns;l++)
              soldiff[l] = sol[l]-fesol[l];
            p3 = sqrt(Ddot(2*N_U,soldiff,soldiff));

            OutPut("iteration " << j << " L2(u): " <<  p1 << " redu ");
            OutPut(" rateL2 " << pow(p1/firsterror,1.0/N_LinIter));
            OutPut(p1/p2 << " l2(u) " << p3 << " redu " << p3/p4);
            OutPut(" ratel2 " << pow(p3/firsterrorl2,1.0/N_LinIter)<< endl);

            p2 = p1;
            p4 = p3;
            lasterror = p1;
            lasterrorl2 = p3;
            // if (res/res2 > convergence_speed)
            // {
            //  OutPut("SLOW !!!!!!!!! " << endl);
            //  break;
            // }
          }                      // endif MEASURE_ERRORS

          for(l=0;l<N_Unknowns;l++)
          {
            p2 = sol[l]-oldsol[l];
            sol[l] = oldsol[l] + omega * p2;
          }

          if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
            IntoL20FEFunction(sol+2*N_U, N_P,PSpaces[mg_level],
              velocity_space_code, pressure_space_code);

          if(TDatabase::ParamDB->P9 == 123456789)
          {
            if (!pre_calculation)
            {
              OutPut("average error reduction rate (L2/l2) " << pow(lasterror/firsterror,1.0/N_LinIter));
              OutPut(" " << pow(lasterrorl2/firsterrorl2,1.0/N_LinIter) << endl);
            }
            if (pre_calculation)
            {
              for (k=0;k<2*N_U+N_P;k++)
              {
                fesol[k] = sol[k];
                sol[k] = 0.0;
              }
              pre_calculation = 0;
              N_LinIter=0;
              TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SADDLE = 1e-10;
            }
          }
        }
        break;
    }

    // ************************************************************* //
    // end of first nonlinear step
    // ************************************************************* //
    // don't know why it does not work
    // if (TDatabase::ParamDB->PRESSURE_SEPARATION<3)
    {
      OutPut("MEMORY: " << setw(10) << GetMemory()/(1048576.0));
      OutPut(" MB" << endl);
    }
    // ************************************************************* //
    // the nonlinear iteration
    // ************************************************************* //

    for(j=1;j<=Max_It;j++)
    {
      if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
        MG->RestrictToAllGrids();

      memcpy(oldsol, sol, SizeOfDouble*N_Unknowns);

      for(k=low;k<=mg_level;k++)
      {
        fesp[0] = USpaces[k];
        fesp[1] = PSpaces[k];

        fefct[0] = U1Array[k];
        fefct[1] = U2Array[k];

        if ((k<i+1)&&(mg_type==1))
          DiscreteForm = DiscreteFormNLUpwind;
        else
          switch(TDatabase::ParamDB->DISCTYPE)
          {
            case GALERKIN:
	    case NSE_RFB:
              DiscreteForm = DiscreteFormNLGalerkin;
              break;

            case SDFEM:
		DiscreteForm = DiscreteFormNLSDFEM;
              break;

            case UPWIND:
              DiscreteForm = DiscreteFormNLUpwind;
              break;

            case SMAGORINSKY:
              DiscreteForm = DiscreteFormNLSmagorinsky;
              break;
	      
	      case VMS_PROJECTION:
		  DiscreteForm = DiscreteFormNLVMSProjection;
		  break;
          }                      // endswitch
        // this can only happen if pressure separation is applied
        if (TDatabase::ParamDB->STOKES_PROBLEM)
          DiscreteForm = DiscreteFormNLUpwind;

        switch(TDatabase::ParamDB->NSTYPE)
        {
          case 1:
            SQMATRICES[0] = MatricesA[k];
            SQMATRICES[0]->Reset();

            N_SquareMatrices = 1;
            N_RectMatrices = 0;

            N_Rhs = 0;
            N_FESpaces = 1;
            if (DiscreteForm == DiscreteFormNLSDFEM)
            {
              N_Rhs = 2;
              rhs = RhsArray[k];
              RHSs[0] = rhs;
              RHSs[1] = rhs + N_Uarray[k];
              memset(rhs, 0, (2*N_Uarray[k]+N_Parray[k])*SizeOfDouble);
              N_FESpaces = 2;
              ferhs[0] = USpaces[k];
              ferhs[1] = USpaces[k];
            }
            if (DiscreteForm == DiscreteFormNLVMSProjection)
            {
              MATRICES[0] =  Matrices_tilde_G11[k];
              MATRICES[1] =  Matrices_tilde_G22[k];

              MATRICES[0]->Reset();
              MATRICES[1]->Reset();
	      N_RectMatrices = 2;

	      N_FESpaces = 2;
	      fesp[1] = ProjectionSpaces[k];
	    }
            break;

          case 2:
            SQMATRICES[0] = MatricesA[k];
            SQMATRICES[0]->Reset();

            N_SquareMatrices = 1;
            if (DiscreteForm == DiscreteFormNLSDFEM)
            {
              N_RectMatrices = 2;
              MATRICES[0] = MatricesB1T[k];
              MATRICES[1] = MatricesB2T[k];

              MATRICES[0]->Reset();
              MATRICES[1]->Reset();

              N_Rhs = 2;
              rhs = RhsArray[k];
              RHSs[0] = rhs;
              RHSs[1] = rhs + N_Uarray[k];
              memset(rhs, 0, (2*N_Uarray[k]+N_Parray[k])*SizeOfDouble);
              N_FESpaces = 2;
              ferhs[0] = USpaces[k];
              ferhs[1] = USpaces[k];
            }
            else
            {
              N_RectMatrices = 0;

              N_Rhs = 0;
              N_FESpaces = 1;
            }
            break;

          case 3:
            if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
            {
              if(DiscreteForm == DiscreteFormNLSmagorinsky)
              {
                SQMATRICES[0] = MatricesA11[k];
                SQMATRICES[1] = MatricesA12[k];
                SQMATRICES[2] = MatricesA21[k];
                SQMATRICES[3] = MatricesA22[k];
                SQMATRICES[0]->Reset();
                SQMATRICES[1]->Reset();
                SQMATRICES[2]->Reset();
                SQMATRICES[3]->Reset();

                N_SquareMatrices = 4;
                N_RectMatrices = 0;
                last_sq = 3;
              }
              else
              {
                SQMATRICES[0] = MatricesA11[k];
                SQMATRICES[1] = MatricesA22[k];
                SQMATRICES[0]->Reset();
                SQMATRICES[1]->Reset();

                N_SquareMatrices = 2;
                N_RectMatrices = 0;
                last_sq = 1;
              }
              N_Rhs = 0;
              N_FESpaces = 1;
            }
            else                 // Newton method
            {
              SQMATRICES[0] = MatricesA11[k];
              SQMATRICES[1] = MatricesA12[k];
              SQMATRICES[2] = MatricesA21[k];
              SQMATRICES[3] = MatricesA22[k];
              SQMATRICES[0]->Reset();
              SQMATRICES[1]->Reset();
              SQMATRICES[2]->Reset();
              SQMATRICES[3]->Reset();

              N_SquareMatrices = 4;
              N_RectMatrices = 0;

              N_Rhs = 2;
              N_FESpaces = 1;
              rhs = RhsArray[k];
              RHSs[0] = rhs;
              RHSs[1] = rhs + N_Uarray[k];
              memset(rhs, 0, (2*N_Uarray[k]+N_Parray[k])*SizeOfDouble);
              ferhs[0] = USpaces[k];
              ferhs[1] = USpaces[k];
              last_sq = 3;
            }
            break;

          case 4:
            if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
            {
              if (DiscreteForm == DiscreteFormNLSDFEM)
              {
                N_SquareMatrices = 2;
                SQMATRICES[0] = MatricesA11[k];
                SQMATRICES[1] = MatricesA22[k];
                SQMATRICES[0]->Reset();
                SQMATRICES[1]->Reset();

                N_RectMatrices = 2;
                MATRICES[0] = MatricesB1T[k];
                MATRICES[1] = MatricesB2T[k];
                MATRICES[0]->Reset();
                MATRICES[1]->Reset();

                N_Rhs = 2;
                rhs = RhsArray[k];
                RHSs[0] = rhs;
                RHSs[1] = rhs + N_Uarray[k];
                memset(rhs, 0, (2*N_Uarray[k]+N_Parray[k])*SizeOfDouble);
                N_FESpaces = 2;

                ferhs[0] = USpaces[k];
                ferhs[1] = USpaces[k];
                last_sq = 1;
              }
              else
              {
                if (DiscreteForm == DiscreteFormNLSmagorinsky)
                {
                  N_RectMatrices = 0;

                  N_SquareMatrices = 4;
                  SQMATRICES[0] = MatricesA11[k];
                  SQMATRICES[1] = MatricesA12[k];
                  SQMATRICES[2] = MatricesA21[k];
                  SQMATRICES[3] = MatricesA22[k];
                  SQMATRICES[0]->Reset();
                  SQMATRICES[1]->Reset();
                  SQMATRICES[2]->Reset();
                  SQMATRICES[3]->Reset();

                  N_Rhs = 0;
                  N_FESpaces = 1;
                  last_sq = 3;
                }
                else             // default
                {
                  N_SquareMatrices = 2;
                  SQMATRICES[0] = MatricesA11[k];
                  SQMATRICES[1] = MatricesA22[k];
                  SQMATRICES[0]->Reset();
                  SQMATRICES[1]->Reset();

                  N_RectMatrices = 0;

                  N_Rhs = 0;
                  N_FESpaces = 1;
                  last_sq = 1;
                }
              }
	      break;
		case 14:
		    SQMATRICES[0] = MatricesA11[k];
		    SQMATRICES[1] = MatricesA12[k];
		    SQMATRICES[2] = MatricesA21[k];
		    SQMATRICES[3] = MatricesA22[k];
		    SQMATRICES[4] = MatricesC[k];
		    MATRICES[0] = MatricesB1[k];
		    MATRICES[1] = MatricesB2[k];
		    MATRICES[2] = MatricesB1T[k];
		    MATRICES[3] = MatricesB2T[k];

		    SQMATRICES[0]->Reset();
		    SQMATRICES[1]->Reset();
		    SQMATRICES[2]->Reset();
		    SQMATRICES[3]->Reset();
		    SQMATRICES[4]->Reset();
		    
		    MATRICES[0]->Reset();
		    MATRICES[1]->Reset();
		    MATRICES[2]->Reset();
		    MATRICES[3]->Reset();

		    N_SquareMatrices = 5;
		    N_RectMatrices = 4;
		    
		    N_Rhs = 3;
		    N_FESpaces = 2;
		    rhs = RhsArray[k];
		    RHSs[0] = rhs;
		    RHSs[1] = rhs + N_Uarray[k];
		    RHSs[2] = rhs + 2*N_Uarray[k];
		    memset(rhs, 0, (2*N_Uarray[k]+N_Parray[k])*SizeOfDouble);
		    ferhs[0] = USpaces[k];
		    ferhs[1] = USpaces[k];
		    ferhs[2] = PSpaces[k];

		    last_sq = 3;
		    break;
            }
            else                 // Newton method
            {
              SQMATRICES[0] = MatricesA11[k];
              SQMATRICES[1] = MatricesA12[k];
              SQMATRICES[2] = MatricesA21[k];
              SQMATRICES[3] = MatricesA22[k];
              SQMATRICES[0]->Reset();
              SQMATRICES[1]->Reset();
              SQMATRICES[2]->Reset();
              SQMATRICES[3]->Reset();

              N_SquareMatrices = 4;
              N_RectMatrices = 0;

              N_Rhs = 2;
              N_FESpaces = 1;
              rhs = RhsArray[k];
              RHSs[0] = rhs;
              RHSs[1] = rhs + N_Uarray[k];
              memset(rhs, 0, (2*N_Uarray[k]+N_Parray[k])*SizeOfDouble);
              ferhs[0] = USpaces[k];
              ferhs[1] = USpaces[k];
              last_sq = 3;

              if (DiscreteForm == DiscreteFormNLSDFEM)
              {
                N_RectMatrices = 2;
                MATRICES[0] = MatricesB1T[k];
                MATRICES[1] = MatricesB2T[k];
                MATRICES[0]->Reset();
                MATRICES[1]->Reset();
                N_FESpaces = 2;
              }
            }
            break;
        }                        // endswitch

        if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
          if ((DiscreteForm == DiscreteFormNLSmagorinsky)||
	      (DiscreteForm == DiscreteFormNLVMSProjection))
        {
          aux =  new TAuxParam2D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo,
            NSN_ParamFctVelo_GradVelo,
            NSN_FEValuesVelo_GradVelo,
            fesp, fefct,
            NSFctVelo_GradVelo,
            NSFEFctIndexVelo_GradVelo, NSFEMultiIndexVelo_GradVelo,
            NSN_ParamsVelo_GradVelo, NSBeginParamVelo_GradVelo);
        }
        else
        {
          aux =  new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo, NSN_ParamFctVelo,
            NSN_FEValuesVelo,
            fesp, fefct,
            NSFctVelo,
            NSFEFctIndexVelo, NSFEMultiIndexVelo,
            NSN_ParamsVelo, NSBeginParamVelo);
        }
        else                     // Newton method
        {
          aux =  new TAuxParam2D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo,
            NSN_ParamFctVelo_GradVelo,
            NSN_FEValuesVelo_GradVelo,
            fesp, fefct,
            NSFctVelo_GradVelo,
            NSFEFctIndexVelo_GradVelo, NSFEMultiIndexVelo_GradVelo,
            NSN_ParamsVelo_GradVelo, NSBeginParamVelo_GradVelo);
        }
        // assembling

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
              UpwindForNavierStokes(Coefficients[0], SQMATRICES[0], U1Array[k], U2Array[k]);
              cout << "UPWINDING DONE : level " << k << endl;
              break;

            case 3:
            case 4:
              // do upwinding with two matrices
              UpwindForNavierStokes(Coefficients[0], SQMATRICES[0], U1Array[k], U2Array[k]);
              UpwindForNavierStokes(Coefficients[0], SQMATRICES[last_sq], U1Array[k], U2Array[k]);
              cout << "UPWINDING DONE(2) : level " << k << endl;
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
          N_Rhs = 2;
          DiscreteForm0 = NULL;

          SQMATRICES[0] = MatricesA11[k];
          SQMATRICES[1] = MatricesA22[k];
          SQMATRICES[2] = MatricesA12[k];
          SQMATRICES[3] = MatricesA21[k];

          fesp[0] = USpaces[k];
          ferhs[0] = USpaces[k];
          ferhs[1] = USpaces[k];

          RHSs[0] = RhsArray[k];
          RHSs[1] = RhsArray[k]+N_Uarray[k];

          Assemble2DSlipBC(N_FESpaces, fesp,
            N_SquareMatrices, SQMATRICES,
            N_RectMatrices, MATRICES,
            N_Rhs, RHSs, ferhs,
            DiscreteForm0,
            BoundaryConditions,
            BoundValues,
            aux,
            U1Array[k],U2Array[k]);
        }
        delete aux;
	
	if (DiscreteForm == DiscreteFormNLVMSProjection)
	{
	    SQMATRICES[0] = MatricesA[k];
	    SQMATRICES[1] = MatricesL[k];
	    MATRICES[2] = Matrices_tilde_G11[k];
	    MATRICES[3] = Matrices_tilde_G22[k];
	    MATRICES[4] = Matrices_G11[k];
	    MATRICES[5] = Matrices_G22[k];

	    LumpMassMatrixToDiag(MatricesL[k]);
	    VMSProjectionUpdateMatrices(N_Uarray[k],USpaces[k]->GetActiveBound(),
					ProjectionSpaces[k]->GetN_DegreesOfFreedom(),
					SQMATRICES,MATRICES);
	}
      }                          // endfor k, assembling done

      if ((TDatabase::ParamDB->PRESSURE_SEPARATION==2)&&(j==1))
      {
        OutPut("apply pressure separation"<<endl);
        // save original pressure
        if ((j==1)||(1))
          memcpy(nosep_p,sol+2*N_U, N_P*SizeOfDouble);
        else
        {
          Daxpy(N_P, 1.0, sol+2*N_U, nosep_p);
          memcpy(sol+2*N_U, nosep_p, N_P*SizeOfDouble);
        }
        // assemble separated pressure entry on right hand side
        // first : interpolate discrete normal pressure to the
        //         pressure separation space
        Prolongate(pressure_space, pressure_separation_space,
          sol+2*N_U, separated_pressure_array,
          separated_pressure_aux_array);

        // second : assemble
        // the gradient of the separated pressure is needed for assembling
        // this has to be said to the assembling routine by an aux object
        fesp[0] = USpaces[mg_level];
        fesp[1] = pressure_separation_space;

        fefct[0] = separated_pressure_fe_funct;

        aux =  new TAuxParam2D(NSN_FESpacesPressSep, NSN_FctPressSep, NSN_ParamFctPressSep,
          NSN_FEValuesPressSep,
          fesp+1, fefct,
          NSFctPressSep,
          NSFEFctIndexPressSep, NSFEMultiIndexPressSep,
          NSN_ParamsPressSep, NSBeginParamPressSep);

        // assemble the right hand side
        N_FESpaces = 2;
        N_SquareMatrices = 0;
        N_RectMatrices = 0;
        N_Rhs = 2;
        RHSs[0] = rhs;
        RHSs[1] = rhs + N_U;
        memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);
        ferhs[0] = USpaces[mg_level];
        ferhs[1] = USpaces[mg_level];
        DiscreteForm = DiscreteFormPressSep;

        Assemble2D(N_FESpaces, fesp,
          N_SquareMatrices, SQMATRICES,
          N_RectMatrices, MATRICES,
          N_Rhs, RHSs, ferhs,
          DiscreteForm,
          BoundaryConditions,
          BoundValues,
          aux);

        // initialize solution array for separated pressure
        memset(sol+2*N_U,0, N_P*SizeOfDouble);
      }

      if ((TDatabase::ParamDB->PRESSURE_SEPARATION==5)&&(j==1))
      {
        OutPut("assemble pressure separation problem"<<endl);
        N_SquareMatrices = 1;
        N_RectMatrices = 0;
        N_Rhs = 1;
        SQMATRICES[0] = sqmatrixPressSep;
        SQMATRICES[0]->Reset();
        fesp[0] = pressure_separation_space;
        ferhs[0] = pressure_separation_space;
        N_FESpaces = 2;
        fesp[1] = USpaces[mg_level];
        fefct[0] = U1Array[mg_level];
        fefct[1] = U2Array[mg_level];
        aux =  new TAuxParam2D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo,
          NSN_ParamFctVelo_GradVelo,
          NSN_FEValuesVelo_GradVelo,
          fesp+1, fefct,
          NSFctVelo_GradVelo,
          NSFEFctIndexVelo_GradVelo, NSFEMultiIndexVelo_GradVelo,
          NSN_ParamsVelo_GradVelo, NSBeginParamVelo_GradVelo);

        RHSs[0] = rhsPressSep;
        memset(rhsPressSep, 0, N_P_sep*SizeOfDouble);
        // assemble
        Assemble2D(N_FESpaces, fesp,
          N_SquareMatrices, SQMATRICES,
          0, NULL,
          N_Rhs, RHSs, ferhs,
          DiscreteFormAuxProbPressSep,
          BoundaryConditionsPressureSeparation,
          BoundaryValuesPressureSeparation,
          aux);

        // solve linear system

        Auxprec = new TMultiGridScaIte(MatVect_Scalar, Defect_Scalar, NULL,
          0, N_P_sep, AuxMG, 0);
        Auxitmethod = new TFgmresIte(MatVect_Scalar, Defect_Scalar, Auxprec,
          0, N_P_sep, 1);
        Auxitmethod_sol = new double[N_P_sep];
        Auxitmethod_rhs = new double[N_P_sep];
        memcpy(Auxitmethod_sol, separated_pressure_array, N_P_sep*SizeOfDouble);
        memcpy(Auxitmethod_rhs, rhsPressSep, N_P_sep*SizeOfDouble);
        SQMATRICES[0] = sqmatrixPressSep;
        Auxitmethod->Iterate(sqmatrices,NULL,Auxitmethod_sol,Auxitmethod_rhs);
        memcpy(separated_pressure_array, Auxitmethod_sol, N_P_sep*SizeOfDouble);
        memcpy(rhsPressSep, Auxitmethod_rhs, N_P_sep*SizeOfDouble);
        delete Auxitmethod_sol;
        delete Auxitmethod_rhs;
        delete Auxprec;
        delete Auxitmethod;

        // store separated pressure in the original pressure space
        Prolongate(pressure_separation_space, pressure_space,
          separated_pressure_array, nosep_p,
          pressure_aux_array);

        // assemble rhs for NSE
        // the gradient of the separated pressure is needed for assembling
        // this has to be said to the assembling routine by an aux object
        fesp[0] = USpaces[mg_level];
        fesp[1] = pressure_separation_space;

        fefct[0] = separated_pressure_fe_funct;

        aux =  new TAuxParam2D(NSN_FESpacesPressSep, NSN_FctPressSep, NSN_ParamFctPressSep,
          NSN_FEValuesPressSep,
          fesp+1, fefct,
          NSFctPressSep,
          NSFEFctIndexPressSep, NSFEMultiIndexPressSep,
          NSN_ParamsPressSep, NSBeginParamPressSep);

        // assemble the right hand side
        N_FESpaces = 2;
        N_SquareMatrices = 0;
        N_RectMatrices = 0;
        N_Rhs = 2;
        RHSs[0] = rhs;
        RHSs[1] = rhs + N_U;
        memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);
        ferhs[0] = USpaces[mg_level];
        ferhs[1] = USpaces[mg_level];
        DiscreteForm = DiscreteFormPressSep;

        Assemble2D(N_FESpaces, fesp,
          N_SquareMatrices, SQMATRICES,
          N_RectMatrices, MATRICES,
          N_Rhs, RHSs, ferhs,
          DiscreteForm,
          BoundaryConditions,
          BoundValues,
          aux);

        // initialize solution array for separated pressure
        memset(sol+2*N_U,0, N_P*SizeOfDouble);
      }

      // RFB stabilization
      if (TDatabase::ParamDB->DISCTYPE==NSE_RFB)
      {
	  // assemble rhs
	  fesp[0] = USpaces[mg_level];
	  aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      
	  // assemble the right hand side
	  N_FESpaces = 1;
	  N_SquareMatrices = 0;
	  N_RectMatrices = 0;
	  N_Rhs = 2;
	  RHSs[0] = rhs;
	  RHSs[1] = rhs + N_U;
	  memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);
	  ferhs[0] = USpaces[mg_level];
	  ferhs[1] = USpaces[mg_level];
	  DiscreteForm = DiscreteFormNSRFBRhs;

	  Assemble2D(N_FESpaces, fesp,
		     N_SquareMatrices, SQMATRICES,
		     N_RectMatrices, MATRICES,
		     N_Rhs, RHSs, ferhs,
		     DiscreteForm,
		     BoundaryConditions,
		     BoundValues,
		     aux);
	  delete aux;
	  // modify rhs with RFB stabilization
	  //ApproximateRFBSolutionQuad_Q2_NSE2D(coll, U1Array[mg_level], U2Array[mg_level],
	  ApproximateRFBSolutionQuadNSE2D(coll, U1Array[mg_level], U2Array[mg_level],
	  				  Coefficients[0], rhs);
	  OutPut("RFB DONE"<<endl);
      }


      // end of assembling

      // reset MATRICES for solver
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          SQMATRICES[0] = MatricesA[mg_level];
          MATRICES[0] = MatricesB1[mg_level];
          MATRICES[1] = MatricesB2[mg_level];
          break;
        case 2:
          SQMATRICES[0] = MatricesA[mg_level];
          MATRICES[0] = MatricesB1[mg_level];
          MATRICES[1] = MatricesB2[mg_level];
          MATRICES[2] = MatricesB1T[mg_level];
          MATRICES[3] = MatricesB2T[mg_level];
          break;
        case 3:
          SQMATRICES[0] = MatricesA11[mg_level];
          SQMATRICES[1] = MatricesA12[mg_level];
          SQMATRICES[2] = MatricesA21[mg_level];
          SQMATRICES[3] = MatricesA22[mg_level];
          MATRICES[0] = MatricesB1[mg_level];
          MATRICES[1] = MatricesB2[mg_level];
          break;
        case 4:
          SQMATRICES[0] = MatricesA11[mg_level];
          SQMATRICES[1] = MatricesA12[mg_level];
          SQMATRICES[2] = MatricesA21[mg_level];
          SQMATRICES[3] = MatricesA22[mg_level];
          MATRICES[0] = MatricesB1[mg_level];
          MATRICES[1] = MatricesB2[mg_level];
          MATRICES[2] = MatricesB1T[mg_level];
          MATRICES[3] = MatricesB2T[mg_level];
          break;
        case 14:
          SQMATRICES[0] = MatricesA11[mg_level];
          SQMATRICES[1] = MatricesA12[mg_level];
          SQMATRICES[2] = MatricesA21[mg_level];
          SQMATRICES[3] = MatricesA22[mg_level];
          SQMATRICES[4] = MatricesC[mg_level];
          MATRICES[0] = MatricesB1[mg_level];
          MATRICES[1] = MatricesB2[mg_level];
          MATRICES[2] = MatricesB1T[mg_level];
          MATRICES[3] = MatricesB2T[mg_level];
          break;
      }

      // set rhs for Dirichlet nodes
      memcpy(sol+N_Active, rhs_high+N_Active, N_NonActive*SizeOfDouble);
      memcpy(sol+N_U+N_Active, rhs_high+N_U+N_Active, N_NonActive*SizeOfDouble);

      // compute defect
      memset(defect,0,N_Unknowns*SizeOfDouble);
      Defect(sqmatrices,matrices,sol,rhs_high,defect);
      if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
        IntoL20Vector2D(defect+2*N_U, N_P,pressure_space_code);
      residual =  Ddot(N_Unknowns,defect,defect);
      impuls_residual = Ddot(2*N_U,defect,defect);
      if (pressure_separation == 1)
        OutPut("p sep ");
      OutPut("nonlinear iteration step " << setw(3) << j);
      OutPut(setw(14) << impuls_residual);
      OutPut(setw(14) << residual-impuls_residual);
      OutPut(setw(14) << sqrt(residual) << endl);

      // check if convergence is too slow
      last_digit_ite = j%10;
      if (sqrt(residual)>convergence_speed*residuals[last_digit_ite])
        slow_conv=1;
      residuals[last_digit_ite]=sqrt(residual);

      // check if convergence criterion is fulfilled
      if ((sqrt(residual)<=limit)||(j==Max_It)||(slow_conv))
      {	      
        total_time2 = GetTime();
        if (pressure_separation == 1)
          OutPut("P SEP ");
        OutPut("ITE : " << setw(3) << j);
        OutPut(" (" << setw(3) << N_LinIter << " LINITE)");
        OutPut(" T/SOLVER : " << solver_time << "s");
        OutPut(" T/TOTAL : " << total_time2-total_time1 << "s");
        OutPut("("<<setprecision(2)<<solver_time/(total_time2-total_time1)<<")");
        OutPut(setprecision(6)<<" RES : " <<  sqrt(residual));
        if (slow_conv)
          OutPut(" SLOW !!!");
        OutPut(endl);
        if (TDatabase::ParamDB->PRESSURE_SEPARATION==0)
          break;
        if (((TDatabase::ParamDB->PRESSURE_SEPARATION==1)||
          (TDatabase::ParamDB->PRESSURE_SEPARATION==2))&&
          (pressure_separation == 1))
        {
          // compute pressure
          memcpy(sol+2*N_U, nosep_p, N_P*SizeOfDouble);
          //Daxpy(N_P,1.0,nosep_p,sol+2*N_U);
          break;
        }
        if ((TDatabase::ParamDB->PRESSURE_SEPARATION==3)||
          (TDatabase::ParamDB->PRESSURE_SEPARATION==4)||
          (TDatabase::ParamDB->PRESSURE_SEPARATION==5))
        {
          Daxpy(N_P,1.0,nosep_p,sol+2*N_U);
          break;
        }
        // save original pressure
        memcpy(nosep_p,sol+2*N_U, N_P*SizeOfDouble);

        // assemble separated pressure entry on right hand side
        // first : interpolate discrete normal pressure to the
        //         pressure separation space
        Prolongate(pressure_space, pressure_separation_space,
          sol+2*N_U, separated_pressure_array,
          separated_pressure_aux_array);

        // TESTING: USE INTERPOLATION
        // separated_pressure_fe_funct->Interpolate(ExactP);
        // second : assemble
        // the gradient of the separated pressure is needed for assembling
        // this has to be said to the assembling routine by an aux object
        fesp[0] = USpaces[mg_level];
        fesp[1] = pressure_separation_space;

        fefct[0] = separated_pressure_fe_funct;

        aux =  new TAuxParam2D(NSN_FESpacesPressSep, NSN_FctPressSep, NSN_ParamFctPressSep,
          NSN_FEValuesPressSep,
          fesp+1, fefct,
          NSFctPressSep,
          NSFEFctIndexPressSep, NSFEMultiIndexPressSep,
          NSN_ParamsPressSep, NSBeginParamPressSep);

        // assemble the right hand side
        N_FESpaces = 2;
        N_SquareMatrices = 0;
        N_RectMatrices = 0;
        N_Rhs = 2;
        RHSs[0] = rhs;
        RHSs[1] = rhs + N_U;
        memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);
        ferhs[0] = USpaces[mg_level];
        ferhs[1] = USpaces[mg_level];
        DiscreteForm = DiscreteFormPressSep;

        Assemble2D(N_FESpaces, fesp,
          N_SquareMatrices, SQMATRICES,
          N_RectMatrices, MATRICES,
          N_Rhs, RHSs, ferhs,
          DiscreteForm,
          BoundaryConditions,
          BoundValues,
          aux);

        // initialize solution array for separated pressure
        memset(sol+2*N_U,0, N_P*SizeOfDouble);
        // compute defect
        memset(defect,0,N_Unknowns*SizeOfDouble);
        Defect(sqmatrices,matrices,sol,rhs_high,defect);
        if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
          IntoL20Vector2D(defect+2*N_U, N_P,pressure_space_code);
        residual =  Ddot(N_Unknowns,defect,defect);
        impuls_residual = Ddot(2*N_U,defect,defect);
        // reset counter
        j = 0;
        pressure_separation = 1;
        OutPut("p sep ");
        OutPut("nonlinear iteration step " << setw(3) << j);
        OutPut(setw(14) << impuls_residual);
        OutPut(setw(14) << residual-impuls_residual);
        OutPut(setw(14) << sqrt(residual) << endl);
      }                          // end (of stopping criterion reached)

      switch(TDatabase::ParamDB->SOLVER_TYPE)
      {
        case AMG:
	    TDatabase::ParamDB->SC_VERBOSE=0;
          TDatabase::ParamDB->CC_VERBOSE=0;
          t1 = GetTime();
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
              Solver(sqmatrixA, matrixB1, matrixB2, rhs_high, sol);
              break;

            case 2:
          #ifdef __MORTAR__
              Solver(sqmatrixA, matrixB1T, matrixB2T,
                matrixB1, matrixB2, matrix_mortar, rhs_high, sol);
          #else
              Solver(sqmatrixA, matrixB1T, matrixB2T,
                matrixB1, matrixB2, rhs_high, sol);
          #endif
              break;

            case 3:
              Solver(sqmatrixA11, sqmatrixA12, sqmatrixA21,
                sqmatrixA22, matrixB1, matrixB2, rhs_high, sol);
              break;

            case 4:
              Solver(sqmatrixA11, sqmatrixA12, sqmatrixA21,
                sqmatrixA22, matrixB1T, matrixB2T,
                matrixB1, matrixB2, rhs_high, sol);
              break;
	    default:
		OutPut("AMG solver for NSTYPE " << TDatabase::ParamDB->NSTYPE 
		       << " not implemented !!!" << endl);
		exit(4711);

          }
          t2 = GetTime();
          solver_time += (t2-t1);
          break;
      case DIRECT:
        t1 = GetTime();
        switch(TDatabase::ParamDB->NSTYPE)
        {
	    case 1: 
		DirectSolver(sqmatrixA, matrixB1, matrixB2, rhs_high, sol);
            break;

	    case 2:
		DirectSolver(sqmatrixA, matrixB1T, matrixB2T,
			     matrixB1, matrixB2, rhs_high, sol);
		break;
	    case 4:
		DirectSolver(sqmatrixA11, sqmatrixA12, sqmatrixA21,
			     sqmatrixA22, matrixB1T, matrixB2T,
			     matrixB1, matrixB2, rhs_high, sol);
		break;
	    case 14:
		DirectSolver(sqmatrixA11, sqmatrixA12, sqmatrixA21,
			     sqmatrixA22, sqmatrixC, matrixB1T, matrixB2T,
			     matrixB1, matrixB2, rhs_high, sol);
		break;
	    default:
		OutPut("Direct solver for NSTYPE " << TDatabase::ParamDB->NSTYPE
		       << " not implemented !!!" << endl);
		exit(4711);
        }
        t2 = GetTime();
        solver_time += (t2-t1);
	break;

        case GMG:
          t1 = GetTime();
          switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
          {
            case 11:
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
              {
                memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
                memcpy(itmethod_rhs, rhs_high, N_Unknowns*SizeOfDouble);
              }
              break;
            case 16:
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
              {
                memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
                memcpy(itmethod_rhs, rhs_high, N_Unknowns*SizeOfDouble);
              }
              break;
          }
          N_LinIter+=itmethod->Iterate(sqmatrices,matrices,itmethod_sol,itmethod_rhs);
          switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
          {
            case 11:
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
              {
                memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
                memcpy(rhs_high, itmethod_rhs, N_Unknowns*SizeOfDouble);
              }
              break;
            case 16:
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
              {
                memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
                memcpy(rhs_high, itmethod_rhs, N_Unknowns*SizeOfDouble);
              }
              break;
          }
          t2 = GetTime();
          solver_time += (t2-t1);
          for(l=0;l<N_Unknowns;l++)
          {
            p2 = sol[l]-oldsol[l];
            sol[l] = oldsol[l] + omega * p2;
          }
          break;
      }
    }                            // endfor
 
    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
      IntoL20FEFunction(sol+2*N_U, N_P,PSpaces[mg_level],
        velocity_space_code, pressure_space_code);

    psi = new double[N_V];
    StreamFunction(velocity_space, sol, sol+N_U, streamfunction_space, psi);
    StreamFct = new TFEFunction2D(streamfunction_space, PsiString, PsiString, psi, N_V);
    // DivU(u1, u2, StreamFct);
    Output->AddFEFunction(StreamFct);
    /*
        for (l=0;l<N_U;l++)
           OutPut("u1(" << l+1 << ") = " << sol[l] << ";" << endl);
        for (l=0;l<N_U;l++)
         OutPut("u2(" << l+1 << ") = " << sol[l+N_U] << ";" << endl);
        for (l=0;l<N_P;l++)
	OutPut("p(" << l+1 << ") = " << sol[2*N_U+l] << ";" << endl);*/
     
    if(TDatabase::ParamDB->WRITE_GRAPE)
    {
      // OutPut("u " << Ddot(2*N_U,sol,sol)<< endl);
      //for (ii=0;ii<N_Unknowns;ii++)
      //  OutPut(ii << " " << sol[ii] << endl);
      // OutPut("p " << Ddot(N_P,sol+2*N_U,sol+2*N_U)<< endl);
      os.seekp(std::ios::beg);
      os << GrapeBaseName << i << ".dat" << ends;
      Output->WriteGrape(os.str().c_str());
    }

    if(TDatabase::ParamDB->WRITE_GNU)
    {
      os.seekp(std::ios::beg);
      os << GnuBaseName << i << ".gnu" << ends;
      Output->WriteGnuplot(os.str().c_str());
    }
    if(TDatabase::ParamDB->WRITE_VTK)
    {
      os.seekp(std::ios::beg);
      os << VtkBaseName << i << ".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
    }
    if(TDatabase::ParamDB->WRITE_GMV)
    {
      os.seekp(std::ios::beg);
      os << GmvBaseName << i << ".gmv" << ends;
      Output->WriteGMV(os.str().c_str());
    }
    if(TDatabase::ParamDB->WRITE_MATLAB)
    {
      os.seekp(std::ios::beg);
      os << MatlabBaseName << i << ".m" << ends;
      Output->WriteMatlab(os.str().c_str());
    }

   // measure errors to known solution
    if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
      fesp[0] = USpaces[mg_level];
      fefct[0] = U1Array[mg_level];
      fefct[1] = U2Array[mg_level];
      aux =  new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo,
        NSN_ParamFctVelo,
        NSN_FEValuesVelo,
        fesp, fefct,
        NSFctVelo,
        NSFEFctIndexVelo, NSFEMultiIndexVelo,
        NSN_ParamsVelo, NSBeginParamVelo);

      // errors in first velocity component
      U1Array[mg_level]->GetErrors(ExactU1, 3, NSAllDerivatives, 2,
        L2H1Errors,
        NULL, aux, 1, USpaces+mg_level, errors);
      l2u1[i] = errors[0];
      h1u1[i] = errors[1];

      // errors in second velocity component
      U2Array[mg_level]->GetErrors(ExactU2, 3, NSAllDerivatives, 2,
        L2H1Errors,
        NULL, aux, 1, USpaces+mg_level, errors);
      l2u2[i] = errors[0];
      h1u2[i] = errors[1];

      // errors in pressure
      PArray[mg_level]->GetErrors(ExactP, 3, NSAllDerivatives, 2,
        L2H1Errors,
        NULL, aux, 1, PSpaces+mg_level, errors);
      l2p[i] = errors[0];
      h1p[i] = errors[1];


      // output of errors
      // coarsest level
      if(i<=FirstSolve)
      {
        p1 = sqrt(l2u1[i]*l2u1[i]+l2u2[i]*l2u2[i]);
        OutPut("L2(u): " <<  p1 << endl);
        p2 = sqrt(h1u1[i]*h1u1[i]+h1u2[i]*h1u2[i]);
        OutPut("H1-semi(u): " <<  p2 << endl);
        OutPut("L2(p): " <<  l2p[i] << endl);
        OutPut("H1-semi(p): " <<  h1p[i] << endl);
      }
      // not coarsest level
      else
      {
        errors[0] = sqrt(l2u1[i]*l2u1[i]+l2u2[i]*l2u2[i]) ;
        errors[1] = sqrt(l2u1[i-1]*l2u1[i-1]+l2u2[i-1]*l2u2[i-1]);
        OutPut("L2(u): " <<  errors[0] << " order ");
        OutPut(log(errors[1]/errors[0])/ln2 << endl);

        errors[0] = sqrt(h1u1[i]*h1u1[i]+h1u2[i]*h1u2[i]) ;
        errors[1] = sqrt(h1u1[i-1]*h1u1[i-1]+h1u2[i-1]*h1u2[i-1]);
        OutPut("H1-semi(u): " << errors[0] << " order ");
        OutPut(log(errors[1]/errors[0])/ln2 << endl);

        errors[0] = l2p[i];
        errors[1] = l2p[i-1];
        OutPut("L2(p): " <<  errors[0] << " order ");
        OutPut(log(errors[1]/errors[0])/ln2 << endl);

        errors[0] = h1p[i];
        errors[1] = h1p[i-1];
        OutPut("H1-semi(p): " << errors[0] << " order ");
        OutPut(log(errors[1]/errors[0])/ln2 << endl);

      }
      UArray[mg_level]->GetDeformationTensorErrors
          (ExactNull, ExactNull,
	   3, NSAllDerivatives,
	   2, DivergenceError,
	   NULL, aux, 1, USpaces+mg_level, errors);
      
      OutPut( "divergence error (L1/L2) : " << errors[0]*errors[0] << " " << errors[1] << endl);
      delete aux;
    }                            // endif MEASURE_ERRORS

    #ifdef __BENCH__
    // compute characteristic values (deltaP, Cd, Cl)
    GetCdCl(U1Array[mg_level], U2Array[mg_level], PArray[mg_level],
      U1Array[mg_level], U2Array[mg_level], Cd, Cl);

    PArray[mg_level]->FindGradient(0.15, 0.2, dP1);
    PArray[mg_level]->FindGradient(0.25, 0.2, dP2);

    OutPut( "C_drag = " << setprecision(16) <<Cd );
    OutPut( " C_lift = " << setprecision(16) << Cl);
    OutPut( " deltaP = " << setprecision(16) << dP1[0] - dP2[0] << endl);
    OutPut( setprecision(7) << endl);
    #endif

    #ifdef __CHANNELSTEP__
    GetReattachmentPoint(U1Array[mg_level], reatt_pt);
    OutPut( "reattachment: " << reatt_pt<< endl);
    #endif

    #ifdef __CHANNEL30__
    GetReattachmentPoint(U1Array[mg_level], reatt_point);
    OutPut( "reattachment: lower " << reatt_point[0] << " left upper " <<
      reatt_point[1] << " right upper " <<  reatt_point[2] << endl);
    #endif

    // in case of convergence, set solution to zero
    // to prevent interpolation of useless data
    if (((slow_conv) || sqrt(residual)>limit )
	&& (!TDatabase::ParamDB->STOKES_PROBLEM))
    {
	//memset(sol, 0, N_Unknowns*SizeOfDouble);
      OutPut("solution on level " << i <<
        " set to zero since method did not converge" << endl);
    }

      fesp[0] = USpaces[mg_level];
      fefct[0] = U1Array[mg_level];
      fefct[1] = U2Array[mg_level];
      aux =  new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo,
        NSN_ParamFctVelo,
        NSN_FEValuesVelo,
        fesp, fefct,
        NSFctVelo,
        NSFEFctIndexVelo, NSFEMultiIndexVelo,
        NSN_ParamsVelo, NSBeginParamVelo);
    U1Array[mg_level]->GetErrors(ExactNull, 3, NSAllDerivatives, 2,
              L2H1Errors,
              NULL, aux, 1, USpaces+mg_level, errors_mg);
      delete aux;
  
    // remove data which will not be used later
    delete oldsol;
    delete defect;
    delete psi;
  
    if ((mg_type==1)||(TDatabase::ParamDB->SOLVER_TYPE == AMG)
	||(TDatabase::ParamDB->SOLVER_TYPE == DIRECT))
    {
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          delete matrixB1;
          delete matrixB2;
          delete sqmatrixA;
	  if (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION)
	  {
	      delete sqmatrixL;
	      delete  matrix_tilde_G11;
	      delete  matrix_tilde_G22;
	      delete  matrix_G11;
	      delete  matrix_G22;
	      delete  sqstructureL;
	      delete  structure_tilde_G;
	      delete  structure_G;
	      delete  projection_space;
	  }
          break;
        case 2:
          delete matrixB1;
          delete matrixB2;
          delete matrixB1T;
          delete matrixB2T;
          delete sqmatrixA;
          break;
        case 3:
          delete matrixB1;
          delete matrixB2;
          delete sqmatrixA11;
          delete sqmatrixA12;
          delete sqmatrixA21;
          delete sqmatrixA22;
          break;
        case 4:
          delete matrixB1;
          delete matrixB2;
          delete matrixB1T;
          delete matrixB2T;
          delete sqmatrixA11;
          delete sqmatrixA12;
          delete sqmatrixA21;
          delete sqmatrixA22;
          break;
      }

      delete structureB;
      delete structureBT;
      delete sqstructureA;
      delete rhs_high;
    }                            // end if (mg_type==1)

    old_sol = sol;
    old_u = u;
    old_p = p;
    old_u_space = velocity_space;
    old_p_space = pressure_space;

    OutPut("memory after: " << setw(10) << GetMemory() << endl);
  }                              // endfor i

  OutPut("used time: " << GetTime() << endl);
  OutPut("used bytes: " << GetMemory() << endl);
  CloseFiles();

  return 0;
}
