// =======================================================================
//
// Purpose:     main program
//
// Author:      Volker John
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
#include <Assemble2D.h>
#include <Output2D.h>
#include <DiscreteForm2D.h>
#include <LinAlg.h>
#include <TNSE2D_ParamRout.h>
#include <Collection.h>
#include <VMS.h>
#include <DirectSolver.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <malloc.h>

#include <Upwind.h>
#include <NSE_MultiGrid.h>
#include <NSE_MGLevel1.h>
#include <NSE_MGLevel2.h>
#include <NSE_MGLevel3.h>
#include <NSE_MGLevel4.h>
#include <Convolution.h>
#include <ItMethod.h>
#include <FixedPointIte.h>
#include <MultiGridIte.h>
#include <FgmresIte.h>

#include <MultiGrid2D.h>
#include <MGLevel2D.h>

#include <MainUtilities.h>
#include <TimeUtilities.h>

#define AMG 0
#define GMG 1

// =======================================================================
// include current example
// =======================================================================
//#include "../Examples/TNSE_2D/Bsp1.h"
// #include "../Examples/TNSE_2D/Bsp2.h"
// #include "../Examples/TNSE_2D/Bsp3.h"
// #include "../Examples/TNSE_2D/Bsp4.h"
// #include "../Examples/TNSE_2D/Bsp5.h"
// #include "../Examples/TNSE_2D/Bsp6.h"
// #include "../Examples/TNSE_2D/Bsp7.h"
// #include "../Examples/TNSE_2D/Bsp8.h"
//#include "../Examples/TNSE_2D/Bsp9.h"
// #include "../Examples/TNSE_2D/Bsp10.h"
// #include "../Examples/TNSE_2D/Benchmark2.new.h"
// #include "../Examples/TNSE_2D/Benchmark2.new_neum.h"
// #include "../Examples/TNSE_2D/Benchmark3.new.h" 
#include "../Examples/TNSE_2D/Benchmark3_Neum.h"
// #include "../Examples/TNSE_2D/Benchmark4.h"
// #include "../Examples/TNSE_2D/Channel.h"
// #include "../Examples/TNSE_2D/Step.h"
//#include "../Examples/TNSE_2D/Driven.h"
// #include "../Examples/TNSE_2D/DrivenTest.h"
//#include "../Examples/TNSE_2D/SinCosExp.h"
// #include "../Examples/TNSE_2D/PoiseuilleSlip.h"
// #include "../Examples/TNSE_2D/ChannelSlip.h"
// #include "../Examples/TNSE_2D/ChannelStep.h"
//#include "../Examples/TNSE_2D/ChannelStepSlip.h"
// #include "../Examples/TNSE_2D/ChannelStep_Bosch.h"
//#include "../Examples/TNSE_2D/SinCosExpT.h"
// #include "../Examples/TNSE_2D/SinCosExpTSma.h"
// #include "../Examples/TNSE_2D/SinCosExpTSma.h"
// #include "../Examples/TNSE_2D/PoiseuilleSma.h"
// #include "../Examples/TNSE_2D/SinCosExpTSmaNabla.h"
// #include "../Examples/TNSE_2D/SinCosExpT_Taylor.h"
// #include "../Examples/TNSE_2D/SinCosExpT_Smago.h"
// #include "../Examples/TNSE_2D/ChannelSlipPara.h"
// #include "../Examples/TNSE_2D/MixingLayer.h"
//#include "../Examples/TNSE_2D/MixingLayerSlip.h"
// #include "../Examples/TNSE_2D/MixingLayerSlipIni01.h"
// #include "../Examples/TNSE_2D/MixingLayerSlipSmallSquare.h"
//#include "Examples/TNSE_2D/UnitSquareHoleDrag.h"
// #include "../Examples/TNSE_2D/RotatingFlow.h"
//#include "../Examples/TNSE_2D/SSMUM_00.h"

//#include "../Examples/TNSE_2D/FF1.h"
//#include "BOSCH/data/bosch_0433175329_00.h"
//#include "BOSCH/data/bosch_0433175329_01.h"
//#include "Examples/NSE_2D/SFBCavity.h"

// ======================================================================
// utilities for main program
// ======================================================================

int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *coll, *mortarcoll = NULL;
  TBaseCell *cell;
  TFESpace2D *velocity_space, *pressure_space, *streamfunction_space, *convolution_space;
  TFESpace2D *vorticity_space, *projection_space;
  TFESpace2D **USpaces, **PSpaces, **PsiSpaces, **duConvSpaces, **uConvSpaces;
  TFESpace2D **VorticitySpaces, **ProjectionSpaces;
  TOutput2D *Output;

  double *B, *rhs, *sol, *oldsol, tol, tolerance, *psi, *defect, *startsol, *frac_step_sol;
  double *app, *oldrhs, *itmethod_sol, *itmethod_rhs, *vorticity, *div, *conv_vort;
  double *solGL00AuxProblem, *rhsGL00AuxProblem, *LESModelRhs, *sol_timestep_m1;
  double umfpack_flag=-2;
  int i,j,k,l,m,n, N_, Len, low;
  int N_Rows, N_Columns, N_U, N_P, N_Unknowns, N_V, N_Vort, img=1;
  double *l2u1, *l2u2, *h1u1, *h1u2;
  double *l2p, *h1p, *sd, *l_inf, vect[3], exactvect[3];
  int which, *permutation, N_GRAPE_images=0,  N_GNU_images=0;
  double DiffL2, DiffH1, t;
  char *PRM, *GEO;
  int LEVELS, BASELEVEL, ORDER, order;
  int ret, pde;
  double negPower;
  double x,y,max,min,sum;
  double RE_NR;
  double tau1, tau2;
  double errors[7], p1, p2, *save_sol[1];
  double t1, t2, res, res2, oldres, solver_time, solver_time_curr, residual, oldresidual;
  double impuls_residual,limit,linredfac, total_time, t3, t4;
  int N_LinIter, N_LinIterCurr, N_LinIterCurrIte, N_SubSteps, N_Active, n_aux;
  double gamma, tau, oldtau;
  int *RowPtr, save_N_Unknowns[1];
  double press_factor1, press_factor2, val0;

  std::ostringstream os;
  char *PsBaseName, *GrapeBaseName, *GnuBaseName, *ReadGrapeBaseName;
  char *GmvBaseName, *VtkBaseName, *SaveDataFileName, *ReadDataFileName;

  double *val, cd, cl;
  TFEFunction2D *u1, *u2, *p, *fefct[7], *StreamFct, *Vorticity, *Divergence, *Conv_Vort;
  TFEFunction2D *du1Conv, *du2Conv, *du3Conv;
  TFEFunction2D *u1Conv, *u2Conv;
  TFEFunction2D **U1Array, **U2Array, **AuxFEFunctArray;
  TFEFunction2D **PArray;
  TFEVectFunct2D *u, **UArray, *uconf, *duConv, **duConvArray;
  TFEVectFunct2D *uConv, **uConvArray, **AuxFEVectFunctArray;
  TFEFunction2D **du1ConvArray, **du2ConvArray, **du3ConvArray;
  TFEFunction2D **u1ConvArray, **u2ConvArray;
  TFEVectFunct2D *GL00AuxProblemSol, **GL00AuxProblemSolArray;
  TFEFunction2D *GL00AuxProblemSol11, *GL00AuxProblemSol12;
  TFEFunction2D *GL00AuxProblemSol22;
  TFEFunction2D **GL00AuxProblemSol11Array, **GL00AuxProblemSol12Array;
  TFEFunction2D **GL00AuxProblemSol22Array;
  TFESpace2D *fesp[4], *ferhs[3];
  double delta, end_time, l_infty_l_2 = 0, l_infty_l_2_time=-4711.0;
  double olderror = 0, l_2_l_2Du=0, l_2_l_2u=0 , olderror_l_2_l_2u=0;
  double l_2_h_1u=0, olderror_l_2_h_1u=0;

  TAuxParam2D *aux;

  TSquareStructure2D *sqstructureA, *sqstructureC;
  TStructure2D *structureB, *structureBT;
  TSquareMatrix2D *sqmatrixA, *SQMATRICES[8];
  TSquareMatrix2D *sqmatrixA11, *sqmatrixA12;
  TSquareMatrix2D *sqmatrixA21, *sqmatrixA22;
  TSquareMatrix2D *sqmatrixM;
  TSquareMatrix2D *sqmatrixM11, *sqmatrixM12;
  TSquareMatrix2D *sqmatrixM21, *sqmatrixM22;
  TSquareMatrix2D **MatricesM11, **MatricesM12;
  TSquareMatrix2D **MatricesM21, **MatricesM22;
  TSquareMatrix2D **MatricesA, **MatricesM;
  TSquareMatrix2D **MatricesA11, **MatricesA12;
  TSquareMatrix2D **MatricesA21, **MatricesA22;
  TSquareMatrix2D *sqmatrixGL00AuxProblem,**MatricesGL00AuxProblem;
  TMatrix2D *matrixB1, *matrixB2, *MATRICES[10];
  TMatrix2D *matrixB1T, *matrixB2T;
  TMatrix2D **MatricesB1, **MatricesB2, **MatricesB1T, **MatricesB2T;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TMatrix **matrices = (TMatrix **)MATRICES;
  MatVecProc *MatVect;
  DefectProc *Defect;

  TSquareStructure2D *sqstructureL;
  TStructure2D *structure_tilde_G, *structure_G;
  TSquareMatrix2D *sqmatrixL, **MatricesL;
  TMatrix2D *matrix_tilde_G11, *matrix_tilde_G22;
  TMatrix2D *matrix_G11, *matrix_G22;
  TMatrix2D **Matrices_tilde_G11, **Matrices_tilde_G22;
  TMatrix2D **Matrices_G11, **Matrices_G22;
  int N_L;

  double **RhsArray;

  TNSE_MGLevel *MGLevel;
  TNSE_MultiGrid *MG;
  TMGLevel2D *MGLevelGL00AuxProblem;
  TMultiGrid2D *MGGL00AuxProblem;

  double *RHSs[4];
  int *N_Uarray, *N_Parray;

  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormColetti;
  TDiscreteForm2D *DiscreteFormGL00Convolution;
  TDiscreteForm2D *DiscreteFormGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormUpwind;
  TDiscreteForm2D *DiscreteFormSmagorinsky;
  TDiscreteForm2D *DiscreteFormVMSProjection;

  TDiscreteForm2D *DiscreteFormNLGalerkin;
  TDiscreteForm2D *DiscreteFormNLUpwind;
  TDiscreteForm2D *DiscreteFormNLColetti;
  TDiscreteForm2D *DiscreteFormNLGL00Convolution;
  TDiscreteForm2D *DiscreteFormNLGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormNLSmagorinsky;
  TDiscreteForm2D *DiscreteFormNLVMSProjection;

  TDiscreteForm2D *DiscreteFormRHS;
  TDiscreteForm2D *DiscreteFormRHSColetti;
  TDiscreteForm2D *DiscreteFormRHSSmagorinskyExpl;
  TDiscreteForm2D *DiscreteFormMatrixGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormGL00AuxProblemRHS;
  TDiscreteForm2D *DiscreteFormRHSLESModel;
  TDiscreteForm2D *DiscreteFormRHSAuxProblemU;
  TDiscreteForm2D *DiscreteFormMatrixAuxProblemU;

  TDiscreteForm2D *DiscreteForm;

  int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;

  BoundCondFunct2D *BoundaryConditions[2], *BoundaryConditionsAuxProblem[3];
  BoundValueFunct2D *BoundValues[2], *BoundValuesAuxProblem[3];
  double average;

  TItMethod *itmethod, *prec;
  int Max_It, FirstSolve;
  double omega, alpha[3], divergence;
  int N_Paramters=10, methods, time_discs, **downwind;
  double Parameters[10], hmin, hmax;

  int N_UConv, level_down, ii, fully_implicit = 0;
  double *u_conv, *auxConv, *du_tensor, *u_uConv, reatt_pt, vort_zero, vort_zero_conv;
  int mg_level,mg_type,CurrentDiscType, last_sq, step_length;
  int velocity_space_code, pressure_space_code;
  int very_first_time=0, zerostart, comp_vort, mixing_layer=0;

  // strings
  char ReadinDat[] = "readin.dat";
  char NameString[] = "name";
  char UString[] = "u";
  char PString[] = "p";
  char PsiString[] = "psi";
  char UConvString[] = "u_conv";
  char UConfString[] = "uconf";
  char AuxProbString[] = "AuxProblem";
  char VorticityString[] = "vorticity";
  char ConvVortString[] = "conv_vort";
  char DivergenceString[] = "divergence";

  #ifdef __BENCH__
  double Cd, Cl, dP1[3], dP2[3], *former_sol;
  TFEFunction2D *U1old, *U2old;
  #endif
  #ifdef __MIXINGLAYERSLIP__
  mixing_layer = 1;
  #endif

  os << " ";

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
  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();
  //======================================================================
  // copy read parameters into local variables
  //======================================================================
  if( (TDatabase::ParamDB->DISCTYPE==2) )
  {
    OutPut("SDFEM does not work!" << endl);
    Error("SDFEM does not work!" << endl);
    exit(4711);
  }
  if( (TDatabase::ParamDB->DISCTYPE==5) )
  {
    OutPut("DISCTYPE 5 NOT IMPLEMENTED!" << endl);
    Error("DISCTYPE 5 NOT IMPLEMENTED!" << endl);
    exit(4711);
  }

  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;

  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  GrapeBaseName = TDatabase::ParamDB->GRAPEBASENAME;
  GnuBaseName = TDatabase::ParamDB->GNUBASENAME;
  ReadGrapeBaseName = TDatabase::ParamDB->READGRAPEBASENAME;
  GmvBaseName = TDatabase::ParamDB->GMVBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  SaveDataFileName = TDatabase::ParamDB->SAVE_DATA_FILENAME;
  ReadDataFileName = TDatabase::ParamDB->READ_DATA_FILENAME;

  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SADDLE;
  if (mg_type)
    mg_level = 1;
  else
    mg_level = 0;

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
  duConvArray = new TFEVectFunct2D*[LEVELS+1];
  du1ConvArray = new TFEFunction2D*[LEVELS+1];
  du2ConvArray = new TFEFunction2D*[LEVELS+1];
  du3ConvArray = new TFEFunction2D*[LEVELS+1];
  uConvArray = new TFEVectFunct2D*[LEVELS+1];
  u1ConvArray = new TFEFunction2D*[LEVELS+1];
  u2ConvArray = new TFEFunction2D*[LEVELS+1];
  GL00AuxProblemSolArray = new TFEVectFunct2D*[LEVELS+1];
  GL00AuxProblemSol11Array = new TFEFunction2D*[LEVELS+1];
  GL00AuxProblemSol12Array = new TFEFunction2D*[LEVELS+1];
  GL00AuxProblemSol22Array = new TFEFunction2D*[LEVELS+1];

  RhsArray = new double* [LEVELS+1];
  N_Uarray = new int[LEVELS+1];
  N_Parray = new int[LEVELS+1];

  USpaces = new TFESpace2D*[LEVELS+1];
  PSpaces = new TFESpace2D*[LEVELS+1];
  PsiSpaces = new TFESpace2D*[LEVELS+1];
  VorticitySpaces = new TFESpace2D*[LEVELS+1];
  ProjectionSpaces = new TFESpace2D*[LEVELS+1];
  duConvSpaces = new TFESpace2D*[LEVELS+1];
  uConvSpaces = new TFESpace2D*[LEVELS+1];

  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      MatricesA = new TSquareMatrix2D*[LEVELS+1];
      MatricesM = new TSquareMatrix2D*[LEVELS+1];

      MatricesB1 = new TMatrix2D*[LEVELS+1];
      MatricesB2 = new TMatrix2D*[LEVELS+1];
      MatVect = MatVect_NSE1;
      Defect = Defect_NSE1;
      break;

    case 2:
      MatricesA = new TSquareMatrix2D*[LEVELS+1];
      MatricesM = new TSquareMatrix2D*[LEVELS+1];

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
      MatricesM11 = new TSquareMatrix2D*[LEVELS+1];
      MatricesM12 = new TSquareMatrix2D*[LEVELS+1];
      MatricesM21 = new TSquareMatrix2D*[LEVELS+1];
      MatricesM22 = new TSquareMatrix2D*[LEVELS+1];

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
      MatricesM11 = new TSquareMatrix2D*[LEVELS+1];
      MatricesM12 = new TSquareMatrix2D*[LEVELS+1];
      MatricesM21 = new TSquareMatrix2D*[LEVELS+1];
      MatricesM22 = new TSquareMatrix2D*[LEVELS+1];

      MatricesB1 = new TMatrix2D*[LEVELS+1];
      MatricesB2 = new TMatrix2D*[LEVELS+1];
      MatricesB1T = new TMatrix2D*[LEVELS+1];
      MatricesB2T = new TMatrix2D*[LEVELS+1];
      MatVect = MatVect_NSE4;
      Defect = Defect_NSE4;
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
  }

  // array of matrices for auxiliary problem in Galdi/Layton model
  if ((TDatabase::ParamDB->DISCTYPE==GL00_AUX_PROBLEM)
    ||(TDatabase::ParamDB->CONVOLUTE_SOLUTION)
    ||(TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE==4))
  {
    MatricesGL00AuxProblem = new TSquareMatrix2D*[LEVELS+1];
  }

  downwind = new int*[LEVELS+1];

  //======================================================================
  // creating discrete forms
  //======================================================================

  InitializeDiscreteForms(DiscreteFormGalerkin,DiscreteFormUpwind,
    DiscreteFormSmagorinsky,DiscreteFormColetti,
    DiscreteFormGL00Convolution,DiscreteFormGL00AuxProblem,
    DiscreteFormVMSProjection,
    DiscreteFormNLGalerkin,
    DiscreteFormNLUpwind, DiscreteFormNLSmagorinsky,
    DiscreteFormNLColetti,DiscreteFormNLGL00Convolution,
    DiscreteFormNLGL00AuxProblem,
    DiscreteFormNLVMSProjection,
    DiscreteFormRHS,
    DiscreteFormRHSColetti,
    DiscreteFormRHSLESModel,
    DiscreteFormMatrixGL00AuxProblem,
    DiscreteFormGL00AuxProblemRHS,
    DiscreteFormRHSSmagorinskyExpl,
    DiscreteFormMatrixAuxProblemU,
    DiscreteFormRHSAuxProblemU,
    LinCoeffs, TDatabase::ParamDB->NSTYPE);

  //======================================================================
  // read boundary parameterization and initialize coarse grid
  //======================================================================
  Domain->Init(PRM, GEO);

  Domain->ConvertQuadToTri(TDatabase::ParamDB->CONVERT_QUAD_TO_TRI);

  BoundaryConditions[0] = BoundCondition;
  BoundaryConditions[1] = BoundCondition;

  BoundValues[0] = U1BoundValue;
  BoundValues[1] = U2BoundValue;

  BoundaryConditionsAuxProblem[0] = BoundConditionAuxProblem;
  BoundaryConditionsAuxProblem[1] = BoundConditionAuxProblem;
  BoundaryConditionsAuxProblem[2] = BoundConditionAuxProblem;

  BoundValuesAuxProblem[0] = BoundValueAuxProblem;
  BoundValuesAuxProblem[1] = BoundValueAuxProblem;
  BoundValuesAuxProblem[2] = BoundValueAuxProblem;

  for(i=0;i<TDatabase::ParamDB->SC_COARSEST_LEVEL_SADDLE;i++)
    Domain->RegRefineAll();

  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
  omega = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SADDLE;
  divergence = TDatabase::ParamDB->SC_DIV_FACTOR;

  Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
  Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
  Parameters[2] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
  Parameters[9] = 0;

  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;
  SetTimeDiscParameters();

  if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    i=1;
    MG = new TNSE_MultiGrid(i, N_Paramters, Parameters);
  }

  // multigrid for Galdi/Layton model with auxiliary problem
  if (TDatabase::ParamDB->DISCTYPE==GL00_AUX_PROBLEM)
  {
    MGGL00AuxProblem = new TMultiGrid2D(i, N_Paramters, Parameters);
  }
  t3 = GetTime();
  total_time = t3 - total_time;
  mg_level = LEVELS+mg_level;
  //======================================================================
  // loop over all levels
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

    if(i && (i<LEVELS)) Domain->RegRefineAll();

    coll=Domain->GetCollection(It_Finest, 0);
    cout << endl << endl;

    if(TDatabase::ParamDB->WRITE_PS)
    {
      // write grid into an Postscript file
      os.seekp(std::ios::beg);
      os << PsBaseName << i << ".ps" << ends;
      Domain->PS(os.str().c_str(),It_Finest,0);
    }

    // get spaces for low order disc on finest geo grid
    if ((mg_type==1)&&(i<mg_level-1))
    {
      velocity_space = new TFESpace2D(coll,NameString, UString, BoundCondition,
        Non_USpace,1, mortarcoll);
      pressure_space = new TFESpace2D(coll,NameString, PString, BoundCondition,
        DiscP_PSpace,0, mortarcoll);
      convolution_space = new TFESpace2D(coll,NameString, UString, BoundConditionAuxProblem,
        Non_USpace,1, mortarcoll);
      velocity_space_code = -1;
      pressure_space_code = 0;
      order = -1;
    }
    // get spaces of high order disc on finest geo grid
    else
    {
      GetVelocityAndPressureSpace(coll,BoundCondition,
        mortarcoll, velocity_space,
        pressure_space, &pressure_space_code,
        TDatabase::ParamDB->VELOCITY_SPACE,
        TDatabase::ParamDB->PRESSURE_SPACE);

      velocity_space_code =   TDatabase::ParamDB->VELOCITY_SPACE;
      TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;

      GetVelocityAndPressureSpace(coll,BoundConditionAuxProblem,
        mortarcoll, convolution_space,
        pressure_space, &pressure_space_code,
        TDatabase::ParamDB->VELOCITY_SPACE,
        TDatabase::ParamDB->PRESSURE_SPACE);

      OutPut("convolution space is velocity space" << endl);
    }

    streamfunction_space = new TFESpace2D(coll,NameString, PsiString, BoundCondition,
      1, NULL);

    vorticity_space = new TFESpace2D(coll,NameString, UString, BoundCondition,
      ContP_USpace,1, mortarcoll);

    // OutPut("convolution space is conforming first order" << endl);
    //convolution_space = new TFESpace2D(coll,NameString, UString, BoundCondition,
    //                                  ContP_USpace,1, mortarcoll);

    // build fespace hierarchy
    // set values and pointers for fe space
    USpaces[i] = velocity_space;
    PSpaces[i] = pressure_space;
    N_U = velocity_space->GetN_DegreesOfFreedom();
    N_P = pressure_space->GetN_DegreesOfFreedom();
    N_Uarray[i] = N_U;
    N_Parray[i] = N_P;
    N_Active = velocity_space->GetActiveBound();
    if (i<LEVELS)
    {
      PsiSpaces[i] = streamfunction_space;
      N_V = streamfunction_space->GetN_DegreesOfFreedom();
      VorticitySpaces[i] = vorticity_space;
      N_Vort = vorticity_space->GetN_DegreesOfFreedom();
    }

    // build matrices
    structureB = new TStructure2D(pressure_space, velocity_space);
    structureBT = new TStructure2D(velocity_space, pressure_space);
    sqstructureA = new TSquareStructure2D(velocity_space);
    sqstructureA->Sort();

    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        matrixB1 = new TMatrix2D(structureB);
        matrixB2 = new TMatrix2D(structureB);

        MatricesB1[i] = matrixB1;
        MatricesB2[i] = matrixB2;

        sqmatrixA = new TSquareMatrix2D(sqstructureA);
        MatricesA[i] = sqmatrixA;

        sqmatrixM = new TSquareMatrix2D(sqstructureA);
        MatricesM[i] = sqmatrixM;
        break;

      case 2:
        matrixB1 = new TMatrix2D(structureB);
        matrixB2 = new TMatrix2D(structureB);
        matrixB1T = new TMatrix2D(structureBT);
        matrixB2T = new TMatrix2D(structureBT);

        MatricesB1[i] = matrixB1;
        MatricesB2[i] = matrixB2;
        MatricesB1T[i] = matrixB1T;
        MatricesB2T[i] = matrixB2T;

        sqmatrixA = new TSquareMatrix2D(sqstructureA);
        MatricesA[i] = sqmatrixA;

        sqmatrixM = new TSquareMatrix2D(sqstructureA);
        MatricesM[i] = sqmatrixM;
        break;

      case 3:
        matrixB1 = new TMatrix2D(structureB);
        matrixB2 = new TMatrix2D(structureB);

        MatricesB1[i] = matrixB1;
        MatricesB2[i] = matrixB2;

        sqmatrixA11 = new TSquareMatrix2D(sqstructureA);
        sqmatrixA12 = new TSquareMatrix2D(sqstructureA);
        sqmatrixA21 = new TSquareMatrix2D(sqstructureA);
        sqmatrixA22 = new TSquareMatrix2D(sqstructureA);

        MatricesA11[i] = sqmatrixA11;
        MatricesA12[i] = sqmatrixA12;
        MatricesA21[i] = sqmatrixA21;
        MatricesA22[i] = sqmatrixA22;

        sqmatrixM11 = new TSquareMatrix2D(sqstructureA);
        sqmatrixM12 = new TSquareMatrix2D(sqstructureA);
        sqmatrixM21 = new TSquareMatrix2D(sqstructureA);
        sqmatrixM22 = new TSquareMatrix2D(sqstructureA);

        MatricesM11[i] = sqmatrixM11;
        MatricesM12[i] = sqmatrixM12;
        MatricesM21[i] = sqmatrixM21;
        MatricesM22[i] = sqmatrixM22;
        break;

      case 4:
        matrixB1 = new TMatrix2D(structureB);
        matrixB2 = new TMatrix2D(structureB);
        matrixB1T = new TMatrix2D(structureBT);
	matrixB2T = new TMatrix2D(structureBT);

        MatricesB1[i] = matrixB1;
        MatricesB2[i] = matrixB2;
        MatricesB1T[i] = matrixB1T;
        MatricesB2T[i] = matrixB2T;

        sqmatrixA11 = new TSquareMatrix2D(sqstructureA);
        sqmatrixA12 = new TSquareMatrix2D(sqstructureA);
        sqmatrixA21 = new TSquareMatrix2D(sqstructureA);
        sqmatrixA22 = new TSquareMatrix2D(sqstructureA);

        MatricesA11[i] = sqmatrixA11;
        MatricesA12[i] = sqmatrixA12;
        MatricesA21[i] = sqmatrixA21;
        MatricesA22[i] = sqmatrixA22;

        sqmatrixM11 = new TSquareMatrix2D(sqstructureA);
        sqmatrixM12 = new TSquareMatrix2D(sqstructureA);
        sqmatrixM21 = new TSquareMatrix2D(sqstructureA);
        sqmatrixM22 = new TSquareMatrix2D(sqstructureA);

        MatricesM11[i] = sqmatrixM11;
        MatricesM12[i] = sqmatrixM12;
        MatricesM21[i] = sqmatrixM21;
        MatricesM22[i] = sqmatrixM22;
        break;
    }

    // allocate matrix for auxiliary problem in Galdi/Layton model
    if (TDatabase::ParamDB->DISCTYPE==GL00_AUX_PROBLEM)
    {
      sqstructureC = new TSquareStructure2D(convolution_space);
      sqstructureC->Sort();
      sqmatrixGL00AuxProblem = new  TSquareMatrix2D(sqstructureC);
      MatricesGL00AuxProblem[i] =  sqmatrixGL00AuxProblem;
    }

    N_Unknowns = 2*N_U + N_P;

    OutPut("dof velocity : "<< setw(10) << 2* N_U << endl);
    OutPut("dof pressure : "<< setw(10) << N_P << endl);
    OutPut("dof all      : "<<  setw(10) << N_Unknowns  << endl);
    // matrices for VMS_PROJECTION
    if (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION)
    {
      if (TDatabase::ParamDB->VMS_LARGE_VELOCITY_SPACE==0)
        projection_space = new TFESpace2D(coll,NameString, UString, BoundCondition,
          DiscP_PSpace,0, mortarcoll);
      else
        projection_space = new TFESpace2D(coll,NameString, UString, BoundCondition,
          DiscP_PSpace,1, mortarcoll);

      ProjectionSpaces[i] = projection_space;
      sqstructureL = new TSquareStructure2D(projection_space);
      sqstructureL->Sort();
      structure_tilde_G = new TStructure2D(velocity_space, projection_space);
      structure_G = new TStructure2D(projection_space, velocity_space);
      sqmatrixL = new TSquareMatrix2D(sqstructureL);
      MatricesL[i] = sqmatrixL;
      matrix_tilde_G11 = new TMatrix2D(structure_tilde_G);
      Matrices_tilde_G11[i] = matrix_tilde_G11;
      matrix_tilde_G22 = new TMatrix2D(structure_tilde_G);
      Matrices_tilde_G22[i] = matrix_tilde_G22;
      matrix_G11 = new TMatrix2D(structure_G);
      Matrices_G11[i] = matrix_G11;
      matrix_G22 = new TMatrix2D(structure_G);
      Matrices_G22[i] = matrix_G22;
      N_L = projection_space->GetN_DegreesOfFreedom();
      OutPut("dof projection : " << setw(10) << N_L << endl);
    }

    rhs = new double[N_Unknowns];
    memset(rhs, 0, N_Unknowns*SizeOfDouble);
    RhsArray[i] = rhs;
    sol = new double[N_Unknowns];
    oldsol = new double[N_Unknowns];
    memset(sol, 0, N_Unknowns*SizeOfDouble);
    memset(oldsol, 0, N_Unknowns*SizeOfDouble);
    if (i==mg_level-1)
      sol_timestep_m1 =  new double[N_Unknowns];

    B = new double [N_Unknowns];

    // ( A B' )
    // ( B 0  )

    if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
    {
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
        if (i==mg_level-1)
        {
          alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
          alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_FINE_SADDLE;
        }
        else
        {
          alpha[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
          alpha[1] = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
        }
      }

      downwind[i] = new int[coll->GetN_Cells()];
      for (j=0;j<coll->GetN_Cells();j++)
        downwind[i][j] = j;
      #ifdef __DOWNWIND__
      DownwindNumberingCells(coll, downwind[i]);
      #endif

      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          MGLevel = new TNSE_MGLevel1(i, sqmatrixM, matrixB1, matrixB2,
            structureBT, B, sol, n_aux,  alpha,
            velocity_space_code ,
            pressure_space_code,NULL,downwind[i]);
          break;

        case 2:
          MGLevel = new TNSE_MGLevel2(i, sqmatrixM, matrixB1, matrixB2,
            matrixB1T, matrixB2T,
            B, sol, n_aux, alpha,
            velocity_space_code ,
            pressure_space_code,NULL,downwind[i]);
          break;

        case 3:
          MGLevel = new TNSE_MGLevel3(i, sqmatrixM11, sqmatrixM12,
            sqmatrixM21, sqmatrixM22,
            matrixB1, matrixB2,
            structureBT,
            B, sol, n_aux, alpha,
            velocity_space_code ,
            pressure_space_code,NULL,downwind[i]);
          break;

        case 4:
          MGLevel = new TNSE_MGLevel4(i, sqmatrixM11,  sqmatrixM12,
            sqmatrixM21, sqmatrixM22,
            matrixB1, matrixB2,
            matrixB1T, matrixB2T,
            B, sol, n_aux, alpha,
            velocity_space_code ,
            pressure_space_code,coll,downwind[i]);
          break;
      }                          // end switch(NSTYPE)
      MG->AddLevel(MGLevel);
    }

    u = new TFEVectFunct2D(velocity_space, UString,  UString,  sol, N_U, 2);
    u1 = u->GetComponent(0);
    u2 = u->GetComponent(1);
    p = new TFEFunction2D(pressure_space, PString,  PString,  sol+2*N_U, N_P);

    U1Array[i] = u1;
    U2Array[i] = u2;
    PArray[i] = p;
    UArray[i] = u;

    u1->Interpolate(InitialU1);
    u2->Interpolate(InitialU2);
    p->Interpolate(InitialP);

    RHSs[0] = rhs;
    RHSs[1] = rhs + N_U;
    RHSs[2] = rhs + 2*N_U;

    memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);

    if (TDatabase::ParamDB->DISCTYPE==CLASSICAL_LES)
    {
      if (i==mg_level-1)
      {
        LESModelRhs =  new double[2*N_U];
        memset(LESModelRhs,0,2*N_U*SizeOfDouble);
      }
    }

    if (TDatabase::ParamDB->DISCTYPE==GL00_CONVOLUTION)
    {
      duConvSpaces[i] = convolution_space;
      // define vector fe function for convolution
      N_UConv = convolution_space->GetN_DegreesOfFreedom();
      // allocate memory for values of convoluted function
      u_conv = new double[3*N_UConv];
      // initialize u_conv to 0
      memset(u_conv,0,3*N_UConv*SizeOfDouble);
      // allocate fe vector function for convolution
      duConv = new TFEVectFunct2D(convolution_space, UConvString, UConvString,
        u_conv, N_UConv, 3);
      // array for duConv for all levels
      duConvArray[i] = duConv;

      // copy the vector fe function to a fe function (only pointers)
      du1Conv = duConv->GetComponent(0);
      du2Conv = duConv->GetComponent(1);
      du3Conv = duConv->GetComponent(2);
      du1ConvArray[i] = du1Conv;
      du2ConvArray[i] = du2Conv;
      du3ConvArray[i] = du3Conv;
      if (i==mg_level-1)
      {
        LESModelRhs =  new double[2*N_U];
        memset(LESModelRhs,0,2*N_U*SizeOfDouble);
      }
    }

    if  (TDatabase::ParamDB->DISCTYPE==GL00_AUX_PROBLEM)
    {
      // allocate arrays for multigrid for auxiliary problem in
      // Galdi/Layton model

      solGL00AuxProblem =  new double[3*N_U];
      memset(solGL00AuxProblem,0,3*N_U*SizeOfDouble);
      rhsGL00AuxProblem =  new double[3*N_U];
      memset(rhsGL00AuxProblem,0,3*N_U*SizeOfDouble);
      if (i==mg_level-1)
      {
        LESModelRhs =  new double[2*N_U];
        memset(LESModelRhs,0,2*N_U*SizeOfDouble);
      }
      // build multigrid for auxiliary problem in Galdi/Layton model
      MGLevelGL00AuxProblem= new TMGLevel2D(i,sqmatrixGL00AuxProblem,
        rhsGL00AuxProblem, solGL00AuxProblem,
        n_aux, NULL);
      MGGL00AuxProblem->AddLevel(MGLevelGL00AuxProblem);

      // allocate fe vector function for solution of auxiliary problem
      GL00AuxProblemSol = new TFEVectFunct2D(convolution_space, AuxProbString, AuxProbString,
        solGL00AuxProblem, N_U, 3);
      // array for duConv for all levels
      GL00AuxProblemSolArray[i] = GL00AuxProblemSol;

      // copy the vector fe function to a fe function (only pointers)
      GL00AuxProblemSol11 = GL00AuxProblemSol->GetComponent(0);
      GL00AuxProblemSol12 = GL00AuxProblemSol->GetComponent(1);
      GL00AuxProblemSol22 = GL00AuxProblemSol->GetComponent(2);
      GL00AuxProblemSol11Array[i] = GL00AuxProblemSol11;
      GL00AuxProblemSol12Array[i] = GL00AuxProblemSol12;
      GL00AuxProblemSol22Array[i] = GL00AuxProblemSol22;
    }

    // turbulent viscosity involving the convolution of the solution
    if ((TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE==4)||
      (TDatabase::ParamDB->CONVOLUTE_SOLUTION))
    {
      // use same space for convolution and velocity
      uConvSpaces[i] = convolution_space;
      // define vector fe function for convolution
      // allocate memory for values of convoluted function
      u_uConv = new double[2*N_U];
      // initialize u_conv to 0
      memset(u_uConv,0,2*N_U*SizeOfDouble);
      // allocate fe vector function for convolution
      uConv = new TFEVectFunct2D(convolution_space, UConvString, UConvString,
        u_uConv, N_U, 2);
      // array for uConv for all levels
      uConvArray[i] = uConv;

      // copy the vector fe function to a fe function (only pointers)
      u1Conv = uConv->GetComponent(0);
      u2Conv = uConv->GetComponent(1);
      u1ConvArray[i] = u1Conv;
      u2ConvArray[i] = u2Conv;
      // compute matrix for auxiliary problem if not yet done
      if (TDatabase::ParamDB->DISCTYPE!=GL00_AUX_PROBLEM)
      {
        if (i==mg_level-1)
        {
          sqstructureC = new TSquareStructure2D(convolution_space);
          sqstructureC->Sort();
          sqmatrixGL00AuxProblem = new  TSquareMatrix2D(sqstructureC);
          MatricesGL00AuxProblem[i] =  sqmatrixGL00AuxProblem;
          rhsGL00AuxProblem =  new double[2*N_U];
          memset(rhsGL00AuxProblem,0,2*N_U*SizeOfDouble);
          // assemble matrix
          DiscreteForm = DiscreteFormMatrixAuxProblemU;
          fesp[0] = convolution_space;
          fesp[1] = velocity_space;

          fefct[0] = u1;
          fefct[1] = u2;

          // assemble matrix
          N_FESpaces = 1;
          N_Rhs = 0;
          N_SquareMatrices = 1;
          N_RectMatrices = 0;

          SQMATRICES[0] = MatricesGL00AuxProblem[mg_level-1];
          SQMATRICES[0]->Reset();
          aux =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
            TimeNSN_ParamFct2,
            TimeNSN_FEValues2,
            fesp+1, fefct,
            TimeNSFct2,
            TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
            TimeNSN_Params2, TimeNSBeginParam2);

          Assemble2D(N_FESpaces, fesp,
            N_SquareMatrices, SQMATRICES,
            N_RectMatrices, MATRICES,
            N_Rhs, RHSs, ferhs,
            DiscreteForm,
            BoundaryConditionsAuxProblem,
            BoundValues,
            aux);
          delete aux;
        }
      }
    }

    #ifdef __BENCH__
    if (i==mg_level-1)
    {
      former_sol =  new double [2*N_U];
      memset(former_sol, 0, 2*N_U*SizeOfDouble);
      U1old = new TFEFunction2D(velocity_space, UString,  UString, former_sol, N_U);
      U2old = new TFEFunction2D(velocity_space, UString,  UString, former_sol+N_U, N_U);
    }
    #endif

    // set discrete forms
    if ((mg_type==1) && (i<mg_level-1))
    {
      DiscreteForm = DiscreteFormUpwind;
      CurrentDiscType =  UPWIND;
    }
    else
      switch(TDatabase::ParamDB->DISCTYPE)
      {
        case GALERKIN:
          DiscreteForm = DiscreteFormGalerkin;
          CurrentDiscType =  GALERKIN;
          break;

        case UPWIND:
          DiscreteForm = DiscreteFormUpwind;
          CurrentDiscType =  UPWIND;
          break;

        case SMAGORINSKY:
          DiscreteForm = DiscreteFormSmagorinsky;
          CurrentDiscType =  SMAGORINSKY;
          break;

        case  CLASSICAL_LES:
          DiscreteForm = DiscreteFormColetti;
          CurrentDiscType =  CLASSICAL_LES ;
          very_first_time=1;
          break;

        case  GL00_CONVOLUTION:
          DiscreteForm = DiscreteFormGL00Convolution;
          CurrentDiscType =  GL00_CONVOLUTION;
          very_first_time=1;
          break;

        case  GL00_AUX_PROBLEM:
          DiscreteForm = DiscreteFormGL00AuxProblem;
          CurrentDiscType =  GL00_AUX_PROBLEM;
          very_first_time=1;
          break;

        case VMS_PROJECTION:
          DiscreteForm = DiscreteFormVMSProjection;
          CurrentDiscType =  VMS_PROJECTION;
          break;

        default:
          OutPut("Unknown DISCTYPE" << endl);
          exit(1);
      }

    // parameters which are the same for all NSTYPEs
    N_Rhs = 2;
    N_FESpaces = 3;

    // set matrices
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        SQMATRICES[0] = MatricesA[i];
        SQMATRICES[1] = MatricesM[i];
        MATRICES[0] = MatricesB1[i];
        MATRICES[1] = MatricesB2[i];

        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        MATRICES[0]->Reset();
        MATRICES[1]->Reset();

        N_SquareMatrices = 2;
        N_RectMatrices = 2;

        if (CurrentDiscType ==  GL00_AUX_PROBLEM)
        {
          N_SquareMatrices = 3;
          SQMATRICES[2] = MatricesGL00AuxProblem[i];
          SQMATRICES[2]->Reset();
        }
        break;

      case 2:
        SQMATRICES[0] = MatricesA[i];
        SQMATRICES[1] = MatricesM[i];
        MATRICES[0] = MatricesB1[i];
        MATRICES[1] = MatricesB2[i];
        MATRICES[2] = MatricesB1T[i];
        MATRICES[3] = MatricesB2T[i];

        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        MATRICES[0]->Reset();
        MATRICES[1]->Reset();
        MATRICES[2]->Reset();
        MATRICES[3]->Reset();

        N_SquareMatrices = 2;
        N_RectMatrices = 4;

        if (CurrentDiscType ==  GL00_AUX_PROBLEM)
        {
          N_SquareMatrices = 3;
          SQMATRICES[2] = MatricesGL00AuxProblem[i];
          SQMATRICES[2]->Reset();
        }
        break;

      case 3:
        SQMATRICES[0] = MatricesA11[i];
        SQMATRICES[1] = MatricesA12[i];
        SQMATRICES[2] = MatricesA21[i];
        SQMATRICES[3] = MatricesA22[i];
        SQMATRICES[4] = MatricesM11[i];
        SQMATRICES[5] = MatricesM22[i];
        MATRICES[0] = MatricesB1[i];
        MATRICES[1] = MatricesB2[i];

        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        SQMATRICES[2]->Reset();
        SQMATRICES[3]->Reset();
        SQMATRICES[4]->Reset();
        SQMATRICES[5]->Reset();
        MATRICES[0]->Reset();
        MATRICES[1]->Reset();

        N_SquareMatrices = 6;
        N_RectMatrices = 2;

        if (CurrentDiscType ==  GL00_AUX_PROBLEM)
        {
          N_SquareMatrices = 7;
          SQMATRICES[6] = MatricesGL00AuxProblem[i];
          SQMATRICES[6]->Reset();
        }
        break;

      case 4:
        SQMATRICES[0] = MatricesA11[i];
        SQMATRICES[1] = MatricesA12[i];
        SQMATRICES[2] = MatricesA21[i];
        SQMATRICES[3] = MatricesA22[i];
        SQMATRICES[4] = MatricesM11[i];
        SQMATRICES[5] = MatricesM22[i];
        MATRICES[0] = MatricesB1[i];
        MATRICES[1] = MatricesB2[i];
        MATRICES[2] = MatricesB1T[i];
        MATRICES[3] = MatricesB2T[i];

        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        SQMATRICES[2]->Reset();
        SQMATRICES[3]->Reset();
        SQMATRICES[4]->Reset();
        SQMATRICES[5]->Reset();
        MATRICES[0]->Reset();
        MATRICES[1]->Reset();
        MATRICES[2]->Reset();
        MATRICES[3]->Reset();

        N_SquareMatrices = 6;
        N_RectMatrices = 4;

        if (CurrentDiscType ==  GL00_AUX_PROBLEM)
        {
          N_SquareMatrices = 7;
          SQMATRICES[6] = MatricesGL00AuxProblem[i];
          SQMATRICES[6]->Reset();
        }
        if (CurrentDiscType == VMS_PROJECTION)
        {
          N_SquareMatrices = 7;
          SQMATRICES[6] =  MatricesL[i];
          SQMATRICES[6]->Reset();
          N_RectMatrices = 8;
          MATRICES[4] = Matrices_tilde_G11[i];
          MATRICES[5] = Matrices_tilde_G22[i];
          MATRICES[6] = Matrices_G11[i];
          MATRICES[7] = Matrices_G22[i];
          MATRICES[4]->Reset();
          MATRICES[5]->Reset();
          MATRICES[6]->Reset();
          MATRICES[7]->Reset();
          N_FESpaces = 4;
          fesp[3] = projection_space;
        }
        break;
    }

    // prepare output
    if (i==mg_level-1)
    {
      // kurz weg Output = new TOutput2D(4, 2, 2, 1, Domain);
      Output = new TOutput2D(5, 4, 1, 1, Domain);
      Output->AddFEVectFunct(u);
      Output->AddFEFunction(p);
      os.seekp(std::ios::beg);
      Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());

      psi = new double[N_V];
      StreamFct = new TFEFunction2D(streamfunction_space, PsiString,  PsiString,  psi, N_V);
//       Output->AddFEFunction(StreamFct);

      vorticity = new double[N_Vort];
      Vorticity = new TFEFunction2D(vorticity_space, VorticityString, VorticityString, vorticity, N_Vort);
      
//       ComputeVorticity(U1Array[mg_level-1], U2Array[mg_level-1],Vorticity);
      
      Output->AddFEFunction(Vorticity);

      if ( TDatabase::ParamDB->CONVOLUTE_SOLUTION && mixing_layer)
      {
        conv_vort = new double[N_Vort];
        Conv_Vort = new TFEFunction2D(vorticity_space, ConvVortString, ConvVortString, conv_vort, N_Vort);
      }
      div = new double[N_Vort];
      Divergence = new TFEFunction2D(vorticity_space, DivergenceString, DivergenceString, div, N_Vort);
//       Output->AddFEFunction(Divergence);

      app = new double[N_V*2];
      uconf = new TFEVectFunct2D(streamfunction_space, UConfString, UConfString, app, N_V, 2);
    }

    // read initial solution of finest level from grape file
    if ((i==mg_level-1)&&(TDatabase::ParamDB->READ_GRAPE_FILE))
    {
      // version without vorticity and divergence
      /*AuxFEFunctArray = new TFEFunction2D*[2];
      AuxFEFunctArray[0] = PArray[mg_level-1];
      AuxFEFunctArray[1] = StreamFct;
      AuxFEVectFunctArray = new TFEVectFunct2D*[2];
      AuxFEVectFunctArray[0] = UArray[mg_level-1];
      AuxFEVectFunctArray[1] = uconf;
      ReadGrapeFile(ReadGrapeBaseName, 2 , 2, AuxFEFunctArray,AuxFEVectFunctArray);
      */
      // version with vorticity and divergence
      AuxFEFunctArray = new TFEFunction2D*[4];
      AuxFEFunctArray[0] = PArray[mg_level-1];
      AuxFEFunctArray[1] = StreamFct;
      AuxFEFunctArray[2] = Vorticity;
      AuxFEFunctArray[3] = Divergence;
      AuxFEVectFunctArray = new TFEVectFunct2D*[1];
      AuxFEVectFunctArray[0] = UArray[mg_level-1];
      //AuxFEVectFunctArray[1] = uconf;
      ReadGrapeFile(ReadGrapeBaseName, 4 , 1, AuxFEFunctArray,AuxFEVectFunctArray);
      // right hand side of first time step will not be completely correct
      // for GL00--discs and Coletti--disc
      // since to reconstruct the rhs of former time step, the initial
      // velocity of this time step is necessary
      // instead final velocity is used now
      if (TDatabase::TimeDB->RESET_CURRENTTIME > 0)
      {
        TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->RESET_CURRENTTIME_STARTTIME;
        OutPut("start time reset to " << TDatabase::TimeDB->CURRENTTIME << endl);
      }
    }

    // read initial solution of finest level from grape file
    if ((i==mg_level-1)&&(TDatabase::ParamDB->READ_DATA))
    {
	save_sol[0] = sol;
	save_N_Unknowns[0] = N_Unknowns;
	ReadData(ReadDataFileName,1,save_sol,save_N_Unknowns);
	//for (j=0;j<N_P;j++)
	//  sol[3*N_U+j] = 0;
    }

    // set rhs
    fesp[0] = velocity_space;
    fesp[1] = pressure_space;
    fesp[2] = convolution_space;

    fefct[0] = u1;
    fefct[1] = u2;
    fefct[2] = du1ConvArray[i];
    fefct[3] = du2ConvArray[i];
    fefct[4] = du3ConvArray[i];

    ferhs[0] = velocity_space;
    ferhs[1] = velocity_space;

    // 3 parameters are needed for assembling
    // which are u1_old, u2_old, norm of grad u_old
    switch(CurrentDiscType)
    {
      // turbulent viscosity must be computed
      case SMAGORINSKY:
      case VMS_PROJECTION:
      case CLASSICAL_LES:
      case GL00_CONVOLUTION:
      case GL00_AUX_PROBLEM:
        aux =  new TAuxParam2D(TimeNSN_FESpacesVelo_GradVelo, TimeNSN_FctVelo_GradVelo,
          TimeNSN_ParamFctVelo_GradVelo,
          TimeNSN_FEValuesVelo_GradVelo,
          fesp, fefct,
          TimeNSFctVelo_GradVelo,
          TimeNSFEFctIndexVelo_GradVelo,
          TimeNSFEMultiIndexVelo_GradVelo,
          TimeNSN_ParamsVelo_GradVelo,
          TimeNSBeginParamVelo_GradVelo);

        break;

      default:
        // 2 parameters are needed for assembling (u1_old, u2_old)
        aux =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
          TimeNSN_ParamFct2,
          TimeNSN_FEValues2,
          fesp, fefct,
          TimeNSFct2,
          TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
          TimeNSN_Params2, TimeNSBeginParam2);
    }

    //======================================================================
    // assembling of matrices for each level
    // A_11 , (A_12), (A_21), (A_22), M_11, (M_22)
    // assembling of rhs not needed at this point
    //======================================================================
    Assemble2D(N_FESpaces, fesp,
      N_SquareMatrices, SQMATRICES,
      N_RectMatrices, MATRICES,
      N_Rhs, RHSs, ferhs,
      DiscreteForm,
      BoundaryConditions,
      BoundValues,
      aux);
    delete aux;

    // copy Dirichlet values from rhs into sol
    memcpy(sol+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
    memcpy(sol+N_Active+N_U, RHSs[1]+N_Active, (N_U-N_Active)*SizeOfDouble);

    // add convective term in the upwind discretizations
    if(DiscreteForm == DiscreteFormUpwind)
    {
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
        case 2:
          // do upwinding with one matrix
          UpwindForNavierStokes(LinCoeffs, SQMATRICES[0], U1Array[i], U2Array[i]);
          //cout << "UPWINDING DONE : level " << i << endl;
          break;

        case 3:
        case 4:
          // do upwinding with two matrices
          UpwindForNavierStokes(LinCoeffs, SQMATRICES[0], U1Array[i], U2Array[i]);
          UpwindForNavierStokes(LinCoeffs, SQMATRICES[3], U1Array[i], U2Array[i]);
          cout << "UPWINDING DONE(2) : level " << i << endl;
          break;
      }                          // endswitch
    }

    // set rows of Dirichlet dof in off diagonal matrix blocks
    // to zero
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 3:
      case 4:
        // get number of active dof
        N_Active = velocity_space->GetActiveBound();
        // get row in off diagonal matrix where the Dirichlet nodes start
        RowPtr = MatricesA12[i]->GetRowPtr();
        // compute number of entries starting from this row to the end
        // of the matrix
        j = RowPtr[N_Active];
        k = RowPtr[N_U]-j;
        // set these entries to zero
        memset(MatricesA12[i]->GetEntries()+j, 0, SizeOfDouble*k);
        memset(MatricesA21[i]->GetEntries()+j, 0, SizeOfDouble*k);
        break;
    }
    // update matrices
    if (CurrentDiscType == VMS_PROJECTION)
    {
      SQMATRICES[0] = MatricesA11[i];
      SQMATRICES[1] = MatricesA12[i];
      SQMATRICES[2] = MatricesA21[i];
      SQMATRICES[3] = MatricesA22[i];
      SQMATRICES[6] =  MatricesL[i];
      MATRICES[2] = Matrices_tilde_G11[i];
      MATRICES[3] = Matrices_tilde_G22[i];
      MATRICES[4] = Matrices_G11[i];
      MATRICES[5] = Matrices_G22[i];

      VMSProjectionUpdateMatrices(N_U,N_Active,N_L,SQMATRICES,MATRICES);
    }
  }                              // endfor i
  t4 =  GetTime();
  total_time += t4 - t3;
  t3 = t4;
  //======================================================================
  // end of space cycle, finest grid reached
  // everything happens on the same grid
  //======================================================================

  // copy sol for extrapolation after time step
  memcpy(sol_timestep_m1,sol,N_Unknowns*SizeOfDouble);

  comp_vort =0;
  if (TDatabase::ParamDB->WRITE_GRAPE|| TDatabase::ParamDB->WRITE_GNU ||
      TDatabase::ParamDB->WRITE_GMV || TDatabase::ParamDB->WRITE_VTK  )
  {
    StreamFunction(USpaces[mg_level-1], sol, sol+N_Uarray[mg_level-1],
      PsiSpaces[LEVELS-1], psi);
    Prolongate(USpaces[mg_level-1], PsiSpaces[LEVELS-1], sol,
      uconf->GetValues(), app);
    Prolongate(USpaces[mg_level-1], PsiSpaces[LEVELS-1],
      sol+N_Uarray[mg_level-1], uconf->GetValues()+N_V, app);
    ComputeVorticityDivergence(USpaces[mg_level-1], U1Array[mg_level-1], U2Array[mg_level-1],
      vorticity_space,vorticity,div);
    comp_vort++;
  }

  #ifdef __MIXINGLAYERSLIP__
  if (!comp_vort)
    ComputeVorticityDivergence(USpaces[mg_level-1], U1Array[mg_level-1], U2Array[mg_level-1],
      vorticity_space,vorticity,div);
  ComputeVorticiyThickness(Vorticity,errors);
  OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
  vort_zero = errors[0];
  OutPut( "vorticity thickness: " << errors[0] << " " <<
    errors[0]/vort_zero << endl);
  comp_vort =0;

  // enstrophy
  aux =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
    TimeNSN_ParamFct2,
    TimeNSN_FEValues2,
    fesp, fefct,
    TimeNSFct2,
    TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
    TimeNSN_Params2, TimeNSBeginParam2);
  Vorticity->GetErrors(ExactNull, 3, TimeNSAllDerivatives,
    2, L2H1Errors,
    NULL, aux, 1,  &vorticity_space, errors);

  OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
  OutPut( "enstrophy " << (errors[0]*errors[0])/2 << endl);
  OutPut(endl);
  delete aux;
  #endif

  if (TDatabase::ParamDB->CONVOLUTE_SOLUTION)
  {
    // prepare auxiliary problem
    MG->RestrictToAllGrids();
    DiscreteForm = DiscreteFormRHSAuxProblemU;
    fesp[0] = USpaces[mg_level-1];

    fefct[0] = U1Array[mg_level-1];
    fefct[1] = U2Array[mg_level-1];

    ferhs[0] = USpaces[mg_level-1];
    ferhs[1] = USpaces[mg_level-1];

    // assemble rhs
    N_FESpaces = 1;
    N_Rhs = 2;
    N_SquareMatrices = 0;
    N_RectMatrices = 0;

    memset(rhsGL00AuxProblem,0,2*N_U*SizeOfDouble);
    RHSs[0] = rhsGL00AuxProblem;
    RHSs[1] = rhsGL00AuxProblem + N_U;

    aux =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
      TimeNSN_ParamFct2,
      TimeNSN_FEValues2,
      fesp, fefct,
      TimeNSFct2,
      TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
      TimeNSN_Params2, TimeNSBeginParam2);

    Assemble2D(N_FESpaces, fesp,
      N_SquareMatrices, SQMATRICES,
      N_RectMatrices, MATRICES,
      N_Rhs, RHSs, ferhs,
      DiscreteForm,
      BoundaryConditionsAuxProblem,
      BoundValuesAuxProblem,
      aux);
    // solve auxiliary problem
    Solver(sqmatrixGL00AuxProblem, RHSs[0], u_uConv,2);

    // error in first component
    u1ConvArray[mg_level-1]->GetErrors(ExactNull, 3, TimeNSAllDerivatives,
      2, L2H1Errors,
      NULL, aux, 1, USpaces+mg_level-1, errors);

    // error in second component
    u2ConvArray[mg_level-1]->GetErrors(ExactNull, 3, TimeNSAllDerivatives,
      2, L2H1Errors,
      NULL, aux, 1, USpaces+mg_level-1, errors+2);

    OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
    OutPut( "conv H1-semi(u): " << sqrt(errors[1]*errors[1]+errors[3]*errors[3]
      ) << endl);
    OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
    OutPut( "conv kinetic energy " << (errors[0]*errors[0]+errors[2]*errors[2]
      )/2<< endl);
    delete aux;

    #ifdef __MIXINGLAYERSLIP__
    // compute vorticity
    ComputeVorticityDivergence(velocity_space, u1ConvArray[mg_level-1],
      u2ConvArray[mg_level-1],
      vorticity_space,conv_vort,div);

    ComputeVorticiyThickness(Conv_Vort,errors);
    OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
    vort_zero_conv =  errors[0];
    OutPut( "vorticity thickness (uconv): " << errors[0] << " " <<
      errors[0]/vort_zero_conv << endl);

    // enstrophy
    aux =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
      TimeNSN_ParamFct2,
      TimeNSN_FEValues2,
      fesp, fefct,
      TimeNSFct2,
      TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
      TimeNSN_Params2, TimeNSBeginParam2);
    Conv_Vort->GetErrors(ExactNull, 3, TimeNSAllDerivatives,
      2, L2H1Errors,
      NULL, aux, 1,  &vorticity_space, errors);

    OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
    OutPut( "enstrophy (uconv) " << (errors[0]*errors[0])/2 << endl);
    delete aux;
    // save convolved vorticity on div for output
    for (k=0;k<N_Vort;k++)
      div[k] = conv_vort[k];
    #endif

  }

  if(TDatabase::ParamDB->WRITE_GRAPE)
  {
    os.seekp(std::ios::beg);
    os << GrapeBaseName << 0 << ".dat" << ends;
    Output->WriteGrape(os.str().c_str());
  }

  if(TDatabase::ParamDB->WRITE_GMV)
  {
      os.seekp(std::ios::beg);
      os << GmvBaseName << 0 << ".gmv" << ends;
      Output->WriteGMV(os.str().c_str());
  }

  if(TDatabase::ParamDB->WRITE_GNU)
  {
    os.seekp(std::ios::beg);
    os << GnuBaseName << 0 << ".gnu" << ends;
    Output->WriteGnuplot(os.str().c_str());

    os.seekp(std::ios::beg);
    os << GnuBaseName << 0 << ".psi" << ends;
    Output->WriteGNU_iso(os.str().c_str(),1);

    N_GNU_images++;
  }

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

  if ((TDatabase::ParamDB->WRITE_GRAPE)||
      (TDatabase::ParamDB->WRITE_GMV) ||
      (TDatabase::ParamDB->WRITE_VTK))
      N_GRAPE_images++;

  defect = new double[N_Unknowns];
  startsol = new double[N_Unknowns];
  frac_step_sol = new double[N_Unknowns];
  oldrhs =  new double[N_Unknowns];

  N_Active = velocity_space->GetActiveBound();

  solver_time = 0.0;
  N_LinIter = 0;

  gamma = 0;
  m = 0;
  N_SubSteps = GetN_SubSteps();
  oldtau = 1;
  end_time = TDatabase::TimeDB->ENDTIME;
  if (TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL)
    time_discs = 2;
  else
    time_discs = 1;
  // initialize solver
  if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
  {
    switch (TDatabase::ParamDB->SC_SOLVER_SADDLE)
    {
      case 11:
        zerostart = 1;
        break;
      case 16:
        zerostart = 0;
        break;
    }
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
        }
        else
        {
          itmethod_sol = sol;
          itmethod_rhs = rhs;
        }
        break;
      case 16:
        itmethod = new TFgmresIte(MatVect, Defect, prec,
          0, N_Unknowns, 0);
        if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
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
  if (fabs(hmin-hmax)<1e-6)
  {
    OutPut("delta " <<  TDatabase::ParamDB->FILTER_WIDTH_CONSTANT *
      pow(hmin,TDatabase::ParamDB->FILTER_WIDTH_POWER) << endl);
  }
  else
  {
    if (TDatabase::ParamDB->FILTER_WIDTH_POWER==0)
    {
      OutPut("delta " <<  TDatabase::ParamDB->FILTER_WIDTH_CONSTANT << endl);
    }
    else
    {
      OutPut("delta is non--constant" << endl);
    }
  }

  //======================================================================
  // start of time cycle
  // everything happens on the same grid
  //======================================================================

  while(TDatabase::TimeDB->CURRENTTIME< end_time)
  {
    // time cycle
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

    for (methods=0;methods<time_discs;methods++)
    {
      if (time_discs==2)
      {
        if (methods==0)          // fractional-step-theta-scheme
        {
          TDatabase::TimeDB->TIME_DISC = 3;
                                 // save start sol (nec. for gl00)
          memcpy(startsol,sol,N_Unknowns*SizeOfDouble);
                                 // save rhs
          memcpy(oldrhs,rhs,N_Unknowns*SizeOfDouble);
        }
        else                     // crank nicolson scheme
        {                        // take solution of first scheme as initial iterate
          TDatabase::TimeDB->TIME_DISC = 2;
          TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->INTERNAL_STARTTIME;
                                 // save solution of fract.step
          memcpy(frac_step_sol,sol,N_Unknowns*SizeOfDouble);
                                 // get former startsol
          memcpy(sol,startsol,N_Unknowns*SizeOfDouble);
                                 // get old rhs
          memcpy(rhs,oldrhs,N_Unknowns*SizeOfDouble);
        }
        N_SubSteps = GetN_SubSteps();
      }

      for(l=0;l<N_SubSteps;l++)  // sub steps of fractional step theta
      {
        if (!very_first_time)
        {
          SetTimeDiscParameters();
	  press_factor1 = 1.0;
	  // set parameters and change matrices for treating
	  // the pressure in the same way as the velocity
	  // in the temporal discretization
	  if (TDatabase::TimeDB->EXTRAPOLATE_PRESSURE==1)
          {
	      press_factor1 = TDatabase::TimeDB->THETA1;
	    press_factor2 = TDatabase::TimeDB->THETA2;
	  }  // end for EXTRAPOLATE_PRESSURE == 1
        }
        if (m==1)
        {
          OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
          OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
          OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
          OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
	  press_factor1 = 1.0;
	  press_factor2 = 0.0;
        }

        tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
        if (!very_first_time)
          TDatabase::TimeDB->CURRENTTIME += tau;
	// only for special discretizations
        if (very_first_time)
          oldtau=tau;
        // working array for rhs is B, initialize B
        memset(B, 0, N_Unknowns*SizeOfDouble);

        // old rhs multiplied with current subtime step and theta3 on B
        Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, rhs, B);
        Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, rhs+N_U, B+N_U);

        //======================================================================
        // prepare input data (parameters) for several discretizations

        // compute convolution of \nabla u \nabla u^T
        if (TDatabase::ParamDB->DISCTYPE==GL00_CONVOLUTION)
        {
          level_down =0;
          // compute convolution
          MG->RestrictToAllGrids();

          // convolute tensor
          // ConvoluteSymmetricTensorFull(UArray[mg_level-1-level_down], duConv);
          ConvoluteSymmetricTensor(UArray[mg_level-1-level_down], duConv);

          auxConv = new double[N_Unknowns];
          // prolongate convoluted function
          for (ii=mg_level-1-level_down ; ii< mg_level-1; ii++)
          {
            Prolongate(duConvSpaces[ii], duConvSpaces[ii+1],
              du1ConvArray[ii]->GetValues(),
              du1ConvArray[ii+1]->GetValues(),
              auxConv);
            Prolongate(duConvSpaces[ii], duConvSpaces[ii+1],
              du2ConvArray[ii]->GetValues(),
              du2ConvArray[ii+1]->GetValues(),
              auxConv);
            Prolongate(duConvSpaces[ii], duConvSpaces[ii+1],
              du3ConvArray[ii]->GetValues(),
              du3ConvArray[ii+1]->GetValues(),
              auxConv);
          }
          delete auxConv ;
        }

        // compute convolution of u for |u-g_\delta\ast(u)|
        if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE==4)
        {
          // prepare auxiliary problem
          MG->RestrictToAllGrids();
          DiscreteForm = DiscreteFormRHSAuxProblemU;
          fesp[0] = USpaces[mg_level-1];

          fefct[0] = U1Array[mg_level-1];
          fefct[1] = U2Array[mg_level-1];

          ferhs[0] = USpaces[mg_level-1];
          ferhs[1] = USpaces[mg_level-1];

          // assemble rhs
          N_FESpaces = 1;
          N_Rhs = 2;
          N_SquareMatrices = 0;
          N_RectMatrices = 0;

          memset(rhsGL00AuxProblem,0,2*N_U*SizeOfDouble);
          RHSs[0] = rhsGL00AuxProblem;
          RHSs[1] = rhsGL00AuxProblem + N_U;

          aux =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
            TimeNSN_ParamFct2,
            TimeNSN_FEValues2,
            fesp, fefct,
            TimeNSFct2,
            TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
            TimeNSN_Params2, TimeNSBeginParam2);

          Assemble2D(N_FESpaces, fesp,
            N_SquareMatrices, SQMATRICES,
            N_RectMatrices, MATRICES,
            N_Rhs, RHSs, ferhs,
            DiscreteForm,
            BoundaryConditionsAuxProblem,
            BoundValuesAuxProblem,
            aux);
          delete aux;
          // solve auxiliary problem
          OutPut("type 4 "<< endl);
          t1 = GetTime();
          Solver(sqmatrixGL00AuxProblem, RHSs[0], u_uConv,2);
          t2 = GetTime();
          OutPut( "time for AMG solving: " << t2-t1 << endl);
          for(ii = mg_level-1 ; ii > 0;ii--)
          {
            RestrictFunction(USpaces[ii-1], USpaces[ii],
              u1ConvArray[ii-1]->GetValues(),
              u1ConvArray[ii]->GetValues(),
              MG->GetLevel(ii-1)->GetAuxVector(0));
            RestrictFunction(USpaces[ii-1], USpaces[ii],
              u2ConvArray[ii-1]->GetValues(),
              u2ConvArray[ii]->GetValues(),
              MG->GetLevel(ii-1)->GetAuxVector(0));
          }
        }
        //======================================================================

        //======================================================================
        // assembling of the rhs of current sub time step
        // only on the finest level necessary
        // there is no assembling of matrices here
        //======================================================================

        fesp[0] = USpaces[mg_level-1];
        RHSs[0] = RhsArray[mg_level-1];
        RHSs[1] = RhsArray[mg_level-1]+N_Uarray[mg_level-1];
        ferhs[0] = USpaces[mg_level-1];
        ferhs[1] = USpaces[mg_level-1];

        // initialize array
        memset(RHSs[0], 0, (2*N_Uarray[mg_level-1]+N_Parray[mg_level-1])*SizeOfDouble);

        N_Rhs = 2;
        N_FESpaces = 1;
        N_SquareMatrices = 0;
        N_RectMatrices = 0;
        fefct[0] = U1Array[mg_level-1];
        fefct[1] = U2Array[mg_level-1];

        switch(TDatabase::ParamDB->DISCTYPE)
        {
          case CLASSICAL_LES :
            fesp[0] = USpaces[mg_level-1];

            fefct[0] = U1Array[mg_level-1];
            fefct[1] = U2Array[mg_level-1];

            ferhs[0] = USpaces[mg_level-1];
            ferhs[1] = USpaces[mg_level-1];

            DiscreteForm = DiscreteFormRHSColetti;
            if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE!=4)
            {
              aux =  new TAuxParam2D(TimeNSN_FESpacesVelo_GradVelo,
                TimeNSN_FctVelo_GradVelo,
                TimeNSN_ParamFctVelo_GradVelo,
                TimeNSN_FEValuesVelo_GradVelo,
                fesp, fefct,
                TimeNSFctVelo_GradVelo,
                TimeNSFEFctIndexVelo_GradVelo,
                TimeNSFEMultiIndexVelo_GradVelo,
                TimeNSN_ParamsVelo_GradVelo,
                TimeNSBeginParamVelo_GradVelo);
            }

            if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 4)
            {
              // convoluted velocity
              fesp[1] = uConvSpaces[mg_level-1];

              fefct[2] = u1ConvArray[mg_level-1];
              fefct[3] = u2ConvArray[mg_level-1];

              aux =  new TAuxParam2D(TimeNSN_FESpacesVelo_GradVeloNuT4,
                TimeNSN_FctVelo_GradVeloNuT4,
                TimeNSN_ParamFctVelo_GradVeloNuT4,
                TimeNSN_FEValuesVelo_GradVeloNuT4,
                fesp, fefct,
                TimeNSFctVelo_GradVeloNuT4,
                TimeNSFEFctIndexVelo_GradVeloNuT4,
                TimeNSFEMultiIndexVelo_GradVeloNuT4,
                TimeNSN_ParamsVelo_GradVeloNuT4,
                TimeNSBeginParamVelo_GradVeloNuT4);

            }

            Assemble2D(N_FESpaces, fesp,
              N_SquareMatrices, SQMATRICES,
              N_RectMatrices, MATRICES,
              N_Rhs, RHSs, ferhs,
              DiscreteForm,
              BoundaryConditions,
              BoundValues,
              aux);
            delete aux;

            memcpy(LESModelRhs, RHSs[0], N_U*SizeOfDouble);
            memcpy(LESModelRhs+N_U, RHSs[1], N_U*SizeOfDouble);

            // assemble rhs from f

            DiscreteForm = DiscreteFormRHS;
            aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
            // initialize array
            memset(RHSs[0], 0,
              (2*N_Uarray[mg_level-1]+N_Parray[mg_level-1])*SizeOfDouble);
            break;

          case GL00_CONVOLUTION :
            // current solution
            fesp[0] = USpaces[mg_level-1];

            fefct[0] = U1Array[mg_level-1];
            fefct[1] = U2Array[mg_level-1];

            // current rhs
            ferhs[0] = USpaces[mg_level-1];
            ferhs[1] = USpaces[mg_level-1];

            // convolution of grad u grad u^T
            fesp[1] = duConvSpaces[mg_level-1];

            fefct[2] = du1ConvArray[mg_level-1];
            fefct[3] = du2ConvArray[mg_level-1];
            fefct[4] = du3ConvArray[mg_level-1];

            // initialize array
            DiscreteForm = DiscreteFormRHSLESModel;
            if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE!=4)
            {
              aux =  new TAuxParam2D(TimeNSN_FESpacesGL00AuxProblem,
                TimeNSN_FctGL00AuxProblem,
                TimeNSN_ParamFctGL00AuxProblem,
                TimeNSN_FEValuesGL00AuxProblem,
                fesp, fefct,
                TimeNSFctGL00AuxProblem,
                TimeNSFEFctIndexGL00AuxProblem,
                TimeNSFEMultiIndexGL00AuxProblem,
                TimeNSN_ParamsGL00AuxProblem,
                TimeNSBeginParamGL00AuxProblem);
              /*  aux =  new TAuxParam2D(TimeNSN_FESpacesRHSGL00Convolution,
                                     TimeNSN_FctRHSGL00Convolution,
                                     TimeNSN_ParamFctRHSGL00Convolution,
                                     TimeNSN_FEValuesRHSGL00Convolution,
                                     fesp, fefct,
                                     TimeNSFctRHSGL00Convolution,
                                     TimeNSFEFctIndexRHSGL00Convolution,
                                     TimeNSFEMultiIndexRHSGL00Convolution,
                                     TimeNSN_ParamsRHSGL00Convolution,
                                     TimeNSBeginParamRHSGL00Convolution);*/
            }
            if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 17)
            {
              // convoluted velocity
              fesp[2] = uConvSpaces[mg_level-1];

              fefct[5] = u1ConvArray[mg_level-1];
              fefct[6] = u2ConvArray[mg_level-1];

              aux =  new TAuxParam2D(TimeNSN_FESpacesRHSGL00ConvolutionNuT4,
                TimeNSN_FctRHSGL00ConvolutionNuT4,
                TimeNSN_ParamFctRHSGL00ConvolutionNuT4,
                TimeNSN_FEValuesRHSGL00ConvolutionNuT4,
                fesp, fefct,
                TimeNSFctRHSGL00ConvolutionNuT4,
                TimeNSFEFctIndexRHSGL00ConvolutionNuT4,
                TimeNSFEMultiIndexRHSGL00ConvolutionNuT4,
                TimeNSN_ParamsRHSGL00ConvolutionNuT4,
                TimeNSBeginParamRHSGL00ConvolutionNuT4);
            }
            //N_FESpaces = 2;
            Assemble2D(N_FESpaces, fesp,
              N_SquareMatrices, SQMATRICES,
              N_RectMatrices, MATRICES,
              N_Rhs, RHSs, ferhs,
              DiscreteForm,
              BoundaryConditions,
              BoundValues,
              aux);

            delete aux;

            memcpy(LESModelRhs, RHSs[0], N_U*SizeOfDouble);
            memcpy(LESModelRhs+N_U, RHSs[1], N_U*SizeOfDouble);

            // assemble rhs from f

            DiscreteForm = DiscreteFormRHS;
            aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
            // initialize array
            memset(RHSs[0], 0,
              (2*N_Uarray[mg_level-1]+N_Parray[mg_level-1])*SizeOfDouble);
            N_FESpaces = 1;

            break;
          case GL00_AUX_PROBLEM :
            // assemble and solve auxiliary problem
            DiscreteForm = DiscreteFormGL00AuxProblemRHS;
            fesp[0] = USpaces[mg_level-1];

            fefct[0] = U1Array[mg_level-1];
            fefct[1] = U2Array[mg_level-1];

            ferhs[0] = USpaces[mg_level-1];
            ferhs[1] = USpaces[mg_level-1];
            ferhs[2] = USpaces[mg_level-1];

            N_FESpaces = 1;
            N_Rhs = 3;
            N_SquareMatrices = 0;
            N_RectMatrices = 0;

            memset(rhsGL00AuxProblem,0,3*N_U*SizeOfDouble);
            RHSs[0] = rhsGL00AuxProblem;
            RHSs[1] = rhsGL00AuxProblem+ N_Uarray[mg_level-1];
            RHSs[2] = rhsGL00AuxProblem+ 2*N_Uarray[mg_level-1];

            aux =  new TAuxParam2D(TimeNSN_FESpacesGrad, TimeNSN_FctGrad,
              TimeNSN_ParamFctGrad,
              TimeNSN_FEValuesGrad,
              fesp, fefct,
              TimeNSFctGrad,
              TimeNSFEFctIndexGrad,
              TimeNSFEMultiIndexGrad,
              TimeNSN_ParamsGrad,
              TimeNSBeginParamGrad);

            Assemble2D(N_FESpaces, fesp,
              N_SquareMatrices, SQMATRICES,
              N_RectMatrices, MATRICES,
              N_Rhs, RHSs, ferhs,
              DiscreteForm,
              BoundaryConditionsAuxProblem,
              BoundValuesAuxProblem,
              aux);

            delete aux;

            // solve auxiliary problem
            t1 = GetTime();
            Solver(sqmatrixGL00AuxProblem, RHSs[0], solGL00AuxProblem,3);
            t2 = GetTime();
            OutPut( "time for AMG solving: " << t2-t1 << endl);

            // assemble term coming from the LES model

            fefct[2] = GL00AuxProblemSol11Array[mg_level-1];
            fefct[3] = GL00AuxProblemSol12Array[mg_level-1];
            fefct[4] = GL00AuxProblemSol22Array[mg_level-1];

            N_Rhs = 2;
            N_SquareMatrices = 0;
            N_RectMatrices = 0;

            RHSs[0] = RhsArray[mg_level-1];
            RHSs[1] = RhsArray[mg_level-1]+N_Uarray[mg_level-1];

            // initialize array
            memset(RHSs[0], 0,
              (2*N_Uarray[mg_level-1]+N_Parray[mg_level-1])*SizeOfDouble);

            DiscreteForm = DiscreteFormRHSLESModel;
            if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE!=4)
            {
              aux =  new TAuxParam2D(TimeNSN_FESpacesGL00AuxProblem,
                TimeNSN_FctGL00AuxProblem,
                TimeNSN_ParamFctGL00AuxProblem,
                TimeNSN_FEValuesGL00AuxProblem,
                fesp, fefct,
                TimeNSFctGL00AuxProblem,
                TimeNSFEFctIndexGL00AuxProblem,
                TimeNSFEMultiIndexGL00AuxProblem,
                TimeNSN_ParamsGL00AuxProblem,
                TimeNSBeginParamGL00AuxProblem);
            }
            Assemble2D(N_FESpaces, fesp,
              N_SquareMatrices, SQMATRICES,
              N_RectMatrices, MATRICES,
              N_Rhs, RHSs, ferhs,
              DiscreteForm,
              BoundaryConditions,
              BoundValues,
              aux);
            delete aux;

            memcpy(LESModelRhs, RHSs[0], N_U*SizeOfDouble);
            memcpy(LESModelRhs+N_U, RHSs[1], N_U*SizeOfDouble);

            // assemble rhs from f

            DiscreteForm = DiscreteFormRHS;
            aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
            // initialize array
            memset(RHSs[0], 0,
              (2*N_Uarray[mg_level-1]+N_Parray[mg_level-1])*SizeOfDouble);
            break;

          case SMAGORINSKY_EXPL :
            fesp[0] = USpaces[mg_level-1];

            fefct[0] = U1Array[mg_level-1];
            fefct[1] = U2Array[mg_level-1];

            ferhs[0] = USpaces[mg_level-1];
            ferhs[1] = USpaces[mg_level-1];

            DiscreteForm = DiscreteFormRHSSmagorinskyExpl;

            aux =  new TAuxParam2D(TimeNSN_FESpacesVelo_GradVelo,
              TimeNSN_FctVelo_GradVelo,
              TimeNSN_ParamFctVelo_GradVelo,
              TimeNSN_FEValuesVelo_GradVelo,
              fesp, fefct,
              TimeNSFctVelo_GradVelo,
              TimeNSFEFctIndexVelo_GradVelo,
              TimeNSFEMultiIndexVelo_GradVelo,
              TimeNSN_ParamsVelo_GradVelo,
              TimeNSBeginParamVelo_GradVelo);

            break;
          default:
            DiscreteForm = DiscreteFormRHS;
            aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
            break;
        }

        // all input data are computed, assemble now rhs
        Assemble2D(N_FESpaces, fesp,
          N_SquareMatrices, SQMATRICES,
          N_RectMatrices, MATRICES,
          N_Rhs, RHSs, ferhs,
          DiscreteForm,
          BoundaryConditions,
          BoundValues,
          aux);

        if (very_first_time==1)
        {
          very_first_time=0;
          l--;
          continue;
        }
        // add rhs from current sub time step to rhs array B
        Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, rhs, B);
        Daxpy(N_Active, tau*TDatabase::TimeDB->THETA4, rhs+N_U, B+N_U);
        if ((TDatabase::ParamDB->DISCTYPE==CLASSICAL_LES)||
          (TDatabase::ParamDB->DISCTYPE==GL00_AUX_PROBLEM)||
          (TDatabase::ParamDB->DISCTYPE==GL00_CONVOLUTION))
        {
          Daxpy(N_U, tau, LESModelRhs, B );
          Daxpy(N_U, tau, LESModelRhs+N_U, B+N_U);
        }

        // do not change rhs during nonlinear iteration !!!
        // needed for next time step !!!

        // slip type bc detected, manipulation of matrices is necessary
        // this is done only at the very beginning
        // the matrices A_12, A_12, M_11, M_12, M_21, M_22, B1T, B2T
        //     stay unchanged during the complete solution process
        // the matrices A_11, A_22 are manipulated after their new
        //     assembling during the nonlinear iteration

        if ((m==1)&& (l==0) &&
          (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1))
        {
          if (TDatabase::ParamDB->NSTYPE <4)
          {
            OutPut("For slip with friction bc NSTYPE 4 is ");
            OutPut("necessary !!!!! " << endl);
            exit(4711);
          }

          // prepare everything for the assembling of slip with friction bc
          // on all levels
          N_FESpaces = 1;
          N_SquareMatrices = 8;
          N_RectMatrices = 2;
          N_Rhs = 2;
          DiscreteForm = NULL;

          for(i=0;i<mg_level;i++)
          {
            SQMATRICES[0] = MatricesA11[i];
            SQMATRICES[1] = MatricesA22[i];
            SQMATRICES[2] = MatricesA12[i];
            SQMATRICES[3] = MatricesA21[i];
            SQMATRICES[4] = MatricesM11[i];
            SQMATRICES[5] = MatricesM22[i];
            SQMATRICES[6] = MatricesM12[i];
            SQMATRICES[7] = MatricesM21[i];

            MATRICES[0] = MatricesB1T[i];
            MATRICES[1] = MatricesB2T[i];

            fesp[0] = USpaces[i];
            ferhs[0] = USpaces[i];
            ferhs[1] = USpaces[i];

            RHSs[0] = RhsArray[i];
            RHSs[1] = RhsArray[i]+N_Uarray[i];

            Assemble2DSlipBC(N_FESpaces, fesp,
              N_SquareMatrices, SQMATRICES,
              N_RectMatrices, MATRICES,
              N_Rhs, RHSs, ferhs,
              DiscreteForm,
              BoundaryConditions,
              BoundValues,
              aux,
              U1Array[i],U2Array[i]);
            TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
          }
        }
        delete aux;

        //======================================================================
        // manipulation of matrices due to current time discretization
        // the stiffness matrix is stored on M11, (M12, M21, M22)
        //======================================================================

        // scale B1T, B2T, B1, B2
	// in the very first step, oldtau = 1 and the matrices are scaled by tau
        for(i=0;i<mg_level;i++)
        {
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
            case 3:
		if (tau/oldtau != 1.0)
		{
		    Dscal(MatricesB1[i]->GetN_Entries(),
			  press_factor1*tau/oldtau,
			  MatricesB1[i]->GetEntries());
		    Dscal(MatricesB2[i]->GetN_Entries(),
			  press_factor1*tau/oldtau,
			  MatricesB2[i]->GetEntries());
		}
              break;

            case 2:
            case 4:
		if (tau/oldtau != 1.0)
		{
		    Dscal(MatricesB1T[i]->GetN_Entries(),
			  press_factor1*tau/oldtau,
			  MatricesB1T[i]->GetEntries());
		    Dscal(MatricesB2T[i]->GetN_Entries(),
			  press_factor1*tau/oldtau,
			  MatricesB2T[i]->GetEntries());
		}
		// scale divergence constraint
		if ((TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT>0)
		    &&(TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT*press_factor1*tau/oldtau != 1.0))
		{
		    Dscal(MatricesB1[i]->GetN_Entries(),
			  TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT*press_factor1*tau/oldtau,
			  MatricesB1[i]->GetEntries());
		    Dscal(MatricesB2[i]->GetN_Entries(),
			  TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT*press_factor1*tau/oldtau,
			  MatricesB2[i]->GetEntries());
		}
              break;
          }
        }                        // endfor

        // update rhs by Laplacian and convective term from previous
        // time step
        // scaled by current sub time step length and theta2
        // currently : M := M + gamma A
        // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2) A

        for(i=0;i<mg_level;i++)
        {
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
            case 2:
              MatAdd(MatricesM[i], MatricesA[i],
                -gamma - tau*TDatabase::TimeDB->THETA2);
              break;

            case 3:
            case 4:
              MatAdd(MatricesM11[i], MatricesA11[i],
                -gamma - tau*TDatabase::TimeDB->THETA2);
              MatAdd(MatricesM12[i], MatricesA12[i],
                -gamma - tau*TDatabase::TimeDB->THETA2);
              MatAdd(MatricesM21[i], MatricesA21[i],
                -gamma - tau*TDatabase::TimeDB->THETA2);
              MatAdd(MatricesM22[i], MatricesA22[i],
                -gamma - tau*TDatabase::TimeDB->THETA2);
              break;
          }                      // endswitch
        }
        // set current factor of steady state matrix
        gamma = -tau*TDatabase::TimeDB->THETA2;

        // defect = M * sol
        // B:= B + defect
        switch(TDatabase::ParamDB->NSTYPE)
        {
          case 1:
          case 2:
            MatVectActive(MatricesM[mg_level-1], sol, defect);
            MatVectActive(MatricesM[mg_level-1], sol+N_U, defect+N_U);
            Daxpy(N_Active, 1, defect, B);
            Daxpy(N_Active, 1, defect+N_U, B+N_U);
            break;

          case 3:
          case 4:
            MatVectActive(MatricesM11[mg_level-1], sol, defect);
            Daxpy(N_Active, 1, defect, B);
            MatVectActive(MatricesM12[mg_level-1], sol+N_U, defect);
            Daxpy(N_Active, 1, defect, B);
            MatVectActive(MatricesM21[mg_level-1], sol, defect+N_U);
            Daxpy(N_Active, 1, defect+N_U, B+N_U);
            MatVectActive(MatricesM22[mg_level-1], sol+N_U, defect+N_U);
            Daxpy(N_Active, 1, defect+N_U, B+N_U);
            break;
        }

	// treating the pressure in the same way as the velocity
	// with respect to the temporal discretization 
        if (TDatabase::TimeDB->EXTRAPOLATE_PRESSURE==1)
        {
	    // the coupling matrices are scaled by press_factor1*tau 
	    // at this moment
	    val0 = -press_factor2/press_factor1;
	    switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
            case 3:
		/*	OutPut("val 0 " << val0 << " " << Ddot(N_P,sol+2*N_U,sol+2*N_U) << 
		       " " <<  Ddot(MatricesB1[mg_level-1]->GetN_Entries(), 
				    MatricesB1[mg_level-1]->GetEntries(), 
				    MatricesB1[mg_level-1]->GetEntries()) << 
		       " " <<  Ddot(MatricesB2[mg_level-1]->GetN_Entries(), 
				    MatricesB2[mg_level-1]->GetEntries(), 
				    MatricesB2[mg_level-1]->GetEntries()) << 
				    endl);*/
		TransMatVect(MatricesB1[mg_level-1], sol+2*N_U, defect);
		TransMatVect(MatricesB2[mg_level-1], sol+2*N_U, defect+N_U);
		Daxpy(N_Active, val0, defect, B);
		Daxpy(N_Active, val0, defect+N_U, B+N_U);
		// OutPut("def1 " << Ddot(N_Active,defect,defect) <<  " " << Ddot(N_Active,B,B) << endl);
		// OutPut("def2 " << Ddot(N_Active,defect+N_U,defect+N_U) <<  " " << 
		// Ddot(N_Active,B+N_U,B+N_U) << endl);
             break;

            case 2:
            case 4:
              MatVect1(MatricesB1T[mg_level-1], sol+2*N_U, defect);
              MatVect1(MatricesB2T[mg_level-1], sol+2*N_U, defect+N_U);
              Daxpy(N_Active, val0, defect, B);
              Daxpy(N_Active, val0, defect+N_U, B+N_U);
              break;
          }
        }

        // extrapolate solution to get starting value at next time step
        // current solution sol
        // solution of last time step sol_timestep_m1

        // save sol
        memcpy(oldsol, sol, SizeOfDouble*N_Unknowns);
        // TDatabase::TimeDB->EXTRAPOLATE_WEIGHT = 1 : linear extrapolation
        // TDatabase::TimeDB->EXTRAPOLATE_WEIGHT = 0 : constant extrapolation
        //         (take solution of previous discrete time as start solution)
        tau2 = TDatabase::TimeDB->EXTRAPOLATE_WEIGHT*tau/oldtau;
        tau1 = 1 + tau2;
        // sol := tau1 *sol - tau2 * sol_timestep_m1
        // at first time step: sol = sol_timestep_m1 -> result is sol
        //for (k=0;k<2*N_U;k++)
        //  sol[k] = tau1*sol[k] - tau2*sol_timestep_m1[k];
        oldtau = tau;
        // save current solution
        memcpy(sol_timestep_m1, oldsol, SizeOfDouble*N_Unknowns);

        // set Dirichlet values
        // RHSs[0] still available from assembling
        memcpy(B+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
        memcpy(B+N_Active+N_U, RHSs[1]+N_Active, (N_U-N_Active)*SizeOfDouble);

        // copy Dirichlet values from rhs into sol
        memcpy(sol+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
        memcpy(sol+N_Active+N_U, RHSs[1]+N_Active, (N_U-N_Active)*SizeOfDouble);

        //========================================================================
        // end assembling of rhs
        //========================================================================

        // M = M + (-gamma ) A
        for(i=0;i<mg_level;i++)
        {
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
            case 2:
              MatAdd(MatricesM[i], MatricesA[i],
                -gamma);
              break;

            case 3:
            case 4:
              MatAdd(MatricesM11[i], MatricesA11[i],
                -gamma);
              MatAdd(MatricesM12[i], MatricesA12[i],
                -gamma);
              MatAdd(MatricesM21[i], MatricesA21[i],
                -gamma);
              MatAdd(MatricesM22[i], MatricesA22[i],
                -gamma);
              break;
          }                      // endswitch
        }
        gamma = 0;

        if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
          MG->RestrictToAllGrids();

        //======================================================================
        // assemble new matrix due to extrapolation of solution
        //======================================================================
        // for all levels
        for(i=0;i<mg_level;i++)
        {
          if ((mg_type==1) && (i<mg_level-1))
          {
            DiscreteForm = DiscreteFormNLUpwind;
            CurrentDiscType =  UPWIND;
          }
          else
            switch(TDatabase::ParamDB->DISCTYPE)
            {
              case GALERKIN:
              case SMAGORINSKY_EXPL:
                DiscreteForm = DiscreteFormNLGalerkin;
                CurrentDiscType =  GALERKIN;
                break;

              case UPWIND:
                DiscreteForm = DiscreteFormNLUpwind;
                CurrentDiscType =  UPWIND;
                break;

              case SMAGORINSKY:
                DiscreteForm = DiscreteFormNLSmagorinsky;
                CurrentDiscType =  SMAGORINSKY;
                break;

              case CLASSICAL_LES:
                DiscreteForm = DiscreteFormNLColetti;
                CurrentDiscType =  CLASSICAL_LES ;
                break;

              case  GL00_CONVOLUTION:
                DiscreteForm = DiscreteFormNLGL00Convolution;
                CurrentDiscType =  GL00_CONVOLUTION;
                break;

              case  GL00_AUX_PROBLEM:
                DiscreteForm = DiscreteFormNLGL00AuxProblem;
                CurrentDiscType =  GL00_AUX_PROBLEM;
                break;

              case VMS_PROJECTION:
                DiscreteForm = DiscreteFormNLVMSProjection;
                CurrentDiscType =  VMS_PROJECTION;
                break;

              default:
                OutPut("Unknown DISCTYPE" << endl);
                exit(1);
            }
          // set pointers to matrices
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
            case 2:
              SQMATRICES[0] = MatricesA[i];
              SQMATRICES[0]->Reset();

              N_SquareMatrices = 1;
              N_RectMatrices = 0;
              N_Rhs = 0;
              N_FESpaces = 1;
              break;

            case 3:
            case 4:
              N_RectMatrices = 0;

              N_Rhs = 0;
              N_FESpaces = 1;

              if (((TDatabase::ParamDB->LAPLACETYPE==1)||
		   (TDatabase::ParamDB->NSE_NONLINEAR_FORM==2))
                &&
                ((CurrentDiscType == SMAGORINSKY) ||
                (CurrentDiscType == VMS_PROJECTION)||
                (CurrentDiscType == CLASSICAL_LES) ||
                (CurrentDiscType == GL00_CONVOLUTION) ||
                (CurrentDiscType == GL00_AUX_PROBLEM)))
              {
                SQMATRICES[0] = MatricesA11[i];
                SQMATRICES[1] = MatricesA12[i];
                SQMATRICES[2] = MatricesA21[i];
                SQMATRICES[3] = MatricesA22[i];
                SQMATRICES[0]->Reset();
                SQMATRICES[1]->Reset();
                SQMATRICES[2]->Reset();
                SQMATRICES[3]->Reset();
                N_SquareMatrices = 4;
                last_sq = 3;
                if (CurrentDiscType == VMS_PROJECTION)
                {
                  N_RectMatrices = 2;
                  MATRICES[0] = Matrices_tilde_G11[i];
                  MATRICES[1] = Matrices_tilde_G22[i];
                  MATRICES[0]->Reset();
                  MATRICES[1]->Reset();
                  N_FESpaces = 4;
                  fesp[3] = ProjectionSpaces[i];
                }
              }
              else
              {
                SQMATRICES[0] = MatricesA11[i];
                SQMATRICES[1] = MatricesA22[i];

                SQMATRICES[0]->Reset();
                SQMATRICES[1]->Reset();

                N_SquareMatrices = 2;
                last_sq = 1;
              }
              break;
          }

          fesp[0] = USpaces[i];

          fefct[0] = U1Array[i];
          fefct[1] = U2Array[i];

          ferhs[0] = USpaces[i];
          ferhs[1] = USpaces[i];

          switch(TDatabase::ParamDB->DISCTYPE)
          {
            // turbulent viscosity must be computed
            case SMAGORINSKY:
            case VMS_PROJECTION:
            case CLASSICAL_LES:
            case GL00_CONVOLUTION:
            case GL00_AUX_PROBLEM:
              // standard turbulent viscosities on coarser grids
              if ((TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE!=4)
                || (i<mg_level -1))
              {
                aux =  new TAuxParam2D(TimeNSN_FESpacesVelo_GradVelo,
                  TimeNSN_FctVelo_GradVelo,
                  TimeNSN_ParamFctVelo_GradVelo,
                  TimeNSN_FEValuesVelo_GradVelo,
                  fesp, fefct,
                  TimeNSFctVelo_GradVelo,
                  TimeNSFEFctIndexVelo_GradVelo,
                  TimeNSFEMultiIndexVelo_GradVelo,
                  TimeNSN_ParamsVelo_GradVelo,
                  TimeNSBeginParamVelo_GradVelo);
                if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 4)
                  TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE=17;
              }
              // u^n - g_\delta * u^0 on finest level
              if ((TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 4)
                && (i==mg_level -1))
              {
                fesp[1] = uConvSpaces[mg_level-1];
                N_FESpaces=2;
                fefct[2] = u1ConvArray[mg_level-1];
                fefct[3] = u2ConvArray[mg_level-1];
                //OutPut("aa " << u1ConvArray[mg_level-1]->GetValues()[0] << endl);
                // OutPut("aa " << (int) fefct[0] << endl);
                // OutPut("aa " << (int) fefct[1] << endl);
                //OutPut("aa " << (int) fefct[2] << endl);
                //OutPut("aa " << (int) fefct[3] << endl);

                aux =  new TAuxParam2D(TimeNSN_FESpacesVelo_GradVelo_ConvVelo,
                  TimeNSN_FctVelo_GradVelo_ConvVelo,
                  TimeNSN_ParamFctVelo_GradVelo_ConvVelo,
                  TimeNSN_FEValuesVelo_GradVelo_ConvVelo,
                  fesp, fefct,
                  TimeNSFctVelo_GradVelo_ConvVelo,
                  TimeNSFEFctIndexVelo_GradVelo_ConvVelo,
                  TimeNSFEMultiIndexVelo_GradVelo_ConvVelo,
                  TimeNSN_ParamsVelo_GradVelo_ConvVelo,
                  TimeNSBeginParamVelo_GradVelo_ConvVelo);
              }
              break;
            default:
              aux =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
                TimeNSN_ParamFct2,
                TimeNSN_FEValues2,
                fesp, fefct,
                TimeNSFct2,
                TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
                TimeNSN_Params2, TimeNSBeginParam2);
              break;
          }

          //======================================================================
          // assembling of matrices for each level due to nonlinearity
          // A_11, (A_22)
          // no assembling of rhs
          //======================================================================

          Assemble2D(N_FESpaces, fesp,
            N_SquareMatrices, SQMATRICES,
            N_RectMatrices, MATRICES,
            N_Rhs, RHSs, ferhs,
            DiscreteForm,
            BoundaryConditions,
            BoundValues,
            aux);

          // reset parameter
          if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 17)
            TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE=4;

          if(DiscreteForm == DiscreteFormNLUpwind)
          {
            switch(TDatabase::ParamDB->NSTYPE)
            {
              case 1:
              case 2:
                // do upwinding with one matrix
                UpwindForNavierStokes(LinCoeffs, SQMATRICES[0], U1Array[i], U2Array[i]);
                //cout << "UPWINDING DONE : level " << i << endl;
                break;

              case 3:
              case 4:
                // do upwinding with two matrices
                UpwindForNavierStokes(LinCoeffs, SQMATRICES[0], U1Array[i], U2Array[i]);
                UpwindForNavierStokes(LinCoeffs, SQMATRICES[last_sq], U1Array[i], U2Array[i]);
                //cout << "UPWINDING DONE(2) : level " << i << endl;
                //cout << "check correct sqmatrix !!!! " << endl;
                break;
            }                    // endswitch
          }

          // slip type bc detected, modify matrices accordingly
          if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
          {

            // prepare everything for the assembling of slip with friction bc
            // on level i
            N_FESpaces = 1;
            N_SquareMatrices = 2;
            N_RectMatrices = 0;
            N_Rhs = 2;
            DiscreteForm = NULL;

            SQMATRICES[0] = MatricesA11[i];
            SQMATRICES[1] = MatricesA22[i];

            fesp[0] = USpaces[i];
            ferhs[0] = USpaces[i];
            ferhs[1] = USpaces[i];

            RHSs[0] = RhsArray[i];
            RHSs[1] = RhsArray[i]+N_Uarray[i];

            Assemble2DSlipBC(N_FESpaces, fesp,
              N_SquareMatrices, SQMATRICES,
              N_RectMatrices, MATRICES,
              N_Rhs, RHSs, ferhs,
              DiscreteForm,
              BoundaryConditions,
              BoundValues,
              aux,
              U1Array[i],U2Array[i]);
          }

          // update matrices
          if (CurrentDiscType == VMS_PROJECTION)
          {
            SQMATRICES[0] = MatricesA11[i];
            SQMATRICES[1] = MatricesA12[i];
            SQMATRICES[2] = MatricesA21[i];
            SQMATRICES[3] = MatricesA22[i];
            SQMATRICES[6] =  MatricesL[i];
            MATRICES[2] = Matrices_tilde_G11[i];
            MATRICES[3] = Matrices_tilde_G22[i];
            MATRICES[4] = Matrices_G11[i];
            MATRICES[5] = Matrices_G22[i];

            VMSProjectionUpdateMatrices(N_Uarray[i], USpaces[i]->GetActiveBound(),
              ProjectionSpaces[i]->GetN_DegreesOfFreedom(),
              SQMATRICES,MATRICES);
          }

          delete aux;
          //======================================================================
          // end of assemble new matrix due to nonlinearity
          //======================================================================

        }                        // endfor i
        //****************
        //========================================================================
        // assembling of system matrix
        //========================================================================

        // M = M + (-gamma + tau*TDatabase::TimeDB->THETA1) A
        for(i=0;i<mg_level;i++)
        {
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
            case 2:
              MatAdd(MatricesM[i], MatricesA[i],
                -gamma + tau*TDatabase::TimeDB->THETA1);
              break;

            case 3:
            case 4:
              MatAdd(MatricesM11[i], MatricesA11[i],
                -gamma + tau*TDatabase::TimeDB->THETA1);
              MatAdd(MatricesM12[i], MatricesA12[i],
                -gamma + tau*TDatabase::TimeDB->THETA1);
              MatAdd(MatricesM21[i], MatricesA21[i],
                -gamma + tau*TDatabase::TimeDB->THETA1);
              MatAdd(MatricesM22[i], MatricesA22[i],
                -gamma + tau*TDatabase::TimeDB->THETA1);
              break;
          }                      // endswitch
        }
        // set current factor of steady state matrix
        gamma = tau*TDatabase::TimeDB->THETA1;

        //========================================================================
        // end assembling of system matrix
        //========================================================================

        OutPut(endl << "CURRENT TIME: ");
        OutPut(TDatabase::TimeDB->CURRENTTIME << endl);
        OutPut(TDatabase::TimeDB->CURRENTTIME << " memory: " << setw(10) << GetMemory() << endl);

	if ((hmin*tau<=1e-5)&&((TDatabase::ParamDB->NSTYPE==1)||(TDatabase::ParamDB->NSTYPE==3)))
	{
	    OutPut("ATTENTION: divergence constraint will be multiplied with small time step !!!"<<endl);
	    OutPut("remove the exit if this is correct " << endl);
	    //exit(4711);
	}

        //======================================================================
        // nonlinear loop
        //======================================================================
        N_LinIterCurr = 0;
        solver_time_curr = 0;

        for(j=0;j<=Max_It;j++)    // solve nonlinear equation
        {
          memcpy(oldsol, sol, SizeOfDouble*N_Unknowns);
          memset(defect, 0, N_Unknowns*SizeOfDouble);

          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
              SQMATRICES[0] = MatricesM[mg_level-1];
              MATRICES[0] = MatricesB1[mg_level-1];
              MATRICES[1] = MatricesB2[mg_level-1];
              break;
            case 2:
              SQMATRICES[0] = MatricesM[mg_level-1];
              MATRICES[0] = MatricesB1[mg_level-1];
              MATRICES[1] = MatricesB2[mg_level-1];
              MATRICES[2] = MatricesB1T[mg_level-1];
              MATRICES[3] = MatricesB2T[mg_level-1];
              break;
            case 3:
              SQMATRICES[0] = MatricesM11[mg_level-1];
              SQMATRICES[1] = MatricesM12[mg_level-1];
              SQMATRICES[2] = MatricesM21[mg_level-1];
              SQMATRICES[3] = MatricesM22[mg_level-1];
              MATRICES[0] = MatricesB1[mg_level-1];
              MATRICES[1] = MatricesB2[mg_level-1];
              break;
            case 4:
              SQMATRICES[0] = MatricesM11[mg_level-1];
              SQMATRICES[1] = MatricesM12[mg_level-1];
              SQMATRICES[2] = MatricesM21[mg_level-1];
              SQMATRICES[3] = MatricesM22[mg_level-1];
              MATRICES[0] = MatricesB1[mg_level-1];
              MATRICES[1] = MatricesB2[mg_level-1];
              MATRICES[2] = MatricesB1T[mg_level-1];
              MATRICES[3] = MatricesB2T[mg_level-1];
              break;
          }
          // compute defect
          Defect(sqmatrices,matrices,sol,B,defect);
	  // scale divergence constraint appropriately
	  /* if ((TDatabase::ParamDB->NSTYPE==1) || (TDatabase::ParamDB->NSTYPE==3))
	  {
	      Dscal(N_P, 1.0/tau, defect+2*N_U);
	      TDatabase::ParamDB->INTERNAL_DIV_SCALE = 1.0/tau;
	      }*/
          if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
            IntoL20Vector2D(defect+2*N_U, N_P,pressure_space_code);
          residual =  Ddot(N_Unknowns, defect, defect);
          impuls_residual = Ddot(2*N_U, defect, defect);
          OutPut("nonlinear step " << setw(3) << j);
          OutPut(setw(14) << impuls_residual);
          OutPut(setw(14) << Ddot(N_P,defect+2*N_U,defect+2*N_U));
          OutPut(setw(14) << sqrt(residual));
          if (j>0)
          {
            OutPut(setw(14) << sqrt(residual)/oldresidual << endl);
          }
          else
          {
            OutPut(endl);
          }
          oldresidual = sqrt(residual);

	  if ((isnan(residual))||(isinf(residual)))
	  {
	      OutPut("Iteration diverged !!!" << endl);
	      exit(4711);
	  }

          if ((((sqrt(residual)<=limit)||(j==Max_It)))
            && (j>=TDatabase::ParamDB->SC_MINIT))
          {
            if (j==Max_It-1)
              j++;
            OutPut("ITE : " << setw(3) << j);
            OutPut(" (" << setw(3) << N_LinIterCurr << "/");
            OutPut(setw(3) << N_LinIter << " LINITE)");
            OutPut("  TIME FOR SOLVER : " << solver_time_curr << "/" << solver_time << "s");
            OutPut("  RES : " <<  sqrt(residual) << endl);
            // count total running time
            t4 =  GetTime();
            total_time += t4 - t3;
            t3 = t4;
            OutPut(TDatabase::TimeDB->CURRENTTIME << " total elapsed running time " <<    total_time << endl);
            break;
          }
 
// 		 		  cout << " Solve Direct  " << TDatabase::ParamDB->SOLVER_TYPE << endl;
// 				  exit(0);
				  
          //======================================================================
          // solve linear system
          //======================================================================
          switch(TDatabase::ParamDB->SOLVER_TYPE)
          {
            case AMG:
              TDatabase::ParamDB->SC_VERBOSE=0;
              t1 = GetTime();
              switch(TDatabase::ParamDB->NSTYPE)
              {
                case 1:
                  Solver(sqmatrixM, matrixB1, matrixB2, B, sol);
                  break;

                case 2:
                  Solver(sqmatrixM, matrixB1T, matrixB2T,
                    matrixB1, matrixB2, B, sol);
                  break;

                case 3:
                  Error("AMG does not work for NSTYPE = 3." << endl);
                  return -1;
                  break;

                case 4:

                  DirectSolver(sqmatrixM11, sqmatrixM12, sqmatrixM21, sqmatrixM22,
                    matrixB1T, matrixB2T, matrixB1, matrixB2, B, sol, 3); 
                  break;
              }
              t2 = GetTime();
              solver_time_curr = t2-t1;
              solver_time += solver_time_curr;
              break;

            case GMG:
              t1 = GetTime();

	      if (umfpack_flag==-2)
	      {
		  MG->SetParam(9,0);
		  umfpack_flag=-1;
	      }
	      else
		  MG->SetParam(9,-1);
	  		  
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
              {
                memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
                memcpy(itmethod_rhs, B, N_Unknowns*SizeOfDouble);
              }
              N_LinIterCurrIte = itmethod->Iterate(sqmatrices,matrices,itmethod_sol,itmethod_rhs);
              N_LinIterCurr += N_LinIterCurrIte;
              N_LinIter += N_LinIterCurrIte;
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE == 5)
              {
                memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
                memcpy(B, itmethod_rhs, N_Unknowns*SizeOfDouble);
              }
              t2 = GetTime();
              solver_time_curr += t2-t1;
              solver_time += solver_time_curr;

              // p1 = 0;
              for(k=0;k<N_Unknowns;k++)
              {
                p2 = sol[k]-oldsol[k];
                sol[k] = oldsol[k] + omega * p2;
                // p1 += p2*p2;
              }
              break;
          }                      // endswitch SOLVER_TYPE
          //======================================================================
          // end solve linear system
          //======================================================================

          // restore mass matrices by subtracting the A-matrices
          for(i=0;i<mg_level;i++)
          {
            switch(TDatabase::ParamDB->NSTYPE)
            {
              case 1:
              case 2:
                MatAdd(MatricesM[i], MatricesA[i], -gamma);
                break;

              case 3:
              case 4:
                MatAdd(MatricesM11[i], MatricesA11[i], -gamma);
                MatAdd(MatricesM12[i], MatricesA12[i], -gamma);
                MatAdd(MatricesM21[i], MatricesA21[i], -gamma);
                MatAdd(MatricesM22[i], MatricesA22[i], -gamma);
                break;
            }                    // endswitch
          }                      // endfor i
          // set current factor of steady state matrix
          gamma = 0;

          if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
            MG->RestrictToAllGrids();

          //======================================================================
          // assemble new matrix due to nonlinearity
          //======================================================================
          // for all levels
          for(i=0;i<mg_level;i++)
          {
            if ((mg_type==1) && (i<mg_level-1))
            {
              DiscreteForm = DiscreteFormNLUpwind;
              CurrentDiscType =  UPWIND;
            }
            else
              switch(TDatabase::ParamDB->DISCTYPE)
              {
                case GALERKIN:
                case SMAGORINSKY_EXPL:
                  DiscreteForm = DiscreteFormNLGalerkin;
                  CurrentDiscType =  GALERKIN;
                  break;

                case UPWIND:
                  DiscreteForm = DiscreteFormNLUpwind;
                  CurrentDiscType =  UPWIND;
                  break;

                case SMAGORINSKY:
                  DiscreteForm = DiscreteFormNLSmagorinsky;
                  CurrentDiscType =  SMAGORINSKY;
                  break;

                case  CLASSICAL_LES:
                  DiscreteForm = DiscreteFormNLColetti;
                  CurrentDiscType =  CLASSICAL_LES ;
                  break;

                case  GL00_CONVOLUTION:
                  DiscreteForm = DiscreteFormNLGL00Convolution;
                  CurrentDiscType =  GL00_CONVOLUTION;
                  break;

                case  GL00_AUX_PROBLEM:
                  DiscreteForm = DiscreteFormNLGL00AuxProblem;
                  CurrentDiscType =  GL00_AUX_PROBLEM;
                  break;

                case VMS_PROJECTION:
                  DiscreteForm = DiscreteFormNLVMSProjection;
                  CurrentDiscType =  VMS_PROJECTION;
                  break;

                default:
                  OutPut("Unknown DISCTYPE" << endl);
                  exit(1);
              }
            // set pointers to matrices
            switch(TDatabase::ParamDB->NSTYPE)
            {
              case 1:
              case 2:
                SQMATRICES[0] = MatricesA[i];
                SQMATRICES[0]->Reset();

                N_SquareMatrices = 1;
                N_RectMatrices = 0;
                N_Rhs = 0;
                N_FESpaces = 1;
                break;

              case 3:
              case 4:
                N_RectMatrices = 0;

                N_Rhs = 0;
                N_FESpaces = 1;

              if (((TDatabase::ParamDB->LAPLACETYPE==1)||
		   (TDatabase::ParamDB->NSE_NONLINEAR_FORM==2))
		  &&
                  ((CurrentDiscType == SMAGORINSKY) ||
                  (CurrentDiscType == VMS_PROJECTION) ||
                  (CurrentDiscType == CLASSICAL_LES) ||
                  (CurrentDiscType == GL00_CONVOLUTION) ||
                  (CurrentDiscType == GL00_AUX_PROBLEM)))
                {
                  SQMATRICES[0] = MatricesA11[i];
                  SQMATRICES[1] = MatricesA12[i];
                  SQMATRICES[2] = MatricesA21[i];
                  SQMATRICES[3] = MatricesA22[i];
                  SQMATRICES[0]->Reset();
                  SQMATRICES[1]->Reset();
                  SQMATRICES[2]->Reset();
                  SQMATRICES[3]->Reset();
                  N_SquareMatrices = 4;
                  last_sq = 3;
                  if (CurrentDiscType == VMS_PROJECTION)
                  {
                    N_RectMatrices = 2;
                    MATRICES[0] = Matrices_tilde_G11[i];
                    MATRICES[1] = Matrices_tilde_G22[i];
                    MATRICES[0]->Reset();
                    MATRICES[1]->Reset();
                    N_FESpaces = 4;
                    fesp[3] = ProjectionSpaces[i];
                  }
                }
                else
                {
                  SQMATRICES[0] = MatricesA11[i];
                  SQMATRICES[1] = MatricesA22[i];

                  SQMATRICES[0]->Reset();
                  SQMATRICES[1]->Reset();

                  N_SquareMatrices = 2;
                  last_sq = 1;
                }
                break;
            }

            fesp[0] = USpaces[i];

            fefct[0] = U1Array[i];
            fefct[1] = U2Array[i];

            ferhs[0] = USpaces[i];
            ferhs[1] = USpaces[i];

            switch(TDatabase::ParamDB->DISCTYPE)
            {
              // turbulent viscosity must be computed
              case SMAGORINSKY:
              case VMS_PROJECTION:
              case CLASSICAL_LES:
              case GL00_CONVOLUTION:
              case GL00_AUX_PROBLEM:
                // standard turbulent viscosities on coarser grids
                if ((TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE!=4)
                  || (i<mg_level -1))
                {
                  aux =  new TAuxParam2D(TimeNSN_FESpacesVelo_GradVelo,
                    TimeNSN_FctVelo_GradVelo,
                    TimeNSN_ParamFctVelo_GradVelo,
                    TimeNSN_FEValuesVelo_GradVelo,
                    fesp, fefct,
                    TimeNSFctVelo_GradVelo,
                    TimeNSFEFctIndexVelo_GradVelo,
                    TimeNSFEMultiIndexVelo_GradVelo,
                    TimeNSN_ParamsVelo_GradVelo,
                    TimeNSBeginParamVelo_GradVelo);
                  if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 4)
                    TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE=17;
                }
                // u^n - g_\delta * u^0 on finest level
                if ((TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 4)
                  && (i==mg_level -1))
                {
                  fesp[1] = uConvSpaces[mg_level-1];
                  N_FESpaces=2;
                  fefct[2] = u1ConvArray[mg_level-1];
                  fefct[3] = u2ConvArray[mg_level-1];
                  //OutPut("aa " << u1ConvArray[mg_level-1]->GetValues()[0] << endl);
                  // OutPut("aa " << (int) fefct[0] << endl);
                  // OutPut("aa " << (int) fefct[1] << endl);
                  //OutPut("aa " << (int) fefct[2] << endl);
                  //OutPut("aa " << (int) fefct[3] << endl);

                  aux =  new TAuxParam2D(TimeNSN_FESpacesVelo_GradVelo_ConvVelo,
                    TimeNSN_FctVelo_GradVelo_ConvVelo,
                    TimeNSN_ParamFctVelo_GradVelo_ConvVelo,
                    TimeNSN_FEValuesVelo_GradVelo_ConvVelo,
                    fesp, fefct,
                    TimeNSFctVelo_GradVelo_ConvVelo,
                    TimeNSFEFctIndexVelo_GradVelo_ConvVelo,
                    TimeNSFEMultiIndexVelo_GradVelo_ConvVelo,
                    TimeNSN_ParamsVelo_GradVelo_ConvVelo,
                    TimeNSBeginParamVelo_GradVelo_ConvVelo);
                }
                break;
              default:
                aux =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
                  TimeNSN_ParamFct2,
                  TimeNSN_FEValues2,
                  fesp, fefct,
                  TimeNSFct2,
                  TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
                  TimeNSN_Params2, TimeNSBeginParam2);
                break;
            }

            //======================================================================
            // assembling of matrices for each level due to nonlinearity
            // A_11, (A_22)
            // no assembling of rhs
            //======================================================================

            Assemble2D(N_FESpaces, fesp,
              N_SquareMatrices, SQMATRICES,
              N_RectMatrices, MATRICES,
              N_Rhs, RHSs, ferhs,
              DiscreteForm,
              BoundaryConditions,
              BoundValues,
              aux);

            // reset parameter
            if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 17)
              TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE=4;

            if(DiscreteForm == DiscreteFormNLUpwind)
            {
              switch(TDatabase::ParamDB->NSTYPE)
              {
                case 1:
                case 2:
                  // do upwinding with one matrix
                  UpwindForNavierStokes(LinCoeffs, SQMATRICES[0], U1Array[i], U2Array[i]);
                  //cout << "UPWINDING DONE : level " << i << endl;
                  break;

                case 3:
                case 4:
                  // do upwinding with two matrices
                  UpwindForNavierStokes(LinCoeffs, SQMATRICES[0], U1Array[i], U2Array[i]);
                  UpwindForNavierStokes(LinCoeffs, SQMATRICES[last_sq], U1Array[i], U2Array[i]);
                  //cout << "UPWINDING DONE(2) : level " << i << endl;
                  //cout << "check correct sqmatrix !!!! " << endl;
                  break;
              }                  // endswitch
            }

            // slip type bc detected, modify matrices accordingly
            if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
            {

              // prepare everything for the assembling of slip with friction bc
              // on level i
              N_FESpaces = 1;
              N_SquareMatrices = 2;
              N_RectMatrices = 0;
              N_Rhs = 2;
              DiscreteForm = NULL;

              SQMATRICES[0] = MatricesA11[i];
              SQMATRICES[1] = MatricesA22[i];

              fesp[0] = USpaces[i];
              ferhs[0] = USpaces[i];
              ferhs[1] = USpaces[i];

              RHSs[0] = RhsArray[i];
              RHSs[1] = RhsArray[i]+N_Uarray[i];

              Assemble2DSlipBC(N_FESpaces, fesp,
                N_SquareMatrices, SQMATRICES,
                N_RectMatrices, MATRICES,
                N_Rhs, RHSs, ferhs,
                DiscreteForm,
                BoundaryConditions,
                BoundValues,
                aux,
                U1Array[i],U2Array[i]);
            }

            // update matrices
            if (CurrentDiscType == VMS_PROJECTION)
            {
              SQMATRICES[0] = MatricesA11[i];
              SQMATRICES[1] = MatricesA12[i];
              SQMATRICES[2] = MatricesA21[i];
              SQMATRICES[3] = MatricesA22[i];
              SQMATRICES[6] =  MatricesL[i];
              MATRICES[2] = Matrices_tilde_G11[i];
              MATRICES[3] = Matrices_tilde_G22[i];
              MATRICES[4] = Matrices_G11[i];
              MATRICES[5] = Matrices_G22[i];

              VMSProjectionUpdateMatrices(N_Uarray[i], USpaces[i]->GetActiveBound(),
                ProjectionSpaces[i]->GetN_DegreesOfFreedom(),
                SQMATRICES,MATRICES);
            }

            delete aux;
            //======================================================================
            // end of assemble new matrix due to nonlinearity
            //======================================================================

            // build stiffness matrix for next nonlinear iteration step
            // stiffness matrix (left upper block) is stored on
            // M11, (M12, M21, M22)
            // M = M +  tau*TDatabase::TimeDB->THETA1 A
            switch(TDatabase::ParamDB->NSTYPE)
            {
              case 1:
              case 2:
                MatAdd(MatricesM[i], MatricesA[i],
                  tau*TDatabase::TimeDB->THETA1);
                break;

              case 3:
              case 4:
                MatAdd(MatricesM11[i], MatricesA11[i],
                  tau*TDatabase::TimeDB->THETA1);
                MatAdd(MatricesM12[i], MatricesA12[i],
                  tau*TDatabase::TimeDB->THETA1);
                MatAdd(MatricesM21[i], MatricesA21[i],
                  tau*TDatabase::TimeDB->THETA1);
                MatAdd(MatricesM22[i], MatricesA22[i],
                  tau*TDatabase::TimeDB->THETA1);
                break;
            }
          }                      // endfor i
          // set current factor of steady state matrix
          gamma = tau*TDatabase::TimeDB->THETA1;

        }                        // endfor Max_It (solution of nonlinear equation)

        //======================================================================
        // end of nonlinear loop
        //======================================================================

        IntoL20FEFunction(PArray[mg_level-1]->GetValues(), N_Parray[mg_level-1],
          PSpaces[mg_level-1], velocity_space_code, pressure_space_code);

	// reset scaling of coupling matrices
        if (TDatabase::TimeDB->EXTRAPOLATE_PRESSURE==1)
        {
	  for(i=0;i<mg_level;i++)
          {
            switch(TDatabase::ParamDB->NSTYPE)
            {
              case 1:
              case 3:
                Dscal(MatricesB1[i]->GetN_Entries(),
                  1.0/press_factor1,
                  MatricesB1[i]->GetEntries());
                Dscal(MatricesB2[i]->GetN_Entries(),
                  1.0/press_factor1,
                  MatricesB2[i]->GetEntries());
                break;

              case 2:
              case 4:
                Dscal(MatricesB1T[i]->GetN_Entries(),
                  1.0/press_factor1,
                  MatricesB1T[i]->GetEntries());
                Dscal(MatricesB2T[i]->GetN_Entries(),
                  1.0/press_factor1,
                  MatricesB2T[i]->GetEntries());
		/*    Dscal(MatricesB1[i]->GetN_Entries(),
			  1.0/press_factor1,
			  MatricesB1[i]->GetEntries());
		    Dscal(MatricesB2[i]->GetN_Entries(),
			  1.0/press_factor1,
			  MatricesB2[i]->GetEntries());
		*/
                break;
            }
          }                                       // endfor
        }

        // measure norms of convoluted velocity
        if ((TDatabase::ParamDB->P7==123456789)&&(l==N_SubSteps-1)&&(m%10 == 0))
        {
          // convolute velocity
          ConvoluteVelocityFull(UArray[mg_level-1], uConv);

          fesp[0] = USpaces[mg_level-1];
          fefct[0] = u1ConvArray[mg_level-1];
          fefct[1] = u2ConvArray[mg_level-1];

          aux =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
            TimeNSN_ParamFct2,
            TimeNSN_FEValues2,
            fesp, fefct,
            TimeNSFct2,
            TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
            TimeNSN_Params2, TimeNSBeginParam2);

          // errors
          u1ConvArray[mg_level-1]->GetErrors(ExactU1, 3, TimeNSAllDerivatives,
            2, L2H1Errors,
            NULL, aux, 1, USpaces+mg_level-1, errors);

          u2ConvArray[mg_level-1]->GetErrors(ExactU2, 3, TimeNSAllDerivatives,
            2, L2H1Errors,
            NULL, aux, 1, USpaces+mg_level-1, errors+2);

          OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
          OutPut( "L2(uconv): " << sqrt(errors[0]*errors[0]+errors[2]*errors[2]));
          OutPut( "   H1-semi(uconv): " << sqrt(errors[1]*errors[1]+errors[3]*errors[3])<<endl);

          OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
          OutPut( "kinetic_energy uconv " << (errors[0]*errors[0]+errors[2]*errors[2])/2 );
          OutPut(endl);
        }
      }                          // endfor l (sub steps of fractional step theta)
    }                            // endfor two time discs of adaptive time step control

    if (time_discs==2)
    {
      // compute difference of solutions
      for (i=0;i<N_Unknowns;i++)
        sol[i]-=frac_step_sol[i];

      // compute norms of difference
      fesp[0] = USpaces[mg_level-1];
      fefct[0] = U1Array[mg_level-1];
      fefct[1] = U2Array[mg_level-1];

      aux =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
        TimeNSN_ParamFct2,
        TimeNSN_FEValues2,
        fesp, fefct,
        TimeNSFct2,
        TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
        TimeNSN_Params2, TimeNSBeginParam2);
      // errors
      U1Array[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1Errors,
        NULL, aux, 1, USpaces+mg_level-1, errors);

      U2Array[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1Errors,
        NULL, aux, 1, USpaces+mg_level-1, errors+2);

      PArray[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1Errors,
        NULL, aux, 1, PSpaces+mg_level-1, errors+4);
      // compute L^2 error
      errors[6] = sqrt(errors[0]*errors[0]+errors[2]*errors[2]+errors[4]*errors[4]);

      if (TDatabase::TimeDB->CURRENTTIME< end_time)
        ComputeNewTimeStep(errors[6]);
      // copy solution of fract. step scheme
      memcpy(sol,frac_step_sol,N_Unknowns*SizeOfDouble);
    }                            // adaptive time step control

    // measure errors
    if (TDatabase::ParamDB->MEASURE_ERRORS)
    {
      // errors in the velocity
      fesp[0] = USpaces[mg_level-1];
      fefct[0] = U1Array[mg_level-1];
      fefct[1] = U2Array[mg_level-1];

      aux =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
        TimeNSN_ParamFct2,
        TimeNSN_FEValues2,
        fesp, fefct,
        TimeNSFct2,
        TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
        TimeNSN_Params2, TimeNSBeginParam2);

      // error in first component
      U1Array[mg_level-1]->GetErrors(ExactU1, 3, TimeNSAllDerivatives,
        2, L2H1Errors,
        NULL, aux, 1, USpaces+mg_level-1, errors);

      // error in second component
      U2Array[mg_level-1]->GetErrors(ExactU2, 3, TimeNSAllDerivatives,
        2, L2H1Errors,
        NULL, aux, 1, USpaces+mg_level-1, errors+2);

      // error in velocity in current time step
      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "L2(u): " << sqrt(errors[0]*errors[0]+errors[2]*errors[2]));
      OutPut( "   H1-semi(u):  " << sqrt(errors[1]*errors[1]+errors[3]*errors[3])<<endl);

      // error in L^infty(0,t,L^2)
      if (sqrt(errors[0]*errors[0]+errors[2]*errors[2])
        > l_infty_l_2)
      {
        l_infty_l_2 = sqrt(errors[0]*errors[0]+errors[2]*errors[2]);
        l_infty_l_2_time =  TDatabase::TimeDB->CURRENTTIME;
      }
      OutPut( l_infty_l_2_time <<  " l_infty(L2(u)) " << l_infty_l_2 << endl);

      // error in L^2(0,t,L^2)
      l_2_l_2u += (errors[0]*errors[0] + errors[2]*errors[2]
        +olderror_l_2_l_2u)*
        TDatabase::TimeDB->TIMESTEPLENGTH/2.0;

      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "L2(0,t,L2)(u) : " <<  sqrt(l_2_l_2u) << endl);

      olderror_l_2_l_2u = errors[0]*errors[0] + errors[2]*errors[2];

      //error in L^2(0,t,H^1)
      l_2_h_1u += (errors[1]*errors[1] + errors[3]*errors[3]
        +olderror_l_2_h_1u)*
        TDatabase::TimeDB->TIMESTEPLENGTH/2.0;

      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "L2(0,t,H1-semi)(u) : " << sqrt(l_2_h_1u) << endl);

      olderror_l_2_h_1u = errors[1]*errors[1] + errors[3]*errors[3];

      // error of deformation tensor in L^2(0,t,L^2)
      errors[0] = 0.0;

      UArray[mg_level-1]->GetDeformationTensorErrors
        (ExactU1, ExactU2,
        3, TimeNSAllDerivatives,
        1, DeformationTensorError ,
        NULL, aux, 1, USpaces+mg_level-1, errors);

      l_2_l_2Du += (errors[0]*errors[0] +olderror * olderror)*
        TDatabase::TimeDB->TIMESTEPLENGTH/2.0;

      olderror = errors[0];

      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "L2(Du) : " << errors[0] << " L2(0,t,L2)(Du) " << sqrt(l_2_l_2Du) << endl);

      // kinetic energy

      U1Array[mg_level-1]->GetErrors(ExactNull, 3, TimeNSAllDerivatives,
        2, L2H1Errors,
        NULL, aux, 1, USpaces+mg_level-1, errors);

      U2Array[mg_level-1]->GetErrors(ExactNull, 3, TimeNSAllDerivatives,
        2, L2H1Errors,
        NULL, aux, 1, USpaces+mg_level-1, errors+2);

      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "kinetic energy " << (errors[0]*errors[0]+errors[2]*errors[2])/2 );
      OutPut(endl);

      // error in pressure
      PArray[mg_level-1]->GetErrors(ExactP, 3, TimeNSAllDerivatives,
        2, L2H1Errors,
        NULL, aux, 1, PSpaces+mg_level-1, errors);
      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "L2(p): " << errors[0]);
      OutPut( "   H1-semi(p): " << errors[1]<<endl);
      //OutPut( "L2(p): " << errors[0] << endl);
      //OutPut( "H1(p): " << errors[1] << endl);

      l2p[mg_level-1] = errors[0];
      h1p[mg_level-1] = errors[1];

      if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 4)
      {
        OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
        OutPut( "subgrid dissipation not implemented !!!" << endl);
      }
      else
      {
        // subgrid dissipation
        UArray[mg_level-1]->GetDeformationTensorErrors
          (ExactU1, ExactU2,
          3, TimeNSAllDerivatives,
          1, SubGridDissipation,
          NULL, aux, 1, USpaces+mg_level-1, errors);

        OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
        OutPut( "subgrid dissipation : " << errors[0] << endl);
      }
      UArray[mg_level-1]->GetDeformationTensorErrors
          (ExactNull, ExactNull,
	   3, TimeNSAllDerivatives,
	   2, DivergenceError,
	   NULL, aux, 1, USpaces+mg_level-1, errors);

      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "divergence error (L1/L2) : " << errors[0]*errors[0]  << 
	      " " << errors[1] << endl);
      delete aux;
      #ifdef __MIXINGLAYERSLIP__
      ComputeVorticityDivergence(USpaces[mg_level-1], U1Array[mg_level-1], U2Array[mg_level-1],
        vorticity_space,vorticity,div);
      ComputeVorticiyThickness(Vorticity,errors);
      comp_vort = 1;
      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "vorticity thickness: " << errors[0] << " " <<
        errors[0]/vort_zero << endl);

      // enstrophy
      aux =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
        TimeNSN_ParamFct2,
        TimeNSN_FEValues2,
        fesp, fefct,
        TimeNSFct2,
        TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
        TimeNSN_Params2, TimeNSBeginParam2);

      Vorticity->GetErrors(ExactNull, 3, TimeNSAllDerivatives,
        2, L2H1Errors,
        NULL, aux, 1,  &vorticity_space, errors);

      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "enstrophy " << (errors[0]*errors[0])/2 << endl);
      OutPut(endl);
      delete aux;
      #endif
    }                            // endif MEASURE_ERRORS

    if (TDatabase::ParamDB->CONVOLUTE_SOLUTION)
    {
      // prepare auxiliary problem
      MG->RestrictToAllGrids();
      DiscreteForm = DiscreteFormRHSAuxProblemU;
      fesp[0] = USpaces[mg_level-1];

      fefct[0] = U1Array[mg_level-1];
      fefct[1] = U2Array[mg_level-1];

      ferhs[0] = USpaces[mg_level-1];
      ferhs[1] = USpaces[mg_level-1];

      // assemble rhs
      N_FESpaces = 1;
      N_Rhs = 2;
      N_SquareMatrices = 0;
      N_RectMatrices = 0;

      memset(u_uConv,0,2*N_U*SizeOfDouble);
      RHSs[0] = rhsGL00AuxProblem;
      RHSs[1] = rhsGL00AuxProblem + N_U;

      aux =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
        TimeNSN_ParamFct2,
        TimeNSN_FEValues2,
        fesp, fefct,
        TimeNSFct2,
        TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
        TimeNSN_Params2, TimeNSBeginParam2);

      Assemble2D(N_FESpaces, fesp,
        N_SquareMatrices, SQMATRICES,
        N_RectMatrices, MATRICES,
        N_Rhs, RHSs, ferhs,
        DiscreteForm,
        BoundaryConditionsAuxProblem,
        BoundValuesAuxProblem,
        aux);

      // solve auxiliary problem
      Solver(sqmatrixGL00AuxProblem, RHSs[0], u_uConv,2);

      // error in first component
      u1ConvArray[mg_level-1]->GetErrors(ExactNull, 3, TimeNSAllDerivatives,
        2, L2H1Errors,
        NULL, aux, 1, USpaces+mg_level-1, errors);

      // error in second component
      u2ConvArray[mg_level-1]->GetErrors(ExactNull, 3, TimeNSAllDerivatives,
        2, L2H1Errors,
        NULL, aux, 1, USpaces+mg_level-1, errors+2);

      OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
      OutPut( "conv H1-semi(u): " << sqrt(errors[1]*errors[1]+errors[3]*errors[3]
        ) << endl);
      OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
      OutPut( "conv kinetic energy " << (errors[0]*errors[0]+errors[2]*errors[2]
        )/2<< endl);

      #ifdef __MIXINGLAYERSLIP__
      // compute vorticity
      ComputeVorticityDivergence(velocity_space, u1ConvArray[mg_level-1],
        u2ConvArray[mg_level-1],
        vorticity_space,conv_vort,div);

      ComputeVorticiyThickness(Conv_Vort,errors);
      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "vorticity thickness (uconv): " << errors[0] << " " <<
        errors[0]/vort_zero_conv << endl);

      Conv_Vort->GetErrors(ExactNull, 3, TimeNSAllDerivatives,
        2, L2H1Errors,
        NULL, aux, 1,  &vorticity_space, errors);

      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "enstrophy (uconv) " << (errors[0]*errors[0])/2 << endl);
      for (k=0;k<N_Vort;k++)
        div[k] = conv_vort[k];
      #endif
      delete aux;
    }

    #ifdef __BENCH__

    // compute characteristic values (deltaP, Cd, Cl)
    GetCdCl(U1Array[mg_level-1], U2Array[mg_level-1], PArray[mg_level-1], U1old, U2old, Cd, Cl);

    PArray[mg_level-1]->FindGradient(0.15, 0.2, dP1);
    PArray[mg_level-1]->FindGradient(0.25, 0.2, dP2);

    OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
    OutPut( "C_drag = " << setprecision(16) <<Cd );
    OutPut( " C_lift = " << setprecision(16) << Cl);
    OutPut( " deltaP = " << setprecision(16) << dP1[0] - dP2[0] << endl);
    OutPut( setprecision(7) << endl);
    memcpy(former_sol,sol,2*N_U*SizeOfDouble);
    #endif

    #ifdef __CHANNELSTEP__
    GetReattachmentPoint(U1Array[mg_level-1], reatt_pt);
    OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
    OutPut( "reattachment: " << reatt_pt<< endl);
    #endif

    #ifdef __CHANNELSTEPSLIP__
    GetReattachmentPoint(U1Array[mg_level-1], reatt_pt);
    OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
    OutPut( "reattachment: " << reatt_pt<< endl);
    #endif

    if ((TDatabase::ParamDB->WRITE_GRAPE)||(TDatabase::ParamDB->WRITE_GNU)
	||(TDatabase::ParamDB->WRITE_GMV)||(TDatabase::ParamDB->WRITE_VTK)
	|| (TDatabase::ParamDB->SAVE_DATA))
    {
      if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      {
        StreamFunction(USpaces[mg_level-1], sol, sol+N_Uarray[mg_level-1],
          PsiSpaces[LEVELS-1], psi);

        Prolongate(USpaces[mg_level-1], PsiSpaces[LEVELS-1], sol,
          uconf->GetValues(), app);

        Prolongate(USpaces[mg_level-1], PsiSpaces[LEVELS-1],
          sol+N_Uarray[mg_level-1], uconf->GetValues()+N_V, app);
        if (!comp_vort)
          ComputeVorticityDivergence(USpaces[mg_level-1], U1Array[mg_level-1], U2Array[mg_level-1],
            vorticity_space,vorticity,div);
      }
    }

    if(TDatabase::ParamDB->WRITE_GRAPE)
    {
      if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      {
        os.seekp(std::ios::beg);
        os << GrapeBaseName << m << ".dat" << ends;
        Output->WriteGrape(os.str().c_str());
      }
    }

    if(TDatabase::ParamDB->WRITE_GMV)
    {
	if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
	{
	    os.seekp(std::ios::beg);
	    os << GmvBaseName << m << ".gmv" << ends;
	    Output->WriteGMV(os.str().c_str());
	}
    }
    
    if(TDatabase::ParamDB->WRITE_GNU)
    {
      if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      {
        os.seekp(std::ios::beg);
        os << GnuBaseName << m << ".gnu" << ends;
        Output->WriteGnuplot(os.str().c_str());

        os.seekp(std::ios::beg);
        os << GnuBaseName << m << ".psi" << ends;
        Output->WriteGNU_iso(os.str().c_str(),1);
        N_GNU_images++;
      }
    }
 
   if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
    {
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
    }
 

    if (TDatabase::ParamDB->SAVE_DATA)
    {
	save_sol[0] = sol;
	save_N_Unknowns[0] = N_Unknowns;
	SaveData(SaveDataFileName,1,save_sol,save_N_Unknowns);
    }	    

    if ((TDatabase::ParamDB->WRITE_GRAPE)||
	(TDatabase::ParamDB->WRITE_GMV) ||
	(TDatabase::ParamDB->WRITE_VTK))
	N_GRAPE_images++;
    /* if(j == 0)
    {
      // stationary limit is reached
      OutPut("stationary limit reached." << endl);
      break;
      }*/
    comp_vort =0;
  }                              // while

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
      os << GmvBaseName << "end." << m << ".gmv" << ends;
      Output->WriteGMV(os.str().c_str());
  }

  if(TDatabase::ParamDB->WRITE_VTK)
  {
      os.seekp(std::ios::beg);
      os << VtkBaseName << "end." << m << ".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
  }

  if(TDatabase::ParamDB->WRITE_GNU)
  {
    os.seekp(std::ios::beg);
    os << GnuBaseName <<  m << ".gnu" << ends;
    Output->WriteGnuplot(os.str().c_str());
    os.seekp(std::ios::beg);
    os << GnuBaseName <<  m << ".psi" << ends;
    Output->WriteGNU_iso(os.str().c_str(),1);
  }
  t4 =  GetTime();
  total_time += t4 - t3;
  OutPut("total running time: " << total_time << endl);
  CloseFiles();
  return 0;
}
