// =======================================================================
//
// Purpose:     main program
//
// Author:      Volker John 2000/08/25
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
#include <TNSE3D_ParamRout.h>
#include <BoundFace.h>
#include <Collection.h>
#include <RefDesc.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include <Upwind3D.h>
#include <NSE_MultiGrid.h>
#include <NSE_MGLevel1.h>
#include <NSE_MGLevel2.h>
#include <NSE_MGLevel3.h>
#include <NSE_MGLevel4.h>
#include <NSE_MGLevel14.h>
#include <Convolution.h>
#include <ItMethod.h>
#include <FixedPointIte.h>
#include <MultiGridIte.h>
#include <FgmresIte.h>

#include <RationalLES.h>
#include <VMS.h>
#include <RFB.h>

#include <MultiGrid3D.h>
#include <MGLevel3D.h>

#include <MainUtilities.h>
#include <TimeUtilities.h>

#define AMG 0
#define GMG 1

// =======================================================================
// include current example
// =======================================================================
//#include "../Examples/NSE_3D/AnsatzQuadLin.h"
//#include "../Examples/NSE_3D/BSExample.h"
//#include "../Examples/TNSE_3D/AnsatzLinConst.h"
//#include "../Examples/TNSE_3D/DrivenCavity3D.h"
//#include "../Examples/TNSE_3D/TimeBSExample.h"
//#include "../Examples/TNSE_3D/Channel3D.h"
//#include "../Examples/TNSE_3D/AnsatzQuadLin_Taylor.h"
//#include "../Examples/TNSE_3D/Polynom_Taylor.h"
//#include "../Examples/TNSE_3D/Polynom_Taylor.01.h"
//#include "../Examples/TNSE_3D/Bench3DQuaderInst1.h"
//#include "../Examples/TNSE_3D/Bench3DQuaderLongInst.h"
//#include "../Examples/TNSE_3D/Bench3DQuaderInstSin.h"
//#include "../Examples/TNSE_3D/Bench3DQuaderInstSinVariableInflow.h"
//#include "../Examples/TNSE_3D/Bench3DQuaderInstSin_NoDrag.h"
//#include "../Examples/TNSE_3D/ChannelQuadLong.h"
//#include "../Examples/TNSE_3D/ChannelQuadLongStart.h"
//#include "../Examples/TNSE_3D/ChannelQuadLongDiri.h"
//#include "../Examples/TNSE_3D/ChannelSlip3D.h"
//#include "../Examples/TNSE_3D/UnitCubeSlip3D.h"
//#include "../Examples/TNSE_3D/LES_Bench3DQuaderNeum.h"
//#include "../Examples/TNSE_3D/LES_Bench3DQuaderNeumStart.h"
//#include "../Examples/TNSE_3D/ChannelStepSlip3D.h"
//#include "../Examples/TNSE_3D/ChannelStepSlipFreeSlip.h"
//#include "../Examples/TNSE_3D/ChannelStepSlipTopFreeSlipLatDiri.h"
//#include "../Examples/TNSE_3D/ChannelSlipDiri3D.h"
//#include "../Examples/TNSE_3D/ChannelSlipDiriNoise3D.h"
//#include "Examples/NSE_3D/Bench3DQuaderStatNeum.h"
//#include "../Examples/TNSE_3D/MixingLayer3D.h"
// #include "../Examples/TNSE_3D/MixingLayerSlip3D.h"
//#include "../Examples/TNSE_3D/MixingLayerDiri3D.h"
// #include "../Examples/TNSE_3D/Channel_Tobias.h"
//#include "../Examples/TNSE_3D/Channel_Tobias.nonzero_rhs.h"
//#include "../Examples/TNSE_3D/Beltrami.h"
//#include "../Examples/TNSE_3D/Channel_Carolina.h"
#include "../Examples/TNSE_3D/ChannelTau180.h"
//#include "../Examples/TNSE_3D/WallMountedCube.h"
//#include "../Examples/TNSE_3D/CylinderSquare22000.h"
//#include "../Examples/TNSE_3D/Calotte.h"
//#include "../Examples/TNSE_3D/windtunnel_00.h"
//#include "../Examples/TNSE_3D/DrivenCavity3D_Bulk.h"
//#include "../Examples/TNSE_3D/windtunnel_fine.h"

// ======================================================================
// utilities for main program
// ======================================================================

int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase3D *FEDatabase = new TFEDatabase3D();
  TCollection *coll, *mortarcoll = NULL, **CollArray, *coll0;
  TBaseCell *cell, *neigh;
  TFESpace3D *velocity_space, *pressure_space, *convolution_space;
  TFESpace3D **USpaces, **PSpaces, **duConvSpaces, **uConvSpaces;
  TFESpace3D **pConvSpaces, **ProjectionSpaces;
  TFESpace3D *vorticity_space, *projection_space, *size_small_scales_fesp, *label_space_fesp;
  TOutput3D *Output;

  double *B, *rhs, *sol, *oldsol, tol, tolerance, *defect, *startsol, *frac_step_sol;
  double *oldrhs, *itmethod_sol, *itmethod_rhs, *sol_timestep_m1;
  double *vorticity, *div, *newton, *old_small_scales;
  double *solGL00AuxProblem, *rhsGL00AuxProblem, *LESModelRhs;
  int i,j,k,l,m,n, N_, Len, low, N_cells;
  int N_Rows, N_Columns, N_U, N_P, N_Unknowns, N_V, N_sub;
  double *l2u1, *l2u2, *l2u3, *h1u1, *h1u2, *h1u3;
  double *l2p, *h1p, *sd, *l_inf, vect[3], exactvect[3];
  int which, *permutation, N_GRAPE_images=0,  N_GNU_images=0;
  double DiffL2, DiffH1, t;
  char *PRM, *GEO;
  int LEVELS, ORDER, order, cell_measure;
  int ret, pde;
  double negPower;
  double x,y,max,min,sum;
  double tau1, tau2;
  double errors[9], p1, p2, *save_sol[1];
  double t1, t2, res, res2, oldres, solver_time, residual, oldresidual;
  double impuls_residual,limit,linredfac, solver_time_curr, t11, t22, t111, t222;
  int N_LinIter, N_LinIterCurr, N_LinIterCurrIte, N_SubSteps, N_Active, n_aux;
  double gamma, tau, oldtau, total_time, t3, t4, t5, ass_time;
  double umfpack_flag=-2;
  int *RowPtr,nonlinite, save_N_Unknowns[1];

  std::ostringstream os;
  char *PsBaseName, *GrapeBaseName, *ReadGrapeBaseName, *GMVBaseName, *VTKBaseName;
  char *SaveDataFileName,  *ReadDataFileName;

  double *val, cd, cl;
  TFEFunction3D *u1, *u2, *u3, *p, *fefct[12], *Approx;
  TFEFunction3D *du11Conv, *du12Conv, *du13Conv, *du22Conv, *du23Conv, *du33Conv;
  TFEFunction3D *u1Conv, *u2Conv, *u3Conv, *u4Conv, *u5Conv, *u6Conv;
  TFEFunction3D **U1Array, **U2Array, **U3Array, **AuxFEFunctArray;
  TFEFunction3D **PArray;
  TFEFunction3D *vort1, *vort2, *vort3, *Divergence, *Vort_z;
  TFEFunction3D *Vort_x, *Vort_y;
  TFEVectFunct3D *u, **UArray, *uconf, *duConv, **duConvArray;
  TFEVectFunct3D *uConv, **uConvArray, **AuxFEVectFunctArray;
  TFEVectFunct3D *Vorticity, *vms_projection_fe;
  TFEFunction3D **du11ConvArray, **du12ConvArray, **du13ConvArray;
  TFEFunction3D **du22ConvArray, **du23ConvArray, **du33ConvArray;
  TFEFunction3D **u1ConvArray, **u2ConvArray, **u3ConvArray;
  TFEFunction3D **u4ConvArray, **u5ConvArray, **u6ConvArray;
  TFEVectFunct3D *GL00AuxProblemSol, **GL00AuxProblemSolArray;
  TFEFunction3D *GL00AuxProblemSol11, *GL00AuxProblemSol12;
  TFEFunction3D *GL00AuxProblemSol13, *GL00AuxProblemSol22;
  TFEFunction3D *GL00AuxProblemSol23, *GL00AuxProblemSol33;
  TFEFunction3D **GL00AuxProblemSol11Array, **GL00AuxProblemSol12Array;
  TFEFunction3D **GL00AuxProblemSol13Array, **GL00AuxProblemSol22Array;
  TFEFunction3D **GL00AuxProblemSol23Array, **GL00AuxProblemSol33Array;
  TFEFunction3D *pConv, **pConvArray;
  TFEFunction3D *vms_proj_11, *vms_proj_12, *vms_proj_13, *vms_proj_22;
  TFEFunction3D *vms_proj_23, *vms_proj_33, *size_small_scales_fefct, *label_space_fefct;
  TFESpace3D *fesp[5], *ferhs[6];
  FE3D *fes, *fes1;
  double delta, end_time, l_infty_l_2 = 0, l_infty_l_2_time=-4711.0;
  double olderror = 0, l_2_l_2Du=0, l_2_l_2u=0 , olderror_l_2_l_2u=0;
  double l_2_h_1u=0, olderror_l_2_h_1u=0;

  TAuxParam3D *aux, *auxn;

  TSquareStructure3D *sqstructureA, *sqstructureC;
  TStructure3D *structureB, *structureBT;
  TSquareMatrix3D *sqmatrixA, *SQMATRICES[18];
  TSquareMatrix3D *sqmatrixA11, *sqmatrixA12, *sqmatrixA13;
  TSquareMatrix3D *sqmatrixA21, *sqmatrixA22, *sqmatrixA23;
  TSquareMatrix3D *sqmatrixA31, *sqmatrixA32, *sqmatrixA33;
  TSquareMatrix3D *sqmatrixM;
  TSquareMatrix3D *sqmatrixM11, *sqmatrixM12, *sqmatrixM13;
  TSquareMatrix3D *sqmatrixM21, *sqmatrixM22, *sqmatrixM23;
  TSquareMatrix3D *sqmatrixM31, *sqmatrixM32, *sqmatrixM33;
  TSquareMatrix3D **MatricesM11, **MatricesM12, **MatricesM13;
  TSquareMatrix3D **MatricesM21, **MatricesM22, **MatricesM23;
  TSquareMatrix3D **MatricesM31, **MatricesM32, **MatricesM33;
  TSquareMatrix3D *sqmatrixC, *sqmatrixK;
  TSquareMatrix3D **MatricesA, **MatricesM;
  TSquareMatrix3D **MatricesK, **MatricesC;
  TSquareMatrix3D *sqmatrixS11, *sqmatrixS12, *sqmatrixS13;
  TSquareMatrix3D *sqmatrixS21, *sqmatrixS22, *sqmatrixS23;
  TSquareMatrix3D *sqmatrixS31, *sqmatrixS32, *sqmatrixS33;
  TSquareMatrix3D **MatricesA11, **MatricesA12, **MatricesA13;
  TSquareMatrix3D **MatricesA21, **MatricesA22, **MatricesA23;
  TSquareMatrix3D **MatricesA31, **MatricesA32, **MatricesA33;
  TSquareMatrix3D **MatricesS11, **MatricesS12, **MatricesS13;
  TSquareMatrix3D **MatricesS21, **MatricesS22, **MatricesS23;
  TSquareMatrix3D **MatricesS31, **MatricesS32, **MatricesS33;
  TSquareMatrix3D *sqmatrixGL00AuxProblem,**MatricesGL00AuxProblem;
  TMatrix3D *matrixB1, *matrixB2, *matrixB3, *MATRICES[12];
  TMatrix3D *matrixB1T, *matrixB2T, *matrixB3T;
  TMatrix3D **MatricesB1, **MatricesB2, **MatricesB1T, **MatricesB2T;
  TMatrix3D **MatricesB3, **MatricesB3T;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TMatrix **matrices = (TMatrix **)MATRICES;

  TSquareStructure3D *sqstructureL;
  TStructure3D *structure_tilde_G, *structure_G;
  TSquareMatrix3D *sqmatrixL, **MatricesL;
  TMatrix3D *matrix_tilde_G11, *matrix_tilde_G22, *matrix_tilde_G33;
  TMatrix3D *matrix_G11, *matrix_G22, *matrix_G33;
  TMatrix3D **Matrices_tilde_G11, **Matrices_tilde_G22, **Matrices_tilde_G33;
  TMatrix3D **Matrices_G11, **Matrices_G22, **Matrices_G33;
  int N_L;
  double *rhs_vms_expl, *vms_projection;

  MatVecProc *MatVect;
  DefectProc *Defect;

  double **RhsArray;

  TNSE_MGLevel *MGLevel;
  TNSE_MultiGrid *MG;
  TMGLevel3D *MGLevelGL00AuxProblem;
  TMultiGrid3D *MGGL00AuxProblem;

  double *RHSs[6];
  int *N_Uarray, *N_Parray;

  TDiscreteForm3D *DiscreteFormGalerkin;
  TDiscreteForm3D *DiscreteFormClassicalLES;
  TDiscreteForm3D *DiscreteFormGL00Convolution;
  TDiscreteForm3D *DiscreteFormGL00AuxProblem;
  TDiscreteForm3D *DiscreteFormUpwind;
  TDiscreteForm3D *DiscreteFormUpwindNC;
  TDiscreteForm3D *DiscreteFormSmagorinsky;
  TDiscreteForm3D *DiscreteFormVMS_Projection;
  TDiscreteForm3D *DiscreteFormVMS_SUPG;

  TDiscreteForm3D *DiscreteFormNLGalerkin;
  TDiscreteForm3D *DiscreteFormNLUpwind;
  TDiscreteForm3D *DiscreteFormNLUpwindNC;
  TDiscreteForm3D *DiscreteFormNLClassicalLES;
  TDiscreteForm3D *DiscreteFormNLGL00Convolution;
  TDiscreteForm3D *DiscreteFormNLGL00AuxProblem;
  TDiscreteForm3D *DiscreteFormNLSmagorinsky;
  TDiscreteForm3D *DiscreteFormNLVMS_Projection;
  TDiscreteForm3D *DiscreteFormNLVMS_ProjectionExpl;
  TDiscreteForm3D *DiscreteFormNLVMS_RFBExplRhs;
  TDiscreteForm3D *DiscreteFormNLVMS_SUPG;

  TDiscreteForm3D *DiscreteFormRHS;
  TDiscreteForm3D *DiscreteFormRHSClassicalLES;
  TDiscreteForm3D *DiscreteFormRHSLESModel;
  TDiscreteForm3D *DiscreteFormRHSSUPG;
  TDiscreteForm3D *DiscreteFormMatrixGL00AuxProblem;
  TDiscreteForm3D *DiscreteFormGL00AuxProblemRHS;
  TDiscreteForm3D *DiscreteFormRHSAuxProblemU;
  TDiscreteForm3D *DiscreteFormMatrixAuxProblemU;

  TDiscreteForm3D *DiscreteFormRHSNewton;
  TDiscreteForm3D *DiscreteFormRHSNewtonNL;

  TDiscreteForm3D *DiscreteFormC;
  TDiscreteForm3D *DiscreteFormJ;

  TDiscreteForm3D *DiscreteForm;

  int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;

  BoundCondFunct3D *BoundaryConditions[6], *BoundaryConditionsAuxProblem[6];
  BoundValueFunct3D *BoundValues[6], *BoundValuesAuxProblem[6];
  CoeffFct3D *Coefficients[1];
  double average, hmin, hmax;

  TItMethod *itmethod, *prec;
  int Max_It, FirstSolve;
  double omega, alpha[2], divergence;
  int N_Paramters=10, methods, time_discs, **downwind;
  int *neum_to_diri, N_neum_to_diri = 0, *neum_bdry, N_neum_bdry = 0;
  double Parameters[10];
  double theta1, theta2, theta3, theta4;

  int N_UConv, level_down, ii, fully_implicit = 0, N_UuConv;
  int N_PpConv, N_Vort;
  double *u_conv, *auxConv, *du_tensor, *u_uConv, *p_pConv, reatt_pt;
  int mg_level,mg_type,CurrentDiscType, last_sq, mid_sq, step_length;
  int velocity_space_code, pressure_space_code;
  int very_first_time=0, very_first_supg = 1, zerostart;
  int memory[4], *left_dof, *right_dof, dof_length, dof_length_new;
  int  *left_dof1, *right_dof1, dof_length_new1,  comp_vort, *direction_dof;
  int mixing_layer_galerkin = 0, N_vort, turb_visc_expl_change = 0;
  double vort_zero, vort_zero_conv, moment_zero, moment_zero_conv;
  double *sol_vort_tmp;
  double *x_dof, *y_dof, *z_dof, *coord_z_layers, *mean_velocity, *rms_velocity;
  double *mean_velocity_u2, *mean_velocity_u3, *dmean_velocity;
  double *rms_velocity2, *rms_velocity3;
  double *rms_velocity1_type1, *rms_velocity2_type1, *rms_velocity3_type1;
  double *R_xx, *R_yy, *R_zz, *R_xy, *R_xz, *R_yz, *u1x, *projection_u1x;
  double *A_M_xx, *A_M_yy, *A_M_zz, *A_M_xy, *A_M_xz, *A_M_yz;
  double *size_small_scales,*label_space;
  double mean, largest_size, mean_time_average, max_time_average;
  int N_z_layers;
  int smoothing_depth;
  struct mallinfo MALLINFO;

#ifdef __BENCH__
  double Cd, Cl, dP1[4], dP2[4], *former_sol, velo_friction[4], *press_cyl, *center_velo;
  double *cyl_velo;
  int fric_count = 0, N_press_cyl, press_count = 0, N_center_velo, center_velo_count = 0;
  int  N_cyl_velo, cyl_velo_count = 0;
  TFEFunction3D *U1old, *U2old;
#endif
#ifdef  __CHANNEL_OBSTACLE__
  double Cd, Cl, dP1[4], dP2[4] *former_sol;
  TFEFunction3D *U1old, *U2old;
#endif

  // Strings
  char ReadinDat[] = "readin.dat";
  char NameString[] = "name";
  char UString[] = "u";
  char PString[] = "p";
  char VortString[] = "vort";
  char UConvString[] = "u_conv";
  char AuxProblemString[] = "AuxProblem";
  char PConvString[] = "p_conv";
  char VorticityString[] = "vorticity";
  char DivergenceString[] = "divergence";
  char SmallScaleString[] = "smallscales";
  char LabelSpaceString[] = "space";

  os << " ";

  //======================================================================
  // read parameter file
  //======================================================================
  total_time = GetTime();
  ass_time = 0;
  if(argc>=2)
    ret=Domain->ReadParam(argv[1]);
  else
    ret=Domain->ReadParam(ReadinDat);

  if(ret==-1)
  {
    OutPut("No data readin file found !!!"<<endl);
    exit(-1);
  }

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  // just to measure the error of the convolved velocity
  if (TDatabase::ParamDB->P7==123456789)
    TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE=4;

  Database->CheckParameterConsistencyNSE();
  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();

  //======================================================================
  // copy read parameters into local variables
  //======================================================================

  if( (TDatabase::ParamDB->DISCTYPE==5) )
  {
    OutPut("DISCTYPE 5 NOT IMPLEMENTED!" << endl);
    exit(4711);
  }

  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;

  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  GrapeBaseName = TDatabase::ParamDB->GRAPEBASENAME;
  ReadGrapeBaseName = TDatabase::ParamDB->READGRAPEBASENAME;
  GMVBaseName = TDatabase::ParamDB->GMVBASENAME;
  VTKBaseName = TDatabase::ParamDB->VTKBASENAME;
  SaveDataFileName = TDatabase::ParamDB->SAVE_DATA_FILENAME;
  ReadDataFileName = TDatabase::ParamDB->READ_DATA_FILENAME;

  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SADDLE;
  if (mg_type)
  {
    if (mg_type==2)
      mg_level = 2;
    else
      mg_level = 1;
  }
  else
    mg_level = 0;

  LEVELS = TDatabase::ParamDB->LEVELS;

  ///////////////////////////////////////////////////////////////////////////////
  //
  // allocation of arrays
  //
  ///////////////////////////////////////////////////////////////////////////////

  l2u1 = new double[LEVELS+1];
  l2u2 = new double[LEVELS+1];
  l2u3 = new double[LEVELS+1];
  l2p = new double[LEVELS+1];
  h1u1 = new double[LEVELS+1];
  h1u2 = new double[LEVELS+1];
  h1u3 = new double[LEVELS+1];
  h1p = new double[LEVELS+1];
  sd = new double[LEVELS+1];
  l_inf = new double[LEVELS+1];

  U1Array = new TFEFunction3D*[LEVELS+1];
  U2Array = new TFEFunction3D*[LEVELS+1];
  U3Array = new TFEFunction3D*[LEVELS+1];
  PArray = new TFEFunction3D*[LEVELS+1];
  UArray = new TFEVectFunct3D*[LEVELS+1];
  duConvArray = new TFEVectFunct3D*[LEVELS+1];
  du11ConvArray = new TFEFunction3D*[LEVELS+1];
  du12ConvArray = new TFEFunction3D*[LEVELS+1];
  du13ConvArray = new TFEFunction3D*[LEVELS+1];
  du22ConvArray = new TFEFunction3D*[LEVELS+1];
  du23ConvArray = new TFEFunction3D*[LEVELS+1];
  du33ConvArray = new TFEFunction3D*[LEVELS+1];
  uConvArray = new TFEVectFunct3D*[LEVELS+1];
  u1ConvArray = new TFEFunction3D*[LEVELS+1];
  u2ConvArray = new TFEFunction3D*[LEVELS+1];
  u3ConvArray = new TFEFunction3D*[LEVELS+1];
  u4ConvArray = new TFEFunction3D*[LEVELS+1];
  u5ConvArray = new TFEFunction3D*[LEVELS+1];
  u6ConvArray = new TFEFunction3D*[LEVELS+1];
  GL00AuxProblemSolArray = new TFEVectFunct3D*[LEVELS+1];
  GL00AuxProblemSol11Array = new TFEFunction3D*[LEVELS+1];
  GL00AuxProblemSol12Array = new TFEFunction3D*[LEVELS+1];
  GL00AuxProblemSol13Array = new TFEFunction3D*[LEVELS+1];
  GL00AuxProblemSol22Array = new TFEFunction3D*[LEVELS+1];
  GL00AuxProblemSol23Array = new TFEFunction3D*[LEVELS+1];
  GL00AuxProblemSol33Array = new TFEFunction3D*[LEVELS+1];
  pConvArray = new TFEFunction3D*[LEVELS+1];

  RhsArray = new double*[LEVELS+1];
  N_Uarray = new int[LEVELS+1];
  N_Parray = new int[LEVELS+1];

  USpaces = new TFESpace3D*[LEVELS+1];
  PSpaces = new TFESpace3D*[LEVELS+1];
  duConvSpaces = new TFESpace3D*[LEVELS+1];
  uConvSpaces = new TFESpace3D*[LEVELS+1];
  pConvSpaces = new TFESpace3D*[LEVELS+1];
  ProjectionSpaces = new TFESpace3D*[LEVELS+1];

  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      MatricesA = new TSquareMatrix3D*[LEVELS+1];
      MatricesM = new TSquareMatrix3D*[LEVELS+1];

      MatricesB1 = new TMatrix3D*[LEVELS+1];
      MatricesB2 = new TMatrix3D*[LEVELS+1];
      MatricesB3 = new TMatrix3D*[LEVELS+1];
      MatVect = MatVect_NSE1;
      Defect = Defect_NSE1;
      break;

    case 2:
      MatricesA = new TSquareMatrix3D*[LEVELS+1];
      MatricesM = new TSquareMatrix3D*[LEVELS+1];

      MatricesB1 = new TMatrix3D*[LEVELS+1];
      MatricesB2 = new TMatrix3D*[LEVELS+1];
      MatricesB3 = new TMatrix3D*[LEVELS+1];
      MatricesB1T = new TMatrix3D*[LEVELS+1];
      MatricesB2T = new TMatrix3D*[LEVELS+1];
      MatricesB3T = new TMatrix3D*[LEVELS+1];
      MatVect = MatVect_NSE2;
      Defect = Defect_NSE2;
      break;

    case 3:
      MatricesA11 = new TSquareMatrix3D*[LEVELS+1];
      MatricesA12 = new TSquareMatrix3D*[LEVELS+1];
      MatricesA13 = new TSquareMatrix3D*[LEVELS+1];
      MatricesA21 = new TSquareMatrix3D*[LEVELS+1];
      MatricesA22 = new TSquareMatrix3D*[LEVELS+1];
      MatricesA23 = new TSquareMatrix3D*[LEVELS+1];
      MatricesA31 = new TSquareMatrix3D*[LEVELS+1];
      MatricesA32 = new TSquareMatrix3D*[LEVELS+1];
      MatricesA33 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM11 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM12 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM13 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM21 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM22 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM23 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM31 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM32 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM33 = new TSquareMatrix3D*[LEVELS+1];

      MatricesB1 = new TMatrix3D*[LEVELS+1];
      MatricesB2 = new TMatrix3D*[LEVELS+1];
      MatricesB3 = new TMatrix3D*[LEVELS+1];
      MatVect = MatVect_NSE3;
      Defect = Defect_NSE3;
      break;

    case 4:
    case 14:
      MatricesA11 = new TSquareMatrix3D*[LEVELS+1];
      MatricesA12 = new TSquareMatrix3D*[LEVELS+1];
      MatricesA13 = new TSquareMatrix3D*[LEVELS+1];
      MatricesA21 = new TSquareMatrix3D*[LEVELS+1];
      MatricesA22 = new TSquareMatrix3D*[LEVELS+1];
      MatricesA23 = new TSquareMatrix3D*[LEVELS+1];
      MatricesA31 = new TSquareMatrix3D*[LEVELS+1];
      MatricesA32 = new TSquareMatrix3D*[LEVELS+1];
      MatricesA33 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM11 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM12 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM13 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM21 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM22 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM23 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM31 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM32 = new TSquareMatrix3D*[LEVELS+1];
      MatricesM33 = new TSquareMatrix3D*[LEVELS+1];

      MatricesB1 = new TMatrix3D*[LEVELS+1];
      MatricesB2 = new TMatrix3D*[LEVELS+1];
      MatricesB3 = new TMatrix3D*[LEVELS+1];
      MatricesB1T = new TMatrix3D*[LEVELS+1];
      MatricesB2T = new TMatrix3D*[LEVELS+1];
      MatricesB3T = new TMatrix3D*[LEVELS+1];
      MatVect = MatVect_NSE4;
      Defect = Defect_NSE4;

      if (TDatabase::ParamDB->NSTYPE == 14)
      {
        MatricesC = new TSquareMatrix3D*[LEVELS+1];
        MatVect = MatVect_NSE14;
        Defect = Defect_NSE14;
      }
      if (TDatabase::ParamDB->DISCTYPE == SDFEM)
      {
        MatricesK = new TSquareMatrix3D*[LEVELS+1];
        MatricesS11 = new TSquareMatrix3D*[LEVELS+1];
        MatricesS12 = new TSquareMatrix3D*[LEVELS+1];
        MatricesS13 = new TSquareMatrix3D*[LEVELS+1];
        MatricesS21 = new TSquareMatrix3D*[LEVELS+1];
        MatricesS22 = new TSquareMatrix3D*[LEVELS+1];
        MatricesS23 = new TSquareMatrix3D*[LEVELS+1];
        MatricesS31 = new TSquareMatrix3D*[LEVELS+1];
        MatricesS32 = new TSquareMatrix3D*[LEVELS+1];
        MatricesS33 = new TSquareMatrix3D*[LEVELS+1];
      }

      break;
  }                                               // endswitch

  // matrices for VMS_PROJECTION
  if ((TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION) ||
    (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION_EXPL) ||
    (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION_SD))
  {
    MatricesL = new TSquareMatrix3D*[LEVELS+1];
    Matrices_tilde_G11 = new TMatrix3D*[LEVELS+1];
    Matrices_tilde_G22 = new TMatrix3D*[LEVELS+1];
    Matrices_tilde_G33 = new TMatrix3D*[LEVELS+1];
    Matrices_G11 = new TMatrix3D*[LEVELS+1];
    Matrices_G22 = new TMatrix3D*[LEVELS+1];
    Matrices_G33 = new TMatrix3D*[LEVELS+1];
  }

  // array of matrices for auxiliary problem in Galdi/Layton model
  if ((TDatabase::ParamDB->DISCTYPE==GL00_AUX_PROBLEM)
    ||(TDatabase::ParamDB->CONVOLUTE_SOLUTION)
    ||(TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE==4))
  {
    MatricesGL00AuxProblem = new TSquareMatrix3D*[LEVELS+1];
  }

  downwind = new int*[LEVELS+1];

  //======================================================================
  // creating discrete forms
  //======================================================================

  InitializeDiscreteForms(DiscreteFormGalerkin,DiscreteFormUpwind,
    DiscreteFormUpwindNC,
    DiscreteFormSmagorinsky,DiscreteFormClassicalLES,
    DiscreteFormGL00Convolution,DiscreteFormGL00AuxProblem,
    DiscreteFormVMS_Projection,
    DiscreteFormVMS_SUPG,
    DiscreteFormNLGalerkin,
    DiscreteFormNLUpwind, DiscreteFormNLUpwindNC, DiscreteFormNLSmagorinsky,
    DiscreteFormNLClassicalLES,
    DiscreteFormNLGL00Convolution,
    DiscreteFormNLGL00AuxProblem,
    DiscreteFormNLVMS_Projection,
    DiscreteFormNLVMS_ProjectionExpl,
    DiscreteFormNLVMS_RFBExplRhs,
    DiscreteFormNLVMS_SUPG,
    DiscreteFormRHS,
    DiscreteFormRHSClassicalLES,
    DiscreteFormRHSLESModel,
    DiscreteFormRHSSUPG,
    DiscreteFormMatrixGL00AuxProblem,
    DiscreteFormGL00AuxProblemRHS,
    DiscreteFormMatrixAuxProblemU,
    DiscreteFormRHSAuxProblemU,
    DiscreteFormRHSNewton,
    DiscreteFormRHSNewtonNL,
    DiscreteFormC, DiscreteFormJ,
    LinCoeffs, TDatabase::ParamDB->NSTYPE);

  // correct assembling routines have been assigned,
  // all other stuff is the same
  if(TDatabase::ParamDB->DISCTYPE==VMS_PROJECTION_SD)
  {
    TDatabase::ParamDB->DISCTYPE=VMS_PROJECTION;
    OutPut("Change internally DISCTYPE to " << TDatabase::ParamDB->DISCTYPE
      << " since correct DiscreteForms have been assigned"<< endl);
  }

  //======================================================================
  // read boundary parameterization and initialize coarse grid
  //======================================================================
  Domain->Init(PRM, GEO);

  //  Domain->ConvertQuadToTri(TDatabase::ParamDB->CONVERT_QUAD_TO_TRI);

  BoundaryConditions[0] = BoundCondition;
  BoundaryConditions[1] = BoundCondition;
  BoundaryConditions[2] = BoundCondition;
  BoundaryConditions[3] = BoundaryConditionNewton;
  BoundaryConditions[4] = BoundaryConditionNewton;
  BoundaryConditions[5] = BoundaryConditionNewton;

  BoundValues[0] = U1BoundValue;
  BoundValues[1] = U2BoundValue;
  BoundValues[2] = U3BoundValue;
  BoundValues[3] = BoundaryValueNewton;
  BoundValues[4] = BoundaryValueNewton;
  BoundValues[5] = BoundaryValueNewton;

  BoundaryConditionsAuxProblem[0] = BoundConditionAuxProblem;
  BoundaryConditionsAuxProblem[1] = BoundConditionAuxProblem;
  BoundaryConditionsAuxProblem[2] = BoundConditionAuxProblem;
  BoundaryConditionsAuxProblem[3] = BoundConditionAuxProblem;
  BoundaryConditionsAuxProblem[4] = BoundConditionAuxProblem;
  BoundaryConditionsAuxProblem[5] = BoundConditionAuxProblem;

  BoundValuesAuxProblem[0] = BoundValueAuxProblem;
  BoundValuesAuxProblem[1] = BoundValueAuxProblem;
  BoundValuesAuxProblem[2] = BoundValueAuxProblem;
  BoundValuesAuxProblem[3] = BoundValueAuxProblem;
  BoundValuesAuxProblem[4] = BoundValueAuxProblem;
  BoundValuesAuxProblem[5] = BoundValueAuxProblem;

  Coefficients[0] = LinCoeffs;
  // create coarsest grid for multigrid method
  for(i=0;i<TDatabase::ParamDB->SC_COARSEST_LEVEL_SADDLE;i++)
  {
    Domain->RegRefineAll();
#ifdef __CALOTTE__
    MakeInAndOutflowCircular(Domain->GetCollection(It_Finest, 0));
#endif
  }
  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
  omega = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SADDLE;
  divergence = TDatabase::ParamDB->SC_DIV_FACTOR;
  nonlinite = TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE;

  Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
  Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
  Parameters[2] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;

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
    MGGL00AuxProblem = new TMultiGrid3D(i, N_Paramters, Parameters);
  }

  mg_level = LEVELS+mg_level;
  /*
    #ifdef __CHANNEL_TAU180__
    // move vertices of the grid
    coll = Domain->GetCollection(It_Finest, 0);
    SetZCoordinates(coll, 0);
    Domain->MakeBdParamsConsistent(coll);
    CheckZCoordinates(coll, 0);
    delete coll;
    #endif

    #ifdef __CHANNEL_CAROLINA__
  // move vertices of the grid
  coll = Domain->GetCollection(It_Finest, 0);
  SetZCoordinates(coll, 0);
  Domain->MakeBdParamsConsistent(coll);
  CheckZCoordinates(coll, 0);
  delete coll;
  #endif
  */
#ifdef __WALL_MOUNTED_CUBE__
  TBaseCell **Cells;
  int N_CellsOld, N_CellsNew, N_Vertices;
  double xm, ym, zm, xp, yp, zp;
#endif
  CollArray = new TCollection*[LEVELS+1];
  coll =  Domain->GetCollection(It_Finest, 0);
  CollArray[0] = coll;

  // refine grids
  for(i=0;i<mg_level;i++)
  {
    if(i && (i<LEVELS))
    {
      OutPut("Refine " << i <<endl);
      Domain->RegRefineAll();
      coll =  Domain->GetCollection(It_Finest, 0);
#ifdef __CALOTTE__
      MakeInAndOutflowCircular(coll);
#endif
      CollArray[i] = coll;
    }
    else
    {
      // multiple discretization multilevel method
      OutPut("mdml method " << i <<endl);
      coll =  Domain->GetCollection(It_Finest, 0);
      CollArray[i] = coll;
    }
#ifdef __WALL_MOUNTED_CUBE__
    coll = Domain->GetCollection(It_Finest, 0);
    N_CellsOld = coll->GetN_Cells();
    N_CellsNew = 0;
    for(j=0;j<N_CellsOld;j++)
    {
      cell = coll->GetCell(j);
      N_Vertices = cell->GetN_Vertices();
      // compute barycenter of mesh cell
      xm = ym = zm = 0.0;
      for(k=0;k<N_Vertices;k++)
      {
        cell->GetVertex(k)->GetCoords(xp, yp, zp);
        xm += xp;
        ym += yp;
        zm += zp;
      }
      xm /= N_Vertices;
      ym /= N_Vertices;
      zm /= N_Vertices;
      // check if mesh cell is in the cube
      if (!((xm > 0.3) && (xm < 0.4) && (ym > 0.3) && (ym < 0.4) && (zm<0.1)))
      {
        N_CellsNew++;
      }
    }                                             // endfor j
    OutPut("N_CellsOld: " << N_CellsOld << " N_CellsNew: " << N_CellsNew << endl);

    Cells = new TBaseCell*[N_CellsNew];
    N_CellsNew = 0;
    for(j=0;j<N_CellsOld;j++)
    {
      cell = coll->GetCell(j);
      N_Vertices = cell->GetN_Vertices();
      xm = ym = zm = 0.0;
      for(k=0;k<N_Vertices;k++)
      {
        cell->GetVertex(k)->GetCoords(xp, yp, zp);
        xm += xp;
        ym += yp;
        zm += zp;
      }
      xm /= N_Vertices;
      ym /= N_Vertices;
      zm /= N_Vertices;

      // check if mesh cell is in the cube
      if (!((xm > 0.3) && (xm < 0.4) && (ym > 0.3) && (ym < 0.4) && (zm<0.1)))
      {
        Cells[N_CellsNew] = cell;
        N_CellsNew++;
      }
      else
        OutPut(xm << " " << ym<< " " << zm <<endl);
    }                                             // endfor j

    delete coll;
    coll = new TCollection(N_CellsNew, Cells);
    CollArray[i] = coll;
#endif
    OutPut("cells " << coll->GetN_Cells()<< endl);

#ifdef __CHANNEL_TAU180__
    coll = Domain->GetCollection(It_Finest, 0);
    SetZCoordinates(coll, i);
    Domain->MakeBdParamsConsistent(coll);
    CheckZCoordinates(coll, i);
    delete coll;
#endif
#ifdef __CHANNEL_CAROLINA__
    coll = Domain->GetCollection(It_Finest, 0);
    SetZCoordinates(coll, i);
    Domain->MakeBdParamsConsistent(coll);
    CheckZCoordinates(coll, i);
    delete coll;
#endif

    if(TDatabase::ParamDB->WRITE_PS)
    {
      // write grid into an Postscript file
      os.seekp(std::ios::beg);
      os << PsBaseName << i << ".ps" << ends;
      //Domain->PS(os.str().c_str(),It_Finest,0);
      Domain->PS(os.str().c_str(),It_EQ,i);
    }
  }

  for (i=0;i<mg_level;i++)
  {
    N_cells = CollArray[i]->GetN_Cells();
    // reset clipboards to -1
    for(j=0;j<N_cells;j++)
    {
      cell=CollArray[i]->GetCell(j);
      l=cell->GetN_Joints();
      for(k=0;k<l;k++)
      {
        neigh=cell->GetJoint(k)->GetNeighbour(cell);
        if(neigh) neigh->SetClipBoard(-1);
      }
      cell->SetClipBoard(-1);
    }                                             // endfor j

    for(j=0;j<N_cells;j++)
      CollArray[i]->GetCell(j)->SetClipBoard(j);

    for(j=0;j<N_cells;j++)
    {
      cell=CollArray[i]->GetCell(j);
      l=cell->GetN_Joints();
      for(k=0;k<l;k++)
      {
        neigh=cell->GetJoint(k)->GetNeighbour(cell);
        if((neigh) && (neigh->GetClipBoard() == -1))
        {
          if (l==6)
            neigh->SetRefDesc(TDatabase::RefDescDB[Hexahedron]);
          else
            neigh->SetRefDesc(TDatabase::RefDescDB[Tetrahedron]);
        }
      }
    }                                             // endfor j

  }                                               // endfor i

#ifdef  __MIXINGLAYERSLIP3D__
  // set periodic bc
  for(i=0;i<mg_level;i++)
  {
    /*if (i==mg_level-1)
    coll=Domain->GetCollection(It_Finest, 0);
    else
    coll=Domain->GetCollection(It_EQ, i);*/
    coll = CollArray[i];
    SetPeriodicFaceJoints(coll);
  }
  TDatabase::ParamDB->COMPUTE_VORTICITY_DIVERGENCE = TRUE;
  if (TDatabase::ParamDB->P9== -1)
  {
    mixing_layer_galerkin = 1;
    smoothing_depth = 0;
    OutPut("Smoothing depth for filtered vorticity is "<< smoothing_depth
      << " levels" << endl);
  }
#endif
#ifdef  __CHANNEL_TOBIAS__
  // set periodic bc
  for(i=0;i<mg_level;i++)
  {
    /*if (i==mg_level-1)
    coll=Domain->GetCollection(It_Finest, 0);
    else
    coll=Domain->GetCollection(It_EQ, i);*/
    coll = CollArray[i];
    SetPeriodicFaceJoints(coll);
  }
#endif

#ifdef  __WINDTUNNEL__
  ReadExperimentalBoundaryConditions();
#endif

  t3 = GetTime();
  total_time = t3 - total_time;
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
      OutPut(LEVELS-1 << "              *******" << endl);
    }
    OutPut("******           MULTIGRID LEVEL ");
    OutPut(i << "              *******" << endl);
    OutPut("*******************************************************" << endl);
    OutPut("memory before: " << setw(10) << GetMemory() << endl);

    //if(i && (i<LEVELS)) Domain->RegRefineAll();
    /*if (i==mg_level-1)
      coll=Domain->GetCollection(It_Finest, 0);
    else
      coll=Domain->GetCollection(It_EQ, i+TDatabase::ParamDB->SC_COARSEST_LEVEL_SADDLE);
    if ((i==mg_level-2)&& (mg_type==2))
    coll=Domain->GetCollection(It_Finest, 0);*/
    j = i;
    if ((mg_type==1)&&(i == mg_level-1))
      j = i-1;
    if ((mg_type==2)&&(i == mg_level-2))
      j = i-1;
    if ((mg_type==2)&&(i == mg_level-1))
      j = i-2;
    coll = CollArray[j];

    TDatabase::ParamDB->INTERNAL_LEVEL = 0;
    if (i== mg_level-1)
      TDatabase::ParamDB->INTERNAL_LEVEL = 1;

    // get spaces for low order disc on finest geo grid
    if (((mg_type==1)&&(i<mg_level-1))||((mg_type==2)&&(i<mg_level-2)))
    {
      velocity_space = new TFESpace3D(coll, NameString, UString, BoundCondition,
        Non_USpace, 1);
      pressure_space = new TFESpace3D(coll, NameString, PString, BoundCondition,
        DiscP_PSpace, 0);
      velocity_space_code = -1;
      pressure_space_code = 0;
      order = -1;
      convolution_space = new TFESpace3D(coll, NameString, UString, BoundConditionAuxProblem,
        Non_USpace, 1);
    }
    // get spaces of high order disc on finest geo grid
    else
    {
      if ((mg_type==2)&&(i<mg_level-1))
      {
        GetVelocityAndPressureSpaceLow3D(coll,BoundCondition,
          velocity_space, &velocity_space_code,
          pressure_space, &pressure_space_code,
          TDatabase::ParamDB->VELOCITY_SPACE,
          TDatabase::ParamDB->PRESSURE_SPACE);
      }
      else
      {
        GetVelocityAndPressureSpace3D(coll,BoundCondition,
          velocity_space,
          pressure_space, &pressure_space_code,
          TDatabase::ParamDB->VELOCITY_SPACE,
          TDatabase::ParamDB->PRESSURE_SPACE);

        velocity_space_code =   TDatabase::ParamDB->VELOCITY_SPACE;
        TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
        GetVelocityAndPressureSpace3D(coll,BoundConditionAuxProblem,
          convolution_space,
          pressure_space, &pressure_space_code,
          TDatabase::ParamDB->VELOCITY_SPACE,
          TDatabase::ParamDB->PRESSURE_SPACE);

        OutPut("convolution space has same order as velocity space" << endl);
      }
    }

    if ((TDatabase::ParamDB->COMPUTE_VORTICITY_DIVERGENCE)&&(i==mg_level-1))
    {
      if (mixing_layer_galerkin)
      {
        OutPut("vorticity space is velocity space" << endl);
        vorticity_space = velocity_space;
      }
      else
      {
        OutPut("vorticity space is of first order" << endl);
        vorticity_space = new TFESpace3D(coll, NameString, VortString, BoundCondition, 1);
      }
      N_Vort = vorticity_space->GetN_DegreesOfFreedom();
    }

    // build fespace hierarchy
    // set values and pointers for fe space
    USpaces[i] = velocity_space;
    PSpaces[i] = pressure_space;
    N_U = velocity_space->GetN_DegreesOfFreedom();
    N_P = pressure_space->GetN_DegreesOfFreedom();
    N_Uarray[i] = N_U;
    N_Parray[i] = N_P;
    N_Active = velocity_space->GetActiveBound();

#ifdef __CHANNEL_TAU180__
    // compute correspondence of dof to coordinates
    if (i==mg_level-1)
    {
      x_dof = new double[3*N_U];
      y_dof = x_dof+N_U;
      z_dof = y_dof+N_U;
      for (j=0;j<N_U;j++)
        x_dof[j] = -4711;
      GetCoordinatesOfDof(coll, velocity_space, x_dof, y_dof, z_dof,
        &N_z_layers, coord_z_layers);
      OutPut("coordinate layers in z directions "<< N_z_layers<< endl);
      for(j=0;j<N_z_layers;j++)
      {
        OutPut(j << " " <<  coord_z_layers[j] << endl);
      }
      // allocate vector for mean velocities
      mean_velocity = new double[4*N_z_layers];
      memset(mean_velocity,0,4*N_z_layers*SizeOfDouble);
      mean_velocity_u2 = mean_velocity + N_z_layers;
      mean_velocity_u3 = mean_velocity_u2 + N_z_layers;
      dmean_velocity = mean_velocity_u3 + N_z_layers;
      // allocate vector for rms velocity
      rms_velocity = new double[6*N_z_layers];
      memset(rms_velocity,0,6*N_z_layers*SizeOfDouble);
      rms_velocity2 = rms_velocity + N_z_layers;
      rms_velocity3 = rms_velocity2 + N_z_layers;
      rms_velocity1_type1 =  rms_velocity3 + N_z_layers;
      rms_velocity2_type1 =  rms_velocity1_type1 + N_z_layers;
      rms_velocity3_type1 =  rms_velocity2_type1 + N_z_layers;
      // allocate vectors for the Reynolds stress tensor
      R_xx = new double[12*N_z_layers];
      memset(R_xx,0,12*N_z_layers*SizeOfDouble);
      R_xy = R_xx + N_z_layers;
      R_xz = R_xy + N_z_layers;
      R_yy = R_xz + N_z_layers;
      R_yz = R_yy + N_z_layers;
      R_zz = R_yz + N_z_layers;
      A_M_xx = R_zz +  N_z_layers;
      A_M_xy = A_M_xx + N_z_layers;
      A_M_xz = A_M_xy + N_z_layers;
      A_M_yy = A_M_xz + N_z_layers;
      A_M_yz = A_M_yy + N_z_layers;
      A_M_zz = A_M_yz + N_z_layers;
      // allocate vector for averaged gradient
      u1x = new double[9*N_U];
      if ((TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION)||
        (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION_EXPL))
        projection_u1x = new double[7*N_U];
      else
        projection_u1x = NULL;
    }
#endif
    // build matrices
    structureB = new TStructure3D(pressure_space, velocity_space);
    structureBT = new TStructure3D(velocity_space, pressure_space);
    sqstructureA = new TSquareStructure3D(velocity_space);
    sqstructureA->Sort();

    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        matrixB1 = new TMatrix3D(structureB);
        matrixB2 = new TMatrix3D(structureB);
        matrixB3 = new TMatrix3D(structureB);

        MatricesB1[i] = matrixB1;
        MatricesB2[i] = matrixB2;
        MatricesB3[i] = matrixB3;

        sqmatrixA = new TSquareMatrix3D(sqstructureA);
        MatricesA[i] = sqmatrixA;

        sqmatrixM = new TSquareMatrix3D(sqstructureA);
        MatricesM[i] = sqmatrixM;
        break;

      case 2:
        matrixB1 = new TMatrix3D(structureB);
        matrixB2 = new TMatrix3D(structureB);
        matrixB3 = new TMatrix3D(structureB);
        matrixB1T = new TMatrix3D(structureBT);
        matrixB2T = new TMatrix3D(structureBT);
        matrixB3T = new TMatrix3D(structureBT);

        MatricesB1[i] = matrixB1;
        MatricesB2[i] = matrixB2;
        MatricesB3[i] = matrixB3;
        MatricesB1T[i] = matrixB1T;
        MatricesB2T[i] = matrixB2T;
        MatricesB3T[i] = matrixB3T;

        sqmatrixA = new TSquareMatrix3D(sqstructureA);
        MatricesA[i] = sqmatrixA;

        sqmatrixM = new TSquareMatrix3D(sqstructureA);
        MatricesM[i] = sqmatrixM;
        break;

      case 3:
        matrixB1 = new TMatrix3D(structureB);
        matrixB2 = new TMatrix3D(structureB);
        matrixB3 = new TMatrix3D(structureB);

        MatricesB1[i] = matrixB1;
        MatricesB2[i] = matrixB2;
        MatricesB3[i] = matrixB3;

        sqmatrixA11 = new TSquareMatrix3D(sqstructureA);
        sqmatrixA12 = new TSquareMatrix3D(sqstructureA);
        sqmatrixA13 = new TSquareMatrix3D(sqstructureA);
        sqmatrixA21 = new TSquareMatrix3D(sqstructureA);
        sqmatrixA22 = new TSquareMatrix3D(sqstructureA);
        sqmatrixA23 = new TSquareMatrix3D(sqstructureA);
        sqmatrixA31 = new TSquareMatrix3D(sqstructureA);
        sqmatrixA32 = new TSquareMatrix3D(sqstructureA);
        sqmatrixA33 = new TSquareMatrix3D(sqstructureA);

        MatricesA11[i] = sqmatrixA11;
        MatricesA12[i] = sqmatrixA12;
        MatricesA13[i] = sqmatrixA13;
        MatricesA21[i] = sqmatrixA21;
        MatricesA22[i] = sqmatrixA22;
        MatricesA23[i] = sqmatrixA23;
        MatricesA31[i] = sqmatrixA31;
        MatricesA32[i] = sqmatrixA32;
        MatricesA33[i] = sqmatrixA33;

        sqmatrixM11 = new TSquareMatrix3D(sqstructureA);
        sqmatrixM12 = new TSquareMatrix3D(sqstructureA);
        sqmatrixM13 = new TSquareMatrix3D(sqstructureA);
        sqmatrixM21 = new TSquareMatrix3D(sqstructureA);
        sqmatrixM22 = new TSquareMatrix3D(sqstructureA);
        sqmatrixM23 = new TSquareMatrix3D(sqstructureA);
        sqmatrixM31 = new TSquareMatrix3D(sqstructureA);
        sqmatrixM32 = new TSquareMatrix3D(sqstructureA);
        sqmatrixM33 = new TSquareMatrix3D(sqstructureA);

        MatricesM11[i] = sqmatrixM11;
        MatricesM12[i] = sqmatrixM12;
        MatricesM13[i] = sqmatrixM13;
        MatricesM21[i] = sqmatrixM21;
        MatricesM22[i] = sqmatrixM22;
        MatricesM23[i] = sqmatrixM23;
        MatricesM31[i] = sqmatrixM31;
        MatricesM32[i] = sqmatrixM32;
        MatricesM33[i] = sqmatrixM33;
        break;

      case 4:
      case 14:
        matrixB1 = new TMatrix3D(structureB);
        matrixB2 = new TMatrix3D(structureB);
        matrixB3 = new TMatrix3D(structureB);
        matrixB1T = new TMatrix3D(structureBT);
        matrixB2T = new TMatrix3D(structureBT);
        matrixB3T = new TMatrix3D(structureBT);

        MatricesB1[i] = matrixB1;
        MatricesB2[i] = matrixB2;
        MatricesB3[i] = matrixB3;
        MatricesB1T[i] = matrixB1T;
        MatricesB2T[i] = matrixB2T;
        MatricesB3T[i] = matrixB3T;

        sqmatrixA11 = new TSquareMatrix3D(sqstructureA);
        sqmatrixA12 = new TSquareMatrix3D(sqstructureA);
        sqmatrixA13 = new TSquareMatrix3D(sqstructureA);
        sqmatrixA21 = new TSquareMatrix3D(sqstructureA);
        sqmatrixA22 = new TSquareMatrix3D(sqstructureA);
        sqmatrixA23 = new TSquareMatrix3D(sqstructureA);
        sqmatrixA31 = new TSquareMatrix3D(sqstructureA);
        sqmatrixA32 = new TSquareMatrix3D(sqstructureA);
        sqmatrixA33 = new TSquareMatrix3D(sqstructureA);

        MatricesA11[i] = sqmatrixA11;
        MatricesA12[i] = sqmatrixA12;
        MatricesA13[i] = sqmatrixA13;
        MatricesA21[i] = sqmatrixA21;
        MatricesA22[i] = sqmatrixA22;
        MatricesA23[i] = sqmatrixA23;
        MatricesA31[i] = sqmatrixA31;
        MatricesA32[i] = sqmatrixA32;
        MatricesA33[i] = sqmatrixA33;

        sqmatrixM11 = new TSquareMatrix3D(sqstructureA);
        sqmatrixM12 = new TSquareMatrix3D(sqstructureA);
        sqmatrixM13 = new TSquareMatrix3D(sqstructureA);
        sqmatrixM21 = new TSquareMatrix3D(sqstructureA);
        sqmatrixM22 = new TSquareMatrix3D(sqstructureA);
        sqmatrixM23 = new TSquareMatrix3D(sqstructureA);
        sqmatrixM31 = new TSquareMatrix3D(sqstructureA);
        sqmatrixM32 = new TSquareMatrix3D(sqstructureA);
        sqmatrixM33 = new TSquareMatrix3D(sqstructureA);

        MatricesM11[i] = sqmatrixM11;
        MatricesM12[i] = sqmatrixM12;
        MatricesM13[i] = sqmatrixM13;
        MatricesM21[i] = sqmatrixM21;
        MatricesM22[i] = sqmatrixM22;
        MatricesM23[i] = sqmatrixM23;
        MatricesM31[i] = sqmatrixM31;
        MatricesM32[i] = sqmatrixM32;
        MatricesM33[i] = sqmatrixM33;

        if (TDatabase::ParamDB->NSTYPE == 14)
        {
          sqstructureC = new TSquareStructure3D(pressure_space);
          sqstructureC->Sort();
          sqmatrixC = new TSquareMatrix3D(sqstructureC);
          sqmatrixK = new TSquareMatrix3D(sqstructureA);
          MatricesC[i] = sqmatrixC;
          MatricesK[i] = sqmatrixK;
        }
        if (TDatabase::ParamDB->DISCTYPE == SDFEM)
        {
          sqmatrixK = new TSquareMatrix3D(sqstructureA);
          MatricesK[i] = sqmatrixK;
          sqmatrixS11 = new TSquareMatrix3D(sqstructureA);
          sqmatrixS12 = new TSquareMatrix3D(sqstructureA);
          sqmatrixS13 = new TSquareMatrix3D(sqstructureA);
          sqmatrixS21 = new TSquareMatrix3D(sqstructureA);
          sqmatrixS22 = new TSquareMatrix3D(sqstructureA);
          sqmatrixS23 = new TSquareMatrix3D(sqstructureA);
          sqmatrixS31 = new TSquareMatrix3D(sqstructureA);
          sqmatrixS32 = new TSquareMatrix3D(sqstructureA);
          sqmatrixS33 = new TSquareMatrix3D(sqstructureA);

          MatricesS11[i] = sqmatrixS11;
          MatricesS12[i] = sqmatrixS12;
          MatricesS13[i] = sqmatrixS13;
          MatricesS21[i] = sqmatrixS21;
          MatricesS22[i] = sqmatrixS22;
          MatricesS23[i] = sqmatrixS23;
          MatricesS31[i] = sqmatrixS31;
          MatricesS32[i] = sqmatrixS32;
          MatricesS33[i] = sqmatrixS33;
        }
        break;
    }
    // allocate matrix for auxiliary problem in Galdi/Layton model
    if (TDatabase::ParamDB->DISCTYPE==GL00_AUX_PROBLEM)
    {
      sqstructureC = new TSquareStructure3D(convolution_space);
      sqstructureC->Sort();
      sqmatrixGL00AuxProblem = new  TSquareMatrix3D(sqstructureC);
      MatricesGL00AuxProblem[i] =  sqmatrixGL00AuxProblem;
    }

    N_Unknowns = 3*N_U + N_P;

    OutPut("dof velocity : "<< setw(10) << 3* N_U << endl);
    OutPut("dof pressure : "<< setw(10) << N_P << endl);
    OutPut("dof all      : "<<  setw(10) << N_Unknowns  << endl);

    if ((TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION)||
      (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION_EXPL))
    {
      switch (TDatabase::ParamDB->VMS_LARGE_VELOCITY_SPACE)
      {
        case 0:
          projection_space = new TFESpace3D(coll, NameString, PString, BoundCondition,
            DiscP_PSpace, 0);
          break;
        case 1:
          projection_space = new TFESpace3D(coll, NameString, PString, BoundCondition,
            DiscP_PSpace, 1);
          break;
        case 2:
          projection_space = new TFESpace3D(coll, NameString, PString, BoundCondition,
            DiscP_PSpace, 2);
          break;
        case 3:
          projection_space = new TFESpace3D(coll, NameString, PString, BoundCondition,
            DiscP_PSpace, 3);
          break;
        case 17:
          fes = new FE3D[coll->GetN_Cells()];
          projection_space = new TFESpace3D(coll, NameString, PString, BoundCondition,
            DiscP_PSpace, 0);
          /*  for (j=0;j<coll->GetN_Cells();j++)
              fes[j] = C_Q0_3D_H_A;
            //AdaptProjectionSpace(projection_space, size_small_scales, fes);
            projection_space = new TFESpace3D(coll, NameString, PString, BoundCondition,
              fes);
          OutPut("dofs " << projection_space->GetN_DegreesOfFreedom()<<endl); */
          break;
        default :
          OutPut("space number not implemented !" << endl);
          exit(4711);
      }

      ProjectionSpaces[i] = projection_space;
      sqstructureL = new TSquareStructure3D(projection_space);
      sqstructureL->Sort();
      structure_tilde_G = new TStructure3D(velocity_space, projection_space);
      structure_G = new TStructure3D(projection_space, velocity_space);
      sqmatrixL = new TSquareMatrix3D(sqstructureL);
      MatricesL[i] = sqmatrixL;
      matrix_tilde_G11 = new TMatrix3D(structure_tilde_G);
      Matrices_tilde_G11[i] = matrix_tilde_G11;
      matrix_tilde_G22 = new TMatrix3D(structure_tilde_G);
      Matrices_tilde_G22[i] = matrix_tilde_G22;
      matrix_tilde_G33 = new TMatrix3D(structure_tilde_G);
      Matrices_tilde_G33[i] = matrix_tilde_G33;
      matrix_G11 = new TMatrix3D(structure_G);
      Matrices_G11[i] = matrix_G11;
      matrix_G22 = new TMatrix3D(structure_G);
      Matrices_G22[i] = matrix_G22;
      matrix_G33 = new TMatrix3D(structure_G);
      Matrices_G33[i] = matrix_G33;
      N_L = projection_space->GetN_DegreesOfFreedom();
      OutPut("dof projection : " << setw(10) << N_L << endl);

      if (i == mg_level -1)
      {
        // fe function for the large scales, only needed on
        // finest level
        vms_projection = new double[6*N_L];
        memset(vms_projection,0,6*N_L*SizeOfDouble);
        vms_projection_fe = new TFEVectFunct3D(projection_space, PString, PString, vms_projection, N_L, 6);
        vms_proj_11 = vms_projection_fe->GetComponent(0);
        vms_proj_12 = vms_projection_fe->GetComponent(1);
        vms_proj_13 = vms_projection_fe->GetComponent(2);
        vms_proj_22 = vms_projection_fe->GetComponent(3);
        vms_proj_23 = vms_projection_fe->GetComponent(4);
        vms_proj_33 = vms_projection_fe->GetComponent(5);
        size_small_scales = new double[coll->GetN_Cells()];
        memset(size_small_scales,0,coll->GetN_Cells()*SizeOfDouble);
        label_space = new double[coll->GetN_Cells()];
        memset(label_space,0,coll->GetN_Cells()*SizeOfDouble);

        size_small_scales_fesp = new TFESpace3D(coll, NameString, SmallScaleString, BoundCondition,
          DiscP_PSpace, 0);
        size_small_scales_fefct = new TFEFunction3D(size_small_scales_fesp, SmallScaleString,
          SmallScaleString,
          size_small_scales, coll->GetN_Cells());

        label_space_fesp = new TFESpace3D(coll, NameString, LabelSpaceString, BoundCondition,
          DiscP_PSpace, 0);
        label_space_fefct = new TFEFunction3D(label_space_fesp, LabelSpaceString,
          LabelSpaceString,
          label_space, coll->GetN_Cells());
      }
    }

    rhs = new double[N_Unknowns];
    memset(rhs, 0, N_Unknowns*SizeOfDouble);
    RhsArray[i] = rhs;
    sol = new double[N_Unknowns];
    oldsol = new double[N_Unknowns];
    memset(sol, 0, N_Unknowns*SizeOfDouble);
    memset(oldsol, 0, N_Unknowns*SizeOfDouble);
    if ((i==mg_level-1)&&(TDatabase::TimeDB->EXTRAPOLATE_VELOCITY))
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
            matrixB3, structureBT, B, sol, n_aux,  alpha,
            velocity_space_code, pressure_space_code,NULL,downwind[i]);
          break;

        case 2:
          MGLevel = new TNSE_MGLevel2(i, sqmatrixM, matrixB1, matrixB2,
            matrixB3,  matrixB1T, matrixB2T, matrixB3T,
            B, sol, n_aux, alpha,
            velocity_space_code, pressure_space_code,NULL,downwind[i]);
          break;

        case 3:
          MGLevel = new TNSE_MGLevel3(i,
            sqmatrixM11, sqmatrixM12, sqmatrixM13,
            sqmatrixM21, sqmatrixM22, sqmatrixM23,
            sqmatrixM31, sqmatrixM32, sqmatrixM33,
            matrixB1, matrixB2, matrixB3,
            structureBT, B, sol, n_aux, alpha,
            velocity_space_code,  pressure_space_code,NULL,downwind[i]);
          break;

        case 4:
          MGLevel = new TNSE_MGLevel4(i,
            sqmatrixM11, sqmatrixM12, sqmatrixM13,
            sqmatrixM21, sqmatrixM22, sqmatrixM23,
            sqmatrixM31, sqmatrixM32, sqmatrixM33,
            matrixB1, matrixB2, matrixB3,
            matrixB1T, matrixB2T, matrixB3T,
            B, sol, n_aux, alpha,
            velocity_space_code, pressure_space_code,coll,downwind[i]);
          break;

        case 14:
          MGLevel = new TNSE_MGLevel14(i,
            sqmatrixM11, sqmatrixM12, sqmatrixM13,
            sqmatrixM21, sqmatrixM22, sqmatrixM23,
            sqmatrixM31, sqmatrixM32, sqmatrixM33,
            sqmatrixC,
            matrixB1, matrixB2, matrixB3,
            matrixB1T, matrixB2T, matrixB3T,
            B, sol, n_aux, alpha,
            velocity_space_code, pressure_space_code,coll,downwind[i]);
          break;
      }                                           // end switch(NSTYPE)
      MG->AddLevel(MGLevel);
    }

    u = new TFEVectFunct3D(velocity_space, UString, UString, sol, N_U, 3);
    u1 = u->GetComponent(0);
    u2 = u->GetComponent(1);
    u3 = u->GetComponent(2);
    p = new TFEFunction3D(pressure_space, PString, PString, sol+3*N_U, N_P);

    U1Array[i] = u1;
    U2Array[i] = u2;
    U3Array[i] = u3;
    PArray[i] = p;
    UArray[i] = u;

    u1->Interpolate(InitialU1);
    u2->Interpolate(InitialU2);
    u3->Interpolate(InitialU3);
    p->Interpolate(InitialP);

    RHSs[0] = rhs;
    RHSs[1] = rhs + N_U;
    RHSs[2] = rhs + 2*N_U;
    RHSs[3] = rhs + 3*N_U;
    memset(rhs, 0, N_Unknowns*SizeOfDouble);

    if (TDatabase::ParamDB->DISCTYPE==CLASSICAL_LES)
    {
      if (i==mg_level-1)
      {
        LESModelRhs =  new double[3*N_U];
        memset(LESModelRhs,0,3*N_U*SizeOfDouble);
      }
    }
    if (TDatabase::ParamDB->DISCTYPE==GL00_CONVOLUTION)
    {
      duConvSpaces[i] = convolution_space;
      // define vector fe function for convolution
      N_UConv = convolution_space->GetN_DegreesOfFreedom();
      // allocate memory for values of convolved function
      u_conv = new double[6*N_UConv];
      // initialize u_conv to 0
      memset(u_conv,0,6*N_UConv*SizeOfDouble);
      // allocate fe vector function for convolution
      duConv = new TFEVectFunct3D(convolution_space, UConvString, UConvString,
        u_conv, N_UConv, 6);
      // array for duConv for all levels
      duConvArray[i] = duConv;

      // copy the vector fe function to a fe function (only pointers)
      du11Conv = duConv->GetComponent(0);
      du12Conv = duConv->GetComponent(1);
      du13Conv = duConv->GetComponent(2);
      du22Conv = duConv->GetComponent(3);
      du23Conv = duConv->GetComponent(4);
      du33Conv = duConv->GetComponent(5);
      du11ConvArray[i] = du11Conv;
      du12ConvArray[i] = du12Conv;
      du13ConvArray[i] = du13Conv;
      du22ConvArray[i] = du22Conv;
      du23ConvArray[i] = du23Conv;
      du33ConvArray[i] = du33Conv;
      if (i==mg_level-1)
      {
        LESModelRhs =  new double[3*N_U];
        memset(LESModelRhs,0,3*N_U*SizeOfDouble);
      }
    }

    if  (TDatabase::ParamDB->DISCTYPE==GL00_AUX_PROBLEM)
    {
      // allocate arrays for multigrid for auxiliary problem in
      // Galdi/Layton model
      solGL00AuxProblem =  new double[6*N_U];
      memset(solGL00AuxProblem,0,6*N_U*SizeOfDouble);
      rhsGL00AuxProblem =  new double[6*N_U];
      memset(rhsGL00AuxProblem,0,6*N_U*SizeOfDouble);
      if (i==mg_level-1)
      {
        LESModelRhs =  new double[3*N_U];
        memset(LESModelRhs,0,3*N_U*SizeOfDouble);
      }

      // build multigrid for auxiliary problem in Galdi/Layton model
      MGLevelGL00AuxProblem= new TMGLevel3D(i,sqmatrixGL00AuxProblem,
        rhsGL00AuxProblem, solGL00AuxProblem,
        n_aux, NULL);
      MGGL00AuxProblem->AddLevel(MGLevelGL00AuxProblem);

      // allocate fe vector function for solution of auxiliary problem
      GL00AuxProblemSol = new TFEVectFunct3D(convolution_space,
        AuxProblemString,
        AuxProblemString,
        solGL00AuxProblem, N_U, 6);
      // array for duConv for all levels
      GL00AuxProblemSolArray[i] = GL00AuxProblemSol;

      // copy the vector fe function to a fe function (only pointers)
      GL00AuxProblemSol11 = GL00AuxProblemSol->GetComponent(0);
      GL00AuxProblemSol12 = GL00AuxProblemSol->GetComponent(1);
      GL00AuxProblemSol13 = GL00AuxProblemSol->GetComponent(2);
      GL00AuxProblemSol22 = GL00AuxProblemSol->GetComponent(3);
      GL00AuxProblemSol23 = GL00AuxProblemSol->GetComponent(4);
      GL00AuxProblemSol33 = GL00AuxProblemSol->GetComponent(5);
      GL00AuxProblemSol11Array[i] = GL00AuxProblemSol11;
      GL00AuxProblemSol12Array[i] = GL00AuxProblemSol12;
      GL00AuxProblemSol13Array[i] = GL00AuxProblemSol13;
      GL00AuxProblemSol22Array[i] = GL00AuxProblemSol22;
      GL00AuxProblemSol23Array[i] = GL00AuxProblemSol23;
      GL00AuxProblemSol33Array[i] = GL00AuxProblemSol33;
    }

    // turbulent viscosity involving the convolution of the solution
    if ((TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE==4)||
      (TDatabase::ParamDB->CONVOLUTE_SOLUTION))
    {
      uConvSpaces[i] = convolution_space;
      // define vector fe function for convolution
      // allocate memory for values of convoluted function
      u_uConv = new double[3*N_U];
      // initialize u_conv to 0
      memset(u_uConv,0,3*N_U*SizeOfDouble);
      // allocate fe vector function for convolution
      uConv = new TFEVectFunct3D(convolution_space, UConvString, UConvString,
        u_uConv, N_U, 3);
      // array for uConv for all levels
      uConvArray[i] = uConv;

      // copy the vector fe function to a fe function (only pointers)
      u1Conv = uConv->GetComponent(0);
      u2Conv = uConv->GetComponent(1);
      u3Conv = uConv->GetComponent(2);
      u1ConvArray[i] = u1Conv;
      u2ConvArray[i] = u2Conv;
      u3ConvArray[i] = u3Conv;

      // compute matrix for auxiliary problem if not yet done
      if (TDatabase::ParamDB->DISCTYPE!=GL00_AUX_PROBLEM)
      {
        if (i==mg_level-1)
        {
          sqstructureC = new TSquareStructure3D(convolution_space);
          sqstructureC->Sort();
          sqmatrixGL00AuxProblem = new  TSquareMatrix3D(sqstructureC);
          MatricesGL00AuxProblem[i] =  sqmatrixGL00AuxProblem;
          rhsGL00AuxProblem =  new double[3*N_U];
          memset(rhsGL00AuxProblem,0,3*N_U*SizeOfDouble);
          // assemble matrix
          DiscreteForm = DiscreteFormMatrixAuxProblemU;
          fesp[0] = convolution_space;
          fesp[1] = velocity_space;

          fefct[0] = u1;
          fefct[1] = u2;
          fefct[2] = u3;

          // assemble matrix
          N_FESpaces = 2;
          N_Rhs = 0;
          N_SquareMatrices = 1;
          N_RectMatrices = 0;

          SQMATRICES[0] = MatricesGL00AuxProblem[mg_level-1];
          SQMATRICES[0]->Reset();
          aux =  new TAuxParam3D(TimeNSN_FESpacesVelo, TimeNSN_FctVelo,
            TimeNSN_ParamFctVelo,
            TimeNSN_FEValuesVelo,
            fesp+1, fefct,
            TimeNSFctVelo,
            TimeNSFEFctIndexVelo, TimeNSFEMultiIndexVelo,
            TimeNSN_ParamsVelo, TimeNSBeginParamVelo);

          Assemble3D(N_FESpaces, fesp,
            N_SquareMatrices, SQMATRICES,
            N_RectMatrices, MATRICES,
            N_Rhs, RHSs, ferhs,
            DiscreteForm,
            BoundaryConditionsAuxProblem,
            BoundValuesAuxProblem,
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
      U1old = new TFEFunction3D(velocity_space, UString,  UString, former_sol, N_U);
      U2old = new TFEFunction3D(velocity_space, UString,  UString, former_sol+N_U, N_U);
      PreparePressureAtCylinder(coll, press_cyl, N_press_cyl);
      PrepareCenterlineVelocities(coll, center_velo, N_center_velo);
      PrepareVelocityAtCylinder(coll, cyl_velo, N_cyl_velo);
    }
#endif

#ifdef __CHANNEL_OBSTACLE__
    if (i==mg_level-1)
    {
      former_sol =  new double [2*N_U];
      memset(former_sol, 0, 2*N_U*SizeOfDouble);
      U1old = new TFEFunction3D(velocity_space, UString,  UString, former_sol, N_U);
      U2old = new TFEFunction3D(velocity_space, UString,  UString, former_sol+N_U, N_U);
    }
#endif

    // set discrete forms
    if (((mg_type==1)&&(i<mg_level-1))||((mg_type==2)&&(i<mg_level-2)))
    {
      DiscreteForm = DiscreteFormUpwindNC;
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
      case VMS_RFB_EXPL:
      case VMS_RFB_EXPL_COUPLED:
        DiscreteForm = DiscreteFormSmagorinsky;
        CurrentDiscType =  SMAGORINSKY;
        break;
      case  CLASSICAL_LES:
        DiscreteForm = DiscreteFormClassicalLES;
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
        DiscreteForm = DiscreteFormVMS_Projection;
        CurrentDiscType =  VMS_PROJECTION;
        break;
      case VMS_PROJECTION_EXPL:
        DiscreteForm = DiscreteFormVMS_Projection;
        CurrentDiscType =  VMS_PROJECTION_EXPL;
        break;
      case SDFEM:
        DiscreteForm = DiscreteFormVMS_SUPG;
        CurrentDiscType =  SDFEM;
        break;
      default:
        OutPut("Unknown DISCTYPE " << TDatabase::ParamDB->DISCTYPE << endl);
        exit(1);
    }
    if (DiscreteForm==NULL)
    {
      OutPut("DiscreteForm not implemented !!!"<< endl);
      exit(4711);
    }
    // parameters which are the minimum for all NSTYPEs
    N_Rhs = 3;
    N_FESpaces = 2;

    // set matrices
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        SQMATRICES[0] = MatricesA[i];
        SQMATRICES[1] = MatricesM[i];
        MATRICES[0] = MatricesB1[i];
        MATRICES[1] = MatricesB2[i];
        MATRICES[2] = MatricesB3[i];

        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        MATRICES[0]->Reset();
        MATRICES[1]->Reset();
        MATRICES[2]->Reset();

        N_SquareMatrices = 2;
        N_RectMatrices = 3;

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
        MATRICES[2] = MatricesB3[i];
        MATRICES[3] = MatricesB1T[i];
        MATRICES[4] = MatricesB2T[i];
        MATRICES[5] = MatricesB3T[i];

        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        MATRICES[0]->Reset();
        MATRICES[1]->Reset();
        MATRICES[2]->Reset();
        MATRICES[3]->Reset();
        MATRICES[4]->Reset();
        MATRICES[5]->Reset();

        N_SquareMatrices = 2;
        N_RectMatrices = 6;

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
        SQMATRICES[2] = MatricesA13[i];
        SQMATRICES[3] = MatricesA21[i];
        SQMATRICES[4] = MatricesA22[i];
        SQMATRICES[5] = MatricesA23[i];
        SQMATRICES[6] = MatricesA31[i];
        SQMATRICES[7] = MatricesA32[i];
        SQMATRICES[8] = MatricesA33[i];
        SQMATRICES[9] = MatricesM11[i];
        SQMATRICES[10] = MatricesM22[i];
        SQMATRICES[11] = MatricesM33[i];
        MATRICES[0] = MatricesB1[i];
        MATRICES[1] = MatricesB2[i];
        MATRICES[2] = MatricesB3[i];

        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        SQMATRICES[2]->Reset();
        SQMATRICES[3]->Reset();
        SQMATRICES[4]->Reset();
        SQMATRICES[5]->Reset();
        SQMATRICES[6]->Reset();
        SQMATRICES[7]->Reset();
        SQMATRICES[8]->Reset();
        SQMATRICES[9]->Reset();
        SQMATRICES[10]->Reset();
        SQMATRICES[11]->Reset();
        MATRICES[0]->Reset();
        MATRICES[1]->Reset();
        MATRICES[2]->Reset();

        N_SquareMatrices = 12;
        N_RectMatrices = 3;

        if (CurrentDiscType ==  GL00_AUX_PROBLEM)
        {
          N_SquareMatrices = 13;
          SQMATRICES[12] = MatricesGL00AuxProblem[i];
          SQMATRICES[12]->Reset();
        }
        break;

      case 4:
      case 14:
        SQMATRICES[0] = MatricesA11[i];
        SQMATRICES[1] = MatricesA12[i];
        SQMATRICES[2] = MatricesA13[i];
        SQMATRICES[3] = MatricesA21[i];
        SQMATRICES[4] = MatricesA22[i];
        SQMATRICES[5] = MatricesA23[i];
        SQMATRICES[6] = MatricesA31[i];
        SQMATRICES[7] = MatricesA32[i];
        SQMATRICES[8] = MatricesA33[i];
        SQMATRICES[9] = MatricesM11[i];
        SQMATRICES[10] = MatricesM22[i];
        SQMATRICES[11] = MatricesM33[i];
        MATRICES[0] = MatricesB1[i];
        MATRICES[1] = MatricesB2[i];
        MATRICES[2] = MatricesB3[i];
        MATRICES[3] = MatricesB1T[i];
        MATRICES[4] = MatricesB2T[i];
        MATRICES[5] = MatricesB3T[i];

        SQMATRICES[0]->Reset();
        SQMATRICES[1]->Reset();
        SQMATRICES[2]->Reset();
        SQMATRICES[3]->Reset();
        SQMATRICES[4]->Reset();
        SQMATRICES[5]->Reset();
        SQMATRICES[6]->Reset();
        SQMATRICES[7]->Reset();
        SQMATRICES[8]->Reset();
        SQMATRICES[9]->Reset();
        SQMATRICES[10]->Reset();
        SQMATRICES[11]->Reset();
        MATRICES[0]->Reset();
        MATRICES[1]->Reset();
        MATRICES[2]->Reset();
        MATRICES[3]->Reset();
        MATRICES[4]->Reset();
        MATRICES[5]->Reset();

        N_SquareMatrices = 12;
        N_RectMatrices = 6;

        if (CurrentDiscType ==  GL00_AUX_PROBLEM)
        {
          N_SquareMatrices += 1;
          SQMATRICES[12] = MatricesGL00AuxProblem[i];
          SQMATRICES[12]->Reset();
        }
        if ((CurrentDiscType == VMS_PROJECTION)||
          (CurrentDiscType == VMS_PROJECTION_EXPL))
        {
          N_SquareMatrices += 1;
          SQMATRICES[N_SquareMatrices-1] =  MatricesL[i];
          SQMATRICES[N_SquareMatrices-1]->Reset();
          N_RectMatrices += 6;
          MATRICES[6] = Matrices_tilde_G11[i];
          MATRICES[7] = Matrices_tilde_G22[i];
          MATRICES[8] = Matrices_tilde_G33[i];
          MATRICES[9] = Matrices_G11[i];
          MATRICES[10] = Matrices_G22[i];
          MATRICES[11] = Matrices_G33[i];
          MATRICES[6]->Reset();
          MATRICES[7]->Reset();
          MATRICES[8]->Reset();
          MATRICES[9]->Reset();
          MATRICES[10]->Reset();
          MATRICES[11]->Reset();
        }
        if (CurrentDiscType == SDFEM)
        {
          N_SquareMatrices += 10;
          SQMATRICES[N_SquareMatrices-10] =  MatricesK[i];
          SQMATRICES[N_SquareMatrices-10]->Reset();
          SQMATRICES[N_SquareMatrices-9] =  MatricesS11[i];
          SQMATRICES[N_SquareMatrices-9]->Reset();
          SQMATRICES[N_SquareMatrices-8] =  MatricesS12[i];
          SQMATRICES[N_SquareMatrices-8]->Reset();
          SQMATRICES[N_SquareMatrices-7] =  MatricesS13[i];
          SQMATRICES[N_SquareMatrices-7]->Reset();
          SQMATRICES[N_SquareMatrices-6] =  MatricesS21[i];
          SQMATRICES[N_SquareMatrices-6]->Reset();
          SQMATRICES[N_SquareMatrices-5] =  MatricesS22[i];
          SQMATRICES[N_SquareMatrices-5]->Reset();
          SQMATRICES[N_SquareMatrices-4] =  MatricesS23[i];
          SQMATRICES[N_SquareMatrices-4]->Reset();
          SQMATRICES[N_SquareMatrices-3] =  MatricesS31[i];
          SQMATRICES[N_SquareMatrices-3]->Reset();
          SQMATRICES[N_SquareMatrices-2] =  MatricesS32[i];
          SQMATRICES[N_SquareMatrices-2]->Reset();
          SQMATRICES[N_SquareMatrices-1] =  MatricesS33[i];
          SQMATRICES[N_SquareMatrices-1]->Reset();
          TDatabase::TimeDB->CURRENTTIMESTEPLENGTH =
            TDatabase::TimeDB->TIMESTEPLENGTH;
        }
        break;
    }

    // prepare output
    if (i==mg_level-1)
    {
      if (TDatabase::ParamDB->COMPUTE_VORTICITY_DIVERGENCE)
      {
        if (mixing_layer_galerkin)
          // 0 - 2 : voricity
          // 3 - 5 : filtered vorticity
          N_vort = 6;
        else
          N_vort = 3;
        vorticity = new double[N_vort*N_Vort];
        Vorticity = new TFEVectFunct3D(vorticity_space, VorticityString,
          VorticityString, vorticity, N_Vort, N_vort);
        Vort_x = Vorticity->GetComponent(0);
        Vort_y = Vorticity->GetComponent(1);
        Vort_z = Vorticity->GetComponent(2);
        div = new double[N_Vort];
        Divergence = new TFEFunction3D(vorticity_space, DivergenceString,
          DivergenceString, div, N_Vort);
        if (mixing_layer_galerkin)
          sol_vort_tmp = new double[3*N_Vort];
      }
      if ((TDatabase::ParamDB->WRITE_GRAPE)||(TDatabase::ParamDB->WRITE_GMV)||
        (TDatabase::ParamDB->WRITE_VTK))
      {
        if (TDatabase::ParamDB->COMPUTE_VORTICITY_DIVERGENCE)
        {
          // output velocity, pressure, divergence, vorticity (z-component)
          Output = new TOutput3D(5, 4 , 1, 1, Domain);
          Output->AddFEVectFunct(u);
          Output->AddFEFunction(p);
          Output->AddFEFunction(Divergence);
          Output->AddFEFunction(Vort_z);
          // convolved vorticity (z-component)
          if (mixing_layer_galerkin)
            Output->AddFEFunction(Vorticity->GetComponent(5));
          os.seekp(std::ios::beg);
          Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());
        }
        else
        {
          // output only velocity and pressure
          Output = new TOutput3D(4, 2 , 1, 1, Domain);
          Output->AddFEVectFunct(u);
          Output->AddFEFunction(p);
          if ((TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION)||
            (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION_EXPL))
          {
            Output->AddFEFunction(size_small_scales_fefct);
            Output->AddFEFunction(label_space_fefct);
          }
          os.seekp(std::ios::beg);
          Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());
        }
      }
    }                                             // end of preparing output

    // read initial solution of finest level from grape file
    if ((i==mg_level-1)&&(TDatabase::ParamDB->READ_GRAPE_FILE))
    {
      // only velocity and pressure is read
      AuxFEFunctArray = new TFEFunction3D*[1];
      AuxFEFunctArray[0] = PArray[mg_level-1];
      AuxFEVectFunctArray = new TFEVectFunct3D*[1];
      AuxFEVectFunctArray[0] = UArray[mg_level-1];
      ReadGrapeFile3D(ReadGrapeBaseName, 1 , 1, AuxFEFunctArray,AuxFEVectFunctArray);
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
      //  sol[3*N_U+j] *= 1.5;
    }
    /* THIS IS FOR PROLONGATION FROM A COARSE GRID */
    // read initial solution of finest level from grape file
    /*if ((i==mg_level-1)&&(TDatabase::ParamDB->READ_DATA))
    {
      auxConv = new double[3*N_U];
        Prolongate(USpaces[i-1], USpaces[i], 3,
          UArray[i-1]->GetValues(),
          UArray[i]->GetValues(),
          auxConv);
    delete auxConv;
    }*/
    t1 = GetTime();
    // set rhs
    fesp[0] = velocity_space;
    fesp[1] = pressure_space;

    fefct[0] = u1;
    fefct[1] = u2;
    fefct[2] = u3;

    ferhs[0] = velocity_space;
    ferhs[1] = velocity_space;
    ferhs[2] = velocity_space;
    ferhs[3] = pressure_space;

    // Newton's method
    if (nonlinite==1)
      aux =  new TAuxParam3D(TimeNSN_FESpacesVelo_GradVelo,
        TimeNSN_FctVelo_GradVelo,
        TimeNSN_ParamFctVelo_GradVelo,
        TimeNSN_FEValuesVelo_GradVelo,
        fesp, fefct,
        TimeNSFctVelo_GradVelo,
        TimeNSFEFctIndexVelo_GradVelo,
        TimeNSFEMultiIndexVelo_GradVelo,
        TimeNSN_ParamsVelo_GradVelo,
        TimeNSBeginParamVelo_GradVelo);
    else
    {
      // 3 parameters are needed for assembling
      // which are u1_old, u2_old, norm of grad u_old
      switch(CurrentDiscType)
      {
        // turbulent viscosity must be computed
        case SMAGORINSKY:
        case CLASSICAL_LES:
        case GL00_CONVOLUTION:
        case GL00_AUX_PROBLEM:
        case SDFEM:
          N_FESpaces = 3;
          fesp[2] = convolution_space;
          fefct[3] = du11ConvArray[i];
          fefct[4] = du12ConvArray[i];
          fefct[5] = du13ConvArray[i];
          fefct[6] = du22ConvArray[i];
          fefct[7] = du23ConvArray[i];
          fefct[8] = du33ConvArray[i];

          aux =  new TAuxParam3D(TimeNSN_FESpacesVelo_GradVelo,
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
        case VMS_PROJECTION:
        case VMS_PROJECTION_EXPL:
          //fesp[2] = velocity_space;

          fefct[0] = u1;
          fefct[1] = u2;
          fefct[2] = u3;
          if (i== mg_level -1)
          {
            OutPut("fine"<<endl);
            fesp[2] = ProjectionSpaces[i];
            fesp[3] = label_space_fesp;
            N_FESpaces = 4;
            TDatabase::ParamDB->INTERNAL_LEVEL = 1;
            fefct[3] = vms_proj_11;
            fefct[4] = vms_proj_12;
            fefct[5] = vms_proj_13;
            fefct[6] = vms_proj_22;
            fefct[7] = vms_proj_23;
            fefct[8] = vms_proj_33;
            fefct[9] = label_space_fefct;
          }
          else
          {                                       // just dummies
            N_FESpaces = 3;
            fesp[2] = ProjectionSpaces[i];
            OutPut("coarse"<<endl);
            TDatabase::ParamDB->INTERNAL_LEVEL = 0;
            fefct[3] = u1;
            fefct[4] = u1;
            fefct[5] = u1;
            fefct[6] = u1;
            fefct[7] = u1;
            fefct[8] = u1;
            fefct[9] = u1;
          }

          aux =  new TAuxParam3D(TimeNSN_FESpacesVelo_GradVelo_LargeScale,
            TimeNSN_FctVelo_GradVelo_LargeScale,
            TimeNSN_ParamFctVelo_GradVelo_LargeScale,
            TimeNSN_FEValuesVelo_GradVelo_LargeScale,
            fesp, fefct,
            TimeNSFctVelo_GradVelo_LargeScale,
            TimeNSFEFctIndexVelo_GradVelo_LargeScale,
            TimeNSFEMultiIndexVelo_GradVelo_LargeScale,
            TimeNSN_ParamsVelo_GradVelo_LargeScale,
            TimeNSBeginParamVelo_GradVelo_LargeScale);
          break;
          // no parameter necessary
        case UPWIND:
          aux = new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
          break;
        default:
          // parameters needed for assembling are (u1_old, u2_old, u3_old)
          aux =  new TAuxParam3D(TimeNSN_FESpacesVelo,
            TimeNSN_FctVelo,
            TimeNSN_ParamFctVelo,
            TimeNSN_FEValuesVelo,
            fesp, fefct,
            TimeNSFctVelo,
            TimeNSFEFctIndexVelo, TimeNSFEMultiIndexVelo,
            TimeNSN_ParamsVelo, TimeNSBeginParamVelo);
      }
    }

    //======================================================================
    // assembling of matrices for each level
    // A_11 , (A_12), (A_21), (A_22), M_11, (M_22)
    // assembling of rhs not needed at this point
    //======================================================================
    Assemble3D(N_FESpaces, fesp,
      N_SquareMatrices, SQMATRICES,
      N_RectMatrices, MATRICES,
      N_Rhs, RHSs, ferhs,
      DiscreteForm,
      BoundaryConditions,
      BoundValues,
      aux);
    t2 = GetTime();
    ass_time += t2-t1;
    OutPut("time for assembling " << ass_time  << "s "<<  endl);
    delete aux;

    // copy Dirichlet values from rhs into sol
    memcpy(sol+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
    memcpy(sol+N_Active+N_U, RHSs[1]+N_Active, (N_U-N_Active)*SizeOfDouble);
    memcpy(sol+N_Active+2*N_U, RHSs[2]+N_Active, (N_U-N_Active)*SizeOfDouble);

    // add convective term in the upwind discretizations
    if ((DiscreteForm == DiscreteFormUpwind)||(DiscreteForm == DiscreteFormUpwindNC))
    {
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
        case 2:
          // do upwinding with one matrix
          UpwindForNavierStokes3D(SQMATRICES[0], U1Array[i], U2Array[i], U3Array[i]);
          //cout << "UPWINDING DONE : level " << i << endl;
          break;

        case 3:
        case 4:
        case 14:
          // do upwinding with two matrices
          UpwindForNavierStokes3D(SQMATRICES[0], U1Array[i], U2Array[i], U3Array[i]);
          UpwindForNavierStokes3D(SQMATRICES[4], U1Array[i], U2Array[i], U3Array[i]);
          UpwindForNavierStokes3D(SQMATRICES[8], U1Array[i], U2Array[i], U3Array[i]);
          cout << "UPWINDING DONE(2) : level " << i << endl;
          break;
      }                                           // endswitch
    }

    // set rows of Dirichlet dof in off diagonal matrix blocks
    // to zero
    switch(TDatabase::ParamDB->NSTYPE)
    {
      case 3:
      case 4:
      case 14:
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
        memset(MatricesA13[i]->GetEntries()+j, 0, SizeOfDouble*k);
        memset(MatricesA21[i]->GetEntries()+j, 0, SizeOfDouble*k);
        memset(MatricesA23[i]->GetEntries()+j, 0, SizeOfDouble*k);
        memset(MatricesA31[i]->GetEntries()+j, 0, SizeOfDouble*k);
        memset(MatricesA32[i]->GetEntries()+j, 0, SizeOfDouble*k);
        break;
    }
    if ((CurrentDiscType == VMS_PROJECTION)||
      (CurrentDiscType == VMS_PROJECTION_EXPL))
    {
      LumpMassMatrixToDiag(MatricesL[i]);
      SQMATRICES[0] = MatricesA11[i];
      SQMATRICES[1] = MatricesA12[i];
      SQMATRICES[2] = MatricesA13[i];
      SQMATRICES[3] = MatricesA21[i];
      SQMATRICES[4] = MatricesA22[i];
      SQMATRICES[5] = MatricesA23[i];
      SQMATRICES[6] = MatricesA31[i];
      SQMATRICES[7] = MatricesA32[i];
      SQMATRICES[8] = MatricesA33[i];
      SQMATRICES[9] =  MatricesL[i];
      MATRICES[0] = Matrices_tilde_G11[i];
      MATRICES[1] = Matrices_tilde_G22[i];
      MATRICES[2] = Matrices_tilde_G33[i];
      MATRICES[3] = Matrices_G11[i];
      MATRICES[4] = Matrices_G22[i];
      MATRICES[5] = Matrices_G33[i];
      if (CurrentDiscType == VMS_PROJECTION)
      {
        if ((i==mg_level - 1)||(!TDatabase::ParamDB->VMS_COARSE_MG_SMAGO))
        {
          /*OutPut( Ddot(Matrices_tilde_G11[i]->GetN_Entries(),
             Matrices_tilde_G11[i]->GetEntries(),Matrices_tilde_G11[i]->GetEntries()) << " ");
          OutPut( Ddot(Matrices_tilde_G22[i]->GetN_Entries(),
             Matrices_tilde_G22[i]->GetEntries(),Matrices_tilde_G22[i]->GetEntries()) << " ");
          OutPut( Ddot(Matrices_tilde_G33[i]->GetN_Entries(),
             Matrices_tilde_G33[i]->GetEntries(),Matrices_tilde_G33[i]->GetEntries()) << " ");
          OutPut( Ddot(Matrices_G11[i]->GetN_Entries(),
             Matrices_G11[i]->GetEntries(),Matrices_G11[i]->GetEntries()) << " ");
          OutPut( Ddot(Matrices_G22[i]->GetN_Entries(),
             Matrices_G22[i]->GetEntries(),Matrices_G22[i]->GetEntries()) << " ");
          OutPut( Ddot(Matrices_G33[i]->GetN_Entries(),
          Matrices_G33[i]->GetEntries(),Matrices_G33[i]->GetEntries()) << endl);*/
          VMS_ProjectionUpdateMatrices(N_U,N_Active,N_L,SQMATRICES,MATRICES);
          OutPut("update done"<<endl);
        }
      }
      else
      {
        if (i==mg_level-1)
        {
          rhs_vms_expl = new double[3*N_U*SizeOfDouble];
          OutPut("update rhs"<<endl);
          VMS_ProjectionExplUpdateRhs(N_U,N_Active, N_L, u,
            SQMATRICES, MATRICES, rhs_vms_expl);
          OutPut("update rhs done " << Ddot(3*N_U,rhs_vms_expl,rhs_vms_expl) << endl);
        }
      }
    }
  }                                               // endfor i  (line 705)
  t4 =  GetTime();
  total_time += t4 - t3;
  t3 = t4;

  if (TDatabase::ParamDB->DISCTYPE == VMS_RFB_EXPL_COUPLED)
  {
    N_sub = TDatabase::ParamDB->RFB_SUBMESH_LAYERS;
    N_sub = 3 * (2*N_sub-1)*(2*N_sub-1)*(2*N_sub-1) + (N_sub+1)*(N_sub+1)*(N_sub+1);
    N_sub *= coll->GetN_Cells();
    OutPut("dof for coupled rfb: " <<  N_sub << endl);
    old_small_scales = new double[N_sub];
    memset(old_small_scales, 0, N_sub*SizeOfDouble);
  }

  if  (TDatabase::ParamDB->DISCTYPE == VMS_RFB_EXPL)
  {
    N_sub = TDatabase::ParamDB->RFB_SUBMESH_LAYERS;
    N_sub = 3 * (N_sub-1)*(N_sub-1)*(N_sub-1);
    N_sub *= coll->GetN_Cells();
    OutPut("dof for rfb: " <<  N_sub << endl);
    old_small_scales = new double[N_sub];
    memset(old_small_scales, 0, N_sub*SizeOfDouble);
  }

  //======================================================================
  // end of space cycle, finest grid reached
  // everything happens on the same grid
  //======================================================================
  // copy sol for extrapolation after time step
  if (TDatabase::TimeDB->EXTRAPOLATE_VELOCITY)
    memcpy(sol_timestep_m1,sol,N_Unknowns*SizeOfDouble);

  comp_vort =0;
  if (TDatabase::ParamDB->COMPUTE_VORTICITY_DIVERGENCE)
  {
    if (!comp_vort)
      ComputeVorticityDivergence(USpaces[mg_level-1], U1Array[mg_level-1],
        U2Array[mg_level-1], U3Array[mg_level-1],
        vorticity_space, Vorticity->GetComponent(0),
        Vorticity->GetComponent(1),
        Vorticity->GetComponent(2),
        div);
    comp_vort++;
  }

#ifdef  __MIXINGLAYERSLIP3D__
  ComputeVorticityThickness3D(Vorticity->GetComponent(2),errors);
  vort_zero = errors[0];
  OutPut( TDatabase::TimeDB->CURRENTTIME << " vorticity thickness: " << errors[0] << " " <<
    errors[0]/vort_zero << endl);

  ComputeMomentumThickness3D(U1Array[mg_level-1],errors);
  moment_zero = errors[0];
  OutPut( TDatabase::TimeDB->CURRENTTIME << " momentum thickness: " << errors[0] << " " <<
    errors[0]/moment_zero << endl);
  comp_vort =0;

  // fesp, fefct not needed here
  aux =  new TAuxParam3D(TimeNSN_FESpacesVelo, TimeNSN_FctVelo,
    TimeNSN_ParamFctVelo,
    TimeNSN_FEValuesVelo,
    fesp, fefct,
    TimeNSFctVelo,
    TimeNSFEFctIndexVelo, TimeNSFEMultiIndexVelo,
    TimeNSN_ParamsVelo, TimeNSBeginParamVelo);
  // enstrophy
  Vorticity->GetComponent(0)->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
    2, L2H1Errors, NULL, aux, 1, &vorticity_space, errors);
  Vorticity->GetComponent(1)->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
    2, L2H1Errors, NULL, aux, 1, &vorticity_space, errors+2);

  Vorticity->GetComponent(2)->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
    2, L2H1Errors, NULL, aux, 1, &vorticity_space, errors+4);
  OutPut( TDatabase::TimeDB->CURRENTTIME << " enstrophy: " << (errors[0]*errors[0]+errors[2]*errors[2]
    +errors[4]*errors[4])/2.0 << endl);
  delete aux;
#endif

  if (TDatabase::ParamDB->CONVOLUTE_SOLUTION)
  {
    ConvolveSolution(MG, USpaces,
      U1Array, U2Array, U3Array,
      DiscreteFormRHSAuxProblemU,
      sqmatrixGL00AuxProblem,
      rhsGL00AuxProblem,
      u_uConv,
      mg_level, N_U);

    aux =  new TAuxParam3D(TimeNSN_FESpacesVelo, TimeNSN_FctVelo,
      TimeNSN_ParamFctVelo,
      TimeNSN_FEValuesVelo,
      fesp, fefct,
      TimeNSFctVelo,
      TimeNSFEFctIndexVelo, TimeNSFEMultiIndexVelo,
      TimeNSN_ParamsVelo, TimeNSBeginParamVelo);

    // error in first component
    u1ConvArray[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
      2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors);

    // error in second component
    u2ConvArray[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
      2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors+2);

    // error in third component
    u3ConvArray[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
      2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors+4);

    OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
    OutPut( "conv H1-semi(u): " << sqrt(errors[1]*errors[1]+errors[3]*errors[3]
      +errors[5]*errors[5]) << endl);
    OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
    OutPut( "conv kinetic energy " << (errors[0]*errors[0]+errors[2]*errors[2]
      +errors[4]*errors[4])/2<< endl);

    if (TDatabase::ParamDB->COMPUTE_VORTICITY_DIVERGENCE)
    {                                             // compute vorticity
      ComputeVorticityDivergence(velocity_space, u1ConvArray[mg_level-1],
        u2ConvArray[mg_level-1], u3ConvArray[mg_level-1],
        vorticity_space,Vorticity->GetComponent(3),
        Vorticity->GetComponent(4),
        Vorticity->GetComponent(5),div);
    }

#ifdef __MIXINGLAYERSLIP3D__
    ComputeVorticityThickness3D(Vorticity->GetComponent(5),errors);
    OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
    vort_zero_conv =  errors[0];
    OutPut( "vorticity thickness (uconv): " << errors[0] << " " <<
      errors[0]/vort_zero_conv << endl);
    ComputeMomentumThickness3D(u1ConvArray[mg_level-1],errors);
    OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
    moment_zero_conv =  errors[0];
    OutPut( "momentum thickness (uconv): " << errors[0] << " " <<
      errors[0]/moment_zero_conv << endl);

    // enstrophy
    Vorticity->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
      2, L2H1Errors,
      NULL, aux, 1,  &vorticity_space, errors);

    OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
    OutPut( "enstrophy (uconv) " << (errors[0]*errors[0]+errors[2]*errors[2]
      +errors[4]*errors[4])/2.0  << endl);

    // smooth filtered vorticity for output
    // copy vorticity to sol to use mg structure
    if (mixing_layer_galerkin==1)
    {
      OutPut("smooth solution " << UArray[mg_level-1]->GetValues()[0] << endl);
      memcpy(sol_vort_tmp,UArray[mg_level-1]->GetValues(),3*N_U*SizeOfDouble);
      memcpy(UArray[mg_level-1]->GetValues(),Vorticity->GetComponent(3)->GetValues(),3*N_U*SizeOfDouble);
      // restrict
      auxConv = new double[3*N_U];
      for(ii = mg_level-1 ; ii > mg_level-1-smoothing_depth;ii--)
      {
        RestrictFunction(USpaces[ii-1], USpaces[ii], 3,
          UArray[ii-1]->GetValues(),
          UArray[ii]->GetValues(),
          auxConv);
      }
      // prolongate restricted vorticity
      for(ii = mg_level-1-smoothing_depth ; ii < mg_level-1;ii++)
      {
        Prolongate(USpaces[ii], USpaces[ii+1], 3,
          UArray[ii]->GetValues(),
          UArray[ii+1]->GetValues(),
          auxConv);

      }
      // copy result to vorticity vector
      memcpy(Vorticity->GetComponent(3)->GetValues(),UArray[mg_level-1]->GetValues(),3*N_U*SizeOfDouble);
      memcpy(UArray[mg_level-1]->GetValues(),sol_vort_tmp,3*N_U*SizeOfDouble);
      // restrict solution to all grids
      MG->RestrictToAllGrids();
      delete auxConv;
    }
#endif

    delete aux;
  }

  if ((TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION) ||
    (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION_EXPL))
  {
    /*for (i=0;i<N_U;i++)
    {
        U1Array[mg_level-1]->GetValues()[i] = x_dof[i]-4*y_dof[i];
        U2Array[mg_level-1]->GetValues()[i] = -x_dof[i]+2*y_dof[i];
        U3Array[mg_level-1]->GetValues()[i] = 0;
        }*/
    ComputeVMSProjection(Matrices_G11[mg_level-1], Matrices_G22[mg_level-1],
      Matrices_G33[mg_level-1], MatricesL[mg_level-1],
      U1Array[mg_level-1], U2Array[mg_level-1],
      U3Array[mg_level-1], vms_projection_fe);
    ComputeSizeOfSmallScales(Matrices_G11[mg_level-1], Matrices_G22[mg_level-1],
      Matrices_G33[mg_level-1], MatricesL[mg_level-1],
      U1Array[mg_level-1], U2Array[mg_level-1],
      U3Array[mg_level-1], vms_projection_fe, size_small_scales);
    MeanAndLargestSize( projection_space, size_small_scales, &mean, &largest_size);
    mean_time_average = mean;
    max_time_average = largest_size;

    if (TDatabase::ParamDB->VMS_LARGE_VELOCITY_SPACE==17)
    {
      AdaptProjectionSpace(projection_space, size_small_scales, fes,
        mean, mean_time_average, largest_size,
        max_time_average, label_space);
      /*for (i=0;i<N_L;i++)
      {
          OutPut(vms_projection_fe->GetValues()[i+1*N_L]<<" ");
      }
      exit(1);*/
    }
  }
#ifdef __CHANNEL_TAU180__
  TDatabase::TimeDB->CURRENTTIME = 0.0;
  OutPut("reset current time to 0"<<endl);

  TDatabase::ParamDB->INTERNAL_MEAN_COMPUTATION = 1;

  ComputeMeanVelocity(coll, USpaces[mg_level-1], projection_space,
    N_z_layers, coord_z_layers,  rms_velocity,
    rms_velocity2, rms_velocity3, rms_velocity1_type1,
    rms_velocity2_type1, rms_velocity3_type1,
    mean_velocity, mean_velocity_u2, mean_velocity_u3,
    dmean_velocity,
    R_xx, R_xy, R_xz, R_yy, R_yz, R_zz,
    A_M_xx, A_M_xy, A_M_xz, A_M_yy, A_M_yz, A_M_zz,
    U1Array[mg_level-1], U2Array[mg_level-1], U3Array[mg_level-1],
    vms_projection_fe,
    x_dof, y_dof, z_dof, u1x, projection_u1x);
#endif

  if ((TDatabase::ParamDB->WRITE_GRAPE)||(TDatabase::ParamDB->WRITE_GMV)
    ||(TDatabase::ParamDB->WRITE_AMIRA)||(TDatabase::ParamDB->WRITE_VTK))
  {
    if (TDatabase::ParamDB->COMPUTE_VORTICITY_DIVERGENCE)
    {
      if (!comp_vort)
        ComputeVorticityDivergence(USpaces[mg_level-1], U1Array[mg_level-1],
          U2Array[mg_level-1], U3Array[mg_level-1],
          vorticity_space, Vorticity->GetComponent(0),
          Vorticity->GetComponent(1),
          Vorticity->GetComponent(2),
          div);
      comp_vort++;
    }

    if (TDatabase::ParamDB->WRITE_GMV)
    {
      os.seekp(std::ios::beg);
      os << GMVBaseName << 0 << ".gmv" << ends;
      Output->WriteGMV(os.str().c_str());
    }
    if (TDatabase::ParamDB->WRITE_VTK)
    {
      os.seekp(std::ios::beg);
      os << VTKBaseName << 0 << ".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
    }
    if (TDatabase::ParamDB->WRITE_AMIRA)
    {
      os.seekp(std::ios::beg);
      os << GrapeBaseName << 0 << ".am" << ends;
      Output->WriteAmira(os.str().c_str());
    }
    if (TDatabase::ParamDB->WRITE_GRAPE)
    {
      os.seekp(std::ios::beg);
      os << GrapeBaseName << 0 << ".dat" << ends;
      Output->WriteGrape(os.str().c_str());
    }
    N_GRAPE_images++;
  }

  defect = new double[N_Unknowns];
  startsol = new double[N_Unknowns];
  frac_step_sol = new double[N_Unknowns];
  oldrhs =  new double[N_Unknowns];

  // Newton's method
  // the entry on the rhs which is due to Newton's method is stored in the
  // array newton
  // this array is assembled separately, added before solving the system to
  // rhs and subtracted afterwards
  if (nonlinite==1)
    newton = new double[3*N_Uarray[mg_level-1]];

  N_Active = velocity_space->GetActiveBound();

  solver_time = 0.0;
  N_LinIter = 0;

  gamma = 0;
  m = 0;
  N_SubSteps = GetN_SubSteps();
  oldtau = 1.0;

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
  //  #ifdef __CALOTTE__
  //  CheckNeumannNodesForVelocity(coll, USpaces[mg_level-1], N_neum_to_diri, neum_to_diri,
  //			       N_neum_bdry, neum_bdry);
  //  #endif

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
      OutPut("delta is non-constant" << endl);
    }
  }

  // check difference of initial condition and prescribed exact solution
  // necessary for updating l_infty_l_2_time, l_2_l_2u, l_2_h_1u if
  // some quantities of the computed solution are measured
  if (TDatabase::ParamDB->MEASURE_ERRORS)
  {
    fesp[0] = USpaces[mg_level-1];
    fefct[0] = U1Array[mg_level-1];
    fefct[1] = U2Array[mg_level-1];
    fefct[2] = U3Array[mg_level-1];

    aux =  new TAuxParam3D(TimeNSN_FESpacesVelo, TimeNSN_FctVelo,
      TimeNSN_ParamFctVelo,
      TimeNSN_FEValuesVelo,
      fesp, fefct,
      TimeNSFctVelo,
      TimeNSFEFctIndexVelo, TimeNSFEMultiIndexVelo,
      TimeNSN_ParamsVelo, TimeNSBeginParamVelo);

    // errors
    // L2: error[0], H1-semi: error[1]
    U1Array[mg_level-1]->GetErrors(ExactU1, 4, TimeNSAllDerivatives,
      2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors);

    // L2: error[2], H1-semi: error[3]
    U2Array[mg_level-1]->GetErrors(ExactU2, 4, TimeNSAllDerivatives,
      2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors+2);

    // L2: error[4], H1-semi: error[5]
    U3Array[mg_level-1]->GetErrors(ExactU3, 4, TimeNSAllDerivatives,
      2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors+4);

    // error in L^infty(0,t,L^2)
    if (sqrt(errors[0]*errors[0]+errors[2]*errors[2]+errors[4]*errors[4])
      > l_infty_l_2)
    {
      l_infty_l_2 = sqrt(errors[0]*errors[0]+errors[2]*errors[2]
        +errors[4]*errors[4]);
      l_infty_l_2_time =  TDatabase::TimeDB->CURRENTTIME;
    }

    // error in L^2(0,t,L^2)
    olderror_l_2_l_2u = errors[0]*errors[0] + errors[2]*errors[2]+errors[4]*errors[4];

    //error in L^2(0,t,H^1)
    olderror_l_2_h_1u = errors[1]*errors[1] + errors[3]*errors[3]+ errors[5]*errors[5];
  }

  //======================================================================
  // start of time cycle
  // everything happens on the same grid
  //======================================================================

  while(TDatabase::TimeDB->CURRENTTIME< end_time)
  {                                               // time cycle

    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

    for (methods=0;methods<time_discs;methods++)
    {
      if (time_discs==2)
      {
        if (methods==0)                           // fractional-step-theta-scheme
        {
          TDatabase::TimeDB->TIME_DISC = 3;
          // save start sol (nec. for gl00)
          memcpy(startsol,sol,N_Unknowns*SizeOfDouble);
          // save rhs
          memcpy(oldrhs,rhs,N_Unknowns*SizeOfDouble);
        }
        else                                      // crank nicolson scheme
        {                                         // take solution of first scheme as initial iterate
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

      for(l=0;l<N_SubSteps;l++)                   // sub steps of fractional step theta
      {
        t11 = GetTime();
        if (!very_first_time)
        {
          SetTimeDiscParameters();
          theta1 = TDatabase::TimeDB->THETA1;
          theta2 = TDatabase::TimeDB->THETA2;
          theta3 = TDatabase::TimeDB->THETA3;
          theta4 = TDatabase::TimeDB->THETA4;
        }
        if (m==1)
        {
          OutPut("Theta1: " << theta1<< endl);
          OutPut("Theta2: " << theta2<< endl);
          OutPut("Theta3: " << theta3<< endl);
          OutPut("Theta4: " << theta4<< endl);
        }
        tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
        if (!very_first_time)
          TDatabase::TimeDB->CURRENTTIME += tau;
        // working array for rhs is B, initialize B
        memset(B, 0, N_Unknowns*SizeOfDouble);

        // old rhs array (f) multiplied with current subtime step and theta3 on B
        Daxpy(N_Active, tau*theta3, rhs, B);
        Daxpy(N_Active, tau*theta3, rhs+N_U, B+N_U);
        Daxpy(N_Active, tau*theta3, rhs+2*N_U, B+2*N_U);

        //======================================================================
        // prepare input data (parameters) for several discretizations

        // compute convolution of \nabla u \nabla u^T
        if (TDatabase::ParamDB->DISCTYPE==GL00_CONVOLUTION)
          ComputeConvolutionOfNabla_uNabla_uTrans3D(MG, UArray, duConv,
            duConvSpaces,
            du11ConvArray,
            du12ConvArray,
            du13ConvArray,
            du22ConvArray,
            du23ConvArray,
            du33ConvArray,
            mg_level, N_Unknowns);

        // compute convolution of u for |u-g_\delta\ast(u)|
        // with auxiliary problem
        if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE==4)
        {
          ComputeConvolutionForTurbVisType4(MG,
            USpaces,
            U1Array, U2Array, U3Array,
            uConvArray,
            u1ConvArray, u2ConvArray, u3ConvArray,
            DiscreteFormRHSAuxProblemU,
            sqmatrixGL00AuxProblem,
            mg_level, N_U);
        }

        //======================================================================
        // assembling of the rhs of current sub time step
        // only on the finest level necessary
        // there is no assembling of matrices here
        //======================================================================

        fesp[0] = USpaces[mg_level-1];
        RHSs[0] = RhsArray[mg_level-1];
        RHSs[1] = RhsArray[mg_level-1]+N_Uarray[mg_level-1];
        RHSs[2] = RhsArray[mg_level-1]+2*N_Uarray[mg_level-1];

        ferhs[0] = USpaces[mg_level-1];
        ferhs[1] = USpaces[mg_level-1];
        ferhs[2] = USpaces[mg_level-1];

        if (nonlinite==0)
          N_Rhs = 3;
        else                                      // Newton's method
        {
          RHSs[3] = newton;
          RHSs[4] = newton+N_Uarray[mg_level-1];
          RHSs[5] = newton+2*N_Uarray[mg_level-1];

          ferhs[3] = USpaces[mg_level-1];
          ferhs[4] = USpaces[mg_level-1];
          ferhs[5] = USpaces[mg_level-1];

          N_Rhs = 6;
        }

        N_FESpaces = 1;
        N_SquareMatrices = 0;
        N_RectMatrices = 0;

        fefct[0] = U1Array[mg_level-1];
        fefct[1] = U2Array[mg_level-1];
        fefct[2] = U3Array[mg_level-1];

        switch(TDatabase::ParamDB->DISCTYPE)
        {
          case CLASSICAL_LES :
          case GL00_CONVOLUTION :
          case GL00_AUX_PROBLEM :
            PrepareRHSLES(USpaces,
              U1Array, U2Array, U3Array,
              uConvSpaces,
              u1ConvArray, u2ConvArray, u3ConvArray,
              duConvSpaces,
              du11ConvArray, du12ConvArray, du13ConvArray,
              du22ConvArray, du23ConvArray, du33ConvArray,
              GL00AuxProblemSol11Array, GL00AuxProblemSol12Array,
              GL00AuxProblemSol13Array, GL00AuxProblemSol22Array,
              GL00AuxProblemSol23Array, GL00AuxProblemSol33Array,
              DiscreteFormRHSClassicalLES,
              DiscreteFormRHSLESModel,
              DiscreteFormGL00AuxProblemRHS,
              BoundaryConditions, BoundValues,
              sqmatrixGL00AuxProblem,
              rhs, solGL00AuxProblem, LESModelRhs,
              mg_level, N_U, N_P);
            break;
          case VMS_PROJECTION_EXPL:
            // assemble Matrices_tilde_G??
            N_FESpaces = 3;
            fesp[0] = USpaces[mg_level-1];
            fesp[1] = ProjectionSpaces[mg_level-1];
            fesp[2] = label_space_fesp;

            fefct[0] = U1Array[mg_level-1];
            fefct[1] = U2Array[mg_level-1];
            fefct[2] = U3Array[mg_level-1];

            TDatabase::ParamDB->INTERNAL_LEVEL = 1;
            fefct[3] = vms_proj_11;
            fefct[4] = vms_proj_12;
            fefct[5] = vms_proj_13;
            fefct[6] = vms_proj_22;
            fefct[7] = vms_proj_23;
            fefct[8] = vms_proj_33;
            fefct[9] = label_space_fefct;
            // compute large scales
            /*for (i=0;i<N_U;i++)
              {
              U1Array[mg_level-1]->GetValues()[i] = x_dof[i]-4*y_dof[i];
              U2Array[mg_level-1]->GetValues()[i] = -x_dof[i]+2*y_dof[i];
              U3Array[mg_level-1]->GetValues()[i] = 0;
              }*/
            ComputeVMSProjection(Matrices_G11[mg_level-1], Matrices_G22[mg_level-1],
              Matrices_G33[mg_level-1], MatricesL[mg_level-1],
              U1Array[mg_level-1], U2Array[mg_level-1],
              U3Array[mg_level-1], vms_projection_fe);

            /*for (i=0;i<N_L;i++)
              {
              OutPut(vms_projection_fe->GetValues()[i+0*N_L]<<" ");
              }*/
            // turbulent viscosity for the large scales
            aux =  new TAuxParam3D(TimeNSN_FESpacesVelo_GradVelo_LargeScale,
              TimeNSN_FctVelo_GradVelo_LargeScale,
              TimeNSN_ParamFctVelo_GradVelo_LargeScale,
              TimeNSN_FEValuesVelo_GradVelo_LargeScale,
              fesp, fefct,
              TimeNSFctVelo_GradVelo_LargeScale,
              TimeNSFEFctIndexVelo_GradVelo_LargeScale,
              TimeNSFEMultiIndexVelo_GradVelo_LargeScale,
              TimeNSN_ParamsVelo_GradVelo_LargeScale,
              TimeNSBeginParamVelo_GradVelo_LargeScale);

            /*aux =  new TAuxParam3D(TimeNSN_FESpacesVelo_GradVelo,
              TimeNSN_FctVelo_GradVelo,
              TimeNSN_ParamFctVelo_GradVelo,
              TimeNSN_FEValuesVelo_GradVelo,
              fesp, fefct,
              TimeNSFctVelo_GradVelo,
              TimeNSFEFctIndexVelo_GradVelo,
              TimeNSFEMultiIndexVelo_GradVelo,
              TimeNSN_ParamsVelo_GradVelo,
              TimeNSBeginParamVelo_GradVelo);*/

            MATRICES[0] = Matrices_tilde_G11[mg_level-1];
            MATRICES[1] = Matrices_tilde_G22[mg_level-1];
            MATRICES[2] = Matrices_tilde_G33[mg_level-1];
            MATRICES[0]->Reset();
            MATRICES[1]->Reset();
            MATRICES[2]->Reset();

            if (turb_visc_expl_change == 1)
            {
              TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE=5;
            }

            Assemble3D(N_FESpaces, fesp,
              0, NULL,
              3, MATRICES,
              0, NULL, ferhs,
              DiscreteFormNLVMS_ProjectionExpl,
              BoundaryConditions,
              BoundValues,
              aux);
            delete aux;

            // small scales
            if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE==5)
            {
              turb_visc_expl_change = 1;
              TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE=1;
            }
            LumpMassMatrixToDiag(MatricesL[mg_level-1]);
            SQMATRICES[9] =  MatricesL[mg_level-1];
            MATRICES[0] = Matrices_tilde_G11[mg_level-1];
            MATRICES[1] = Matrices_tilde_G22[mg_level-1];
            MATRICES[2] = Matrices_tilde_G33[mg_level-1];
            MATRICES[3] = Matrices_G11[mg_level-1];
            MATRICES[4] = Matrices_G22[mg_level-1];
            MATRICES[5] = Matrices_G33[mg_level-1];
            VMS_ProjectionExplUpdateRhs(N_U,N_Active, N_L, UArray[mg_level-1],
              SQMATRICES, MATRICES, rhs_vms_expl);
            OutPut("update rhs done " << Ddot(3*N_U,rhs_vms_expl,rhs_vms_expl) << endl);
            break;
        }

        // assemble rhs from f
        DiscreteForm = DiscreteFormRHS;
        aux = new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

        // initialize array
        memset(RHSs[0], 0,
          (3*N_Uarray[mg_level-1]+N_Parray[mg_level-1])*SizeOfDouble);
        // all input data are computed, assemble now rhs

        // only for fixed point iteration
        if (TDatabase::ParamDB->DISCTYPE == SDFEM)
        {
          fesp[1] = PSpaces[mg_level-1];
          RHSs[3] = RhsArray[mg_level-1]+3*N_Uarray[mg_level-1];
          ferhs[3] = PSpaces[mg_level-1];
          fefct[3] = PArray[mg_level-1];
          N_Rhs = 4;
          N_FESpaces = 2;
          aux =  new TAuxParam3D(TimeNSN_FESpacesVelo_GradVelo,
            TimeNSN_FctVelo_GradVelo,
            TimeNSN_ParamFctVelo_GradVelo,
            TimeNSN_FEValuesVelo_GradVelo,
            fesp, fefct,
            TimeNSFctVelo_GradVelo,
            TimeNSFEFctIndexVelo_GradVelo,
            TimeNSFEMultiIndexVelo_GradVelo,
            TimeNSN_ParamsVelo_GradVelo,
            TimeNSBeginParamVelo_GradVelo);
          DiscreteForm = DiscreteFormRHSSUPG;
        }

        if (nonlinite==1)                         // Newton's method
        {
          DiscreteForm = DiscreteFormRHSNewton;
          aux =  new TAuxParam3D(TimeNSN_FESpacesVelo_GradVelo,
            TimeNSN_FctVelo_GradVelo,
            TimeNSN_ParamFctVelo_GradVelo,
            TimeNSN_FEValuesVelo_GradVelo,
            fesp, fefct,
            TimeNSFctVelo_GradVelo,
            TimeNSFEFctIndexVelo_GradVelo,
            TimeNSFEMultiIndexVelo_GradVelo,
            TimeNSN_ParamsVelo_GradVelo,
            TimeNSBeginParamVelo_GradVelo);

          memset(newton, 0, 3*N_Uarray[mg_level-1]*SizeOfDouble);
        }

        Assemble3D(N_FESpaces, fesp,
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
        Daxpy(N_Active, tau*theta4, rhs, B);
        Daxpy(N_Active, tau*theta4, rhs+N_U, B+N_U);
        Daxpy(N_Active, tau*theta4, rhs+2*N_U, B+2*N_U);

        // update rhs in LES models
        if ((TDatabase::ParamDB->DISCTYPE==CLASSICAL_LES)||
          (TDatabase::ParamDB->DISCTYPE==GL00_AUX_PROBLEM)||
          (TDatabase::ParamDB->DISCTYPE==GL00_CONVOLUTION))
          Daxpy(3*N_U, tau, LESModelRhs, B);

        // update rhs in explicit projection-based VMS
        if (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION_EXPL)
          Daxpy(3*N_U, tau, rhs_vms_expl, B);

        // update rhs in explicit bubble VMS
        if (TDatabase::ParamDB->DISCTYPE == VMS_RFB_EXPL)
        {
          t111 = GetTime();
          if (TDatabase::ParamDB->P5!=10)
          {

            ApproximateTimeRFBSolutionQuadNSE3D(coll, U1Array[mg_level-1], U2Array[mg_level-1],
              U3Array[mg_level-1], PArray[mg_level-1],
              Coefficients[0], B);
            t222 = GetTime();
            OutPut(" RFB DONE, time " << t222-t111 << "s"<<endl);
          }
          else
          {
            ApproximateTimeRFBSolutionQuad_cn_NSE3D(coll, U1Array[mg_level-1], U2Array[mg_level-1],
              U3Array[mg_level-1], PArray[mg_level-1],
              Coefficients[0], old_small_scales, B);
            t22 = GetTime();
            OutPut(" RFB CN DONE, time  "<< t222-t111 << "s"<< endl);
          }
        }
        if (TDatabase::ParamDB->DISCTYPE == VMS_RFB_EXPL_COUPLED)
        {
          if (TDatabase::ParamDB->P5!=10)
          {
            ApproximateTimeRFB_coupled_SolutionQuad_Q2_NSE3D(coll, U1Array[mg_level-1], U2Array[mg_level-1],
              U3Array[mg_level-1], PArray[mg_level-1],
              Coefficients[0], B);
            OutPut(" RFB_COUPLED_RED DONE"<<endl);
          }
          else
          {
            ApproximateTimeRFB_coupled_cn_SolutionQuad_Q2_NSE3D(coll, U1Array[mg_level-1], U2Array[mg_level-1],
              U3Array[mg_level-1], PArray[mg_level-1],
              Coefficients[0], old_small_scales, B);
            OutPut(" RFB_COUPLED_CN DONE"<<endl);
          }
        }

        // Newton's method, add entries which come from Newton's method to rhs
        if (nonlinite==1)
        {
          Daxpy(N_Active, tau*(theta1+theta2), newton, B);
          Daxpy(N_Active, tau*(theta1+theta2), newton+N_U, B+N_U);
          Daxpy(N_Active, tau*(theta1+theta2), newton+2*N_U, B+2*N_U);
        }

        // do not change the array rhs during nonlinear iteration !!!
        // needed for next time step !!!

        // slip type bc detected, manipulation of matrices is necessary
        // this is done only at the very beginning
        // the matrices A_12, A_21, ..., M_11, M_12, M_21, M_22, ...,
        //     B1T, B2T, B3T
        //     stay unchanged during the complete solution process
        // the matrices A_11, A_22, A33 are manipulated after their new
        //     assembling during the nonlinear iteration

        if ((m==1)&& (l==0) &&
          (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION == 1))
        {
          if (TDatabase::ParamDB->NSTYPE <4)
          {
            OutPut("For slip with friction bc NSTYPE 4 or 14 is ");
            OutPut("necessary !!!!! " << endl);
            exit(4711);
          }
          // prepare everything for the assembling of slip with friction bc
          // on all levels

          N_FESpaces = 1;
          N_SquareMatrices = 18;
          N_RectMatrices = 3;
          N_Rhs = 3;
          DiscreteForm = NULL;

          for(i=0;i<mg_level;i++)
          {
            SQMATRICES[0] = MatricesA11[i];
            SQMATRICES[1] = MatricesA22[i];
            SQMATRICES[2] = MatricesA33[i];
            SQMATRICES[3] = MatricesA12[i];
            SQMATRICES[4] = MatricesA13[i];
            SQMATRICES[5] = MatricesA21[i];
            SQMATRICES[6] = MatricesA23[i];
            SQMATRICES[7] = MatricesA31[i];
            SQMATRICES[8] = MatricesA32[i];

            SQMATRICES[9] = MatricesM11[i];
            SQMATRICES[10] = MatricesM22[i];
            SQMATRICES[11] = MatricesM33[i];
            SQMATRICES[12] = MatricesM12[i];
            SQMATRICES[13] = MatricesM13[i];
            SQMATRICES[14] = MatricesM21[i];
            SQMATRICES[15] = MatricesM23[i];
            SQMATRICES[16] = MatricesM31[i];
            SQMATRICES[17] = MatricesM32[i];

            MATRICES[0] = MatricesB1T[i];
            MATRICES[1] = MatricesB2T[i];
            MATRICES[2] = MatricesB3T[i];

            fesp[0] = USpaces[i];
            ferhs[0] = USpaces[i];
            ferhs[1] = USpaces[i];
            ferhs[2] = USpaces[i];

            RHSs[0] = RhsArray[i];
            RHSs[1] = RhsArray[i]+N_Uarray[i];
            RHSs[2] = RhsArray[i]+2*N_Uarray[i];

            Assemble3DSlipBC(N_FESpaces, fesp,
              N_SquareMatrices, SQMATRICES,
              N_RectMatrices, MATRICES,
              N_Rhs, RHSs, ferhs,
              DiscreteForm,
              BoundaryConditions,
              BoundValues,
              aux);
            TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;

#ifdef __BENCH__
            SetNoPenetrationValues(SQMATRICES, UArray[i], hmin);
            // copy sol for extrapolation after time step
            if (TDatabase::TimeDB->EXTRAPOLATE_VELOCITY)
              memcpy(sol_timestep_m1,sol,N_Unknowns*SizeOfDouble);
#endif
          }
        }
        delete aux;

        //======================================================================
        // manipulation of matrices due to current time discretization
        // the stiffness matrix is stored on M11, (M12, M21, M22)
        //======================================================================

        // scale B1T, B2T, B3T, B1, B2, B3
        for(i=0;i<mg_level;i++)
        {
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
            case 3:
              if (tau/oldtau != 1.0)
              {
                Dscal(MatricesB1[i]->GetN_Entries(),
                  tau/oldtau,
                  MatricesB1[i]->GetEntries());
                Dscal(MatricesB2[i]->GetN_Entries(),
                  tau/oldtau,
                  MatricesB2[i]->GetEntries());
                Dscal(MatricesB3[i]->GetN_Entries(),
                  tau/oldtau,
                  MatricesB3[i]->GetEntries());
              }
              break;

            case 2:
            case 4:
            case 14:
              if (tau/oldtau != 1.0)
              {
                Dscal(MatricesB1T[i]->GetN_Entries(),
                  tau/oldtau,
                  MatricesB1T[i]->GetEntries());
                Dscal(MatricesB2T[i]->GetN_Entries(),
                  tau/oldtau,
                  MatricesB2T[i]->GetEntries());
                Dscal(MatricesB3T[i]->GetN_Entries(),
                  tau/oldtau,
                  MatricesB3T[i]->GetEntries());
              }
              // scale divergence constraint
              if ((TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT>0)
                &&(TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT*tau/oldtau != 1.0))
              {
                if (i==mg_level-1)
                  OutPut("scale " << tau/oldtau << endl);
                Dscal(MatricesB1[i]->GetN_Entries(),
                  TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT*tau/oldtau,
                  MatricesB1[i]->GetEntries());
                Dscal(MatricesB2[i]->GetN_Entries(),
                  TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT*tau/oldtau,
                  MatricesB2[i]->GetEntries());
                Dscal(MatricesB3[i]->GetN_Entries(),
                  TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT*tau/oldtau,
                  MatricesB3[i]->GetEntries());

                if (TDatabase::ParamDB->NSTYPE == 14)
                {
                  Dscal(MatricesC[i]->GetN_Entries(),
                    TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT*tau/oldtau,
                    MatricesC[i]->GetEntries());
                  Dscal(N_P, TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT*tau/oldtau, RHSs[3]);
                }
                if (i == mg_level -1)
                  TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT = 1.0;
              }
              break;
          }
        }                                         // endfor
        oldtau = tau;

        // update rhs by Laplacian and convective term from previous
        // time step
        // scaled by current sub time step length and theta2
        // currently : M := M + gamma A
        // M = M + (-gamma - tau*theta2) A
        // NOTE: if (methods==0) A is defined by u of the start time step
        //       if (methods==1) A is defined by u of the final time step
        //                       of the first method
        // ==> one does not get the same numbers if the methods are changed

        for(i=0;i<mg_level;i++)
        {
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
            case 2:
              MatAdd(MatricesM[i], MatricesA[i],
                -gamma - tau*theta2);
              break;

            case 3:
            case 4:
              MatAdd(MatricesM11[i], MatricesA11[i],
                -gamma - tau*theta2);
              MatAdd(MatricesM12[i], MatricesA12[i],
                -gamma - tau*theta2);
              MatAdd(MatricesM13[i], MatricesA13[i],
                -gamma - tau*theta2);
              MatAdd(MatricesM21[i], MatricesA21[i],
                -gamma - tau*theta2);
              MatAdd(MatricesM22[i], MatricesA22[i],
                -gamma - tau*theta2);
              MatAdd(MatricesM23[i], MatricesA23[i],
                -gamma - tau*theta2);
              MatAdd(MatricesM31[i], MatricesA31[i],
                -gamma - tau*theta2);
              MatAdd(MatricesM32[i], MatricesA32[i],
                -gamma - tau*theta2);
              MatAdd(MatricesM33[i], MatricesA33[i],
                -gamma - tau*theta2);
              // remove contributions from stabilizations
              if ((TDatabase::ParamDB->DISCTYPE == SDFEM)&&(!very_first_supg))
              {
                MatAdd(MatricesM11[i], MatricesS11[i], -tau * theta1);
                MatAdd(MatricesM12[i], MatricesS12[i], -tau * theta1);
                MatAdd(MatricesM13[i], MatricesS13[i], -tau * theta1);
                MatAdd(MatricesM21[i], MatricesS21[i], -tau * theta1);
                MatAdd(MatricesM22[i], MatricesS22[i], -tau * theta1);
                MatAdd(MatricesM23[i], MatricesS23[i], -tau * theta1);
                MatAdd(MatricesM31[i], MatricesS31[i], -tau * theta1);
                MatAdd(MatricesM32[i], MatricesS32[i], -tau * theta1);
                MatAdd(MatricesM33[i], MatricesS33[i], -tau * theta1);
                MatAdd(MatricesM11[i], MatricesK[i], -1);
                MatAdd(MatricesM22[i], MatricesK[i], -1);
                MatAdd(MatricesM33[i], MatricesK[i], -1);
              }
              break;
            case 14:
              MatAdd(MatricesM11[i], MatricesA11[i],
                -gamma - tau*theta2);
              MatAdd(MatricesM12[i], MatricesA12[i],
                -gamma - tau*theta2);
              MatAdd(MatricesM13[i], MatricesA13[i],
                -gamma - tau*theta2);
              MatAdd(MatricesM21[i], MatricesA21[i],
                -gamma - tau*theta2);
              MatAdd(MatricesM22[i], MatricesA22[i],
                -gamma - tau*theta2);
              MatAdd(MatricesM23[i], MatricesA23[i],
                -gamma - tau*theta2);
              MatAdd(MatricesM31[i], MatricesA31[i],
                -gamma - tau*theta2);
              MatAdd(MatricesM32[i], MatricesA32[i],
                -gamma - tau*theta2);
              MatAdd(MatricesM33[i], MatricesA33[i],
                -gamma - tau*theta2);
              break;
          }                                       // endswitch
        }                                         // end i
        // set current factor of steady state matrix
        gamma = -tau*theta2;
        very_first_supg = 0;

        // defect = M * sol
        // B:= B + defect
        switch(TDatabase::ParamDB->NSTYPE)
        {
          case 1:
          case 2:
            MatVectActive(MatricesM[mg_level-1], sol, defect);
            Daxpy(N_Active, 1, defect, B);
            MatVectActive(MatricesM[mg_level-1], sol+N_U, defect+N_U);
            Daxpy(N_Active, 1, defect+N_U, B+N_U);
            MatVectActive(MatricesM[mg_level-1], sol+2*N_U, defect+2*N_U);
            Daxpy(N_Active, 1, defect+2*N_U, B+2*N_U);
            break;

          case 3:
          case 4:
          case 14:
            MatVectActive(MatricesM11[mg_level-1], sol, defect);
            Daxpy(N_Active, 1, defect, B);
            MatVectActive(MatricesM12[mg_level-1], sol+N_U, defect);
            Daxpy(N_Active, 1, defect, B);
            MatVectActive(MatricesM13[mg_level-1], sol+2*N_U, defect);
            Daxpy(N_Active, 1, defect, B);
            MatVectActive(MatricesM21[mg_level-1], sol, defect+N_U);
            Daxpy(N_Active, 1, defect+N_U, B+N_U);
            MatVectActive(MatricesM22[mg_level-1], sol+N_U, defect+N_U);
            Daxpy(N_Active, 1, defect+N_U, B+N_U);
            MatVectActive(MatricesM23[mg_level-1], sol+2*N_U, defect+N_U);
            Daxpy(N_Active, 1, defect+N_U, B+N_U);
            MatVectActive(MatricesM31[mg_level-1], sol, defect+2*N_U);
            Daxpy(N_Active, 1, defect+2*N_U, B+2*N_U);
            MatVectActive(MatricesM32[mg_level-1], sol+N_U, defect+2*N_U);
            Daxpy(N_Active, 1, defect+2*N_U, B+2*N_U);
            MatVectActive(MatricesM33[mg_level-1], sol+2*N_U, defect+2*N_U);
            Daxpy(N_Active, 1, defect+2*N_U, B+2*N_U);
            break;
        }

        //========================================================================
        // end assembling of rhs
        //========================================================================

        // extrapolate solution to get starting value at next time step
        // equidistant time steps assumed
        // current solution sol
        // solution of last time step sol_timestep_m1

        // save sol
        memcpy(oldsol, sol, SizeOfDouble*N_Unknowns);
        if (TDatabase::TimeDB->EXTRAPOLATE_VELOCITY)
        {
          // sol := tau1 *sol - tau2 * sol_timestep_m1
          // at first time step: sol = sol_timestep_m1 -> result is sol
          tau2 = tau/oldtau;
          tau1 = 1 + tau2;
          for (k=0;k<3*N_U;k++)
            sol[k] = tau1*sol[k] - tau2*sol_timestep_m1[k];
          // save current solution
          memcpy(sol_timestep_m1, oldsol, SizeOfDouble*N_Unknowns);
        }

        // set Dirichlet values
        // RHSs[0] still available from assembling
        memcpy(B+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
        memcpy(B+N_Active+N_U, RHSs[1]+N_Active, (N_U-N_Active)*SizeOfDouble);
        memcpy(B+N_Active+2*N_U, RHSs[2]+N_Active, (N_U-N_Active)*SizeOfDouble);

        // copy Dirichlet values from rhs into sol
        memcpy(sol+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
        memcpy(sol+N_Active+N_U, RHSs[1]+N_Active, (N_U-N_Active)*SizeOfDouble);
        memcpy(sol+N_Active+2*N_U, RHSs[2]+N_Active, (N_U-N_Active)*SizeOfDouble);

        //========================================================================
        // assembling of system matrix
        //========================================================================

        // M = M + (-gamma + tau*theta1) A
        for(i=0;i<mg_level;i++)
        {
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
            case 2:
              MatAdd(MatricesM[i], MatricesA[i],
                -gamma + tau*theta1);
              break;

            case 3:
            case 4:
              MatAdd(MatricesM11[i], MatricesA11[i],
                -gamma + tau*theta1);
              MatAdd(MatricesM12[i], MatricesA12[i],
                -gamma + tau*theta1);
              MatAdd(MatricesM13[i], MatricesA13[i],
                -gamma + tau*theta1);
              MatAdd(MatricesM21[i], MatricesA21[i],
                -gamma + tau*theta1);
              MatAdd(MatricesM22[i], MatricesA22[i],
                -gamma + tau*theta1);
              MatAdd(MatricesM23[i], MatricesA23[i],
                -gamma + tau*theta1);
              MatAdd(MatricesM31[i], MatricesA31[i],
                -gamma + tau*theta1);
              MatAdd(MatricesM32[i], MatricesA32[i],
                -gamma + tau*theta1);
              MatAdd(MatricesM33[i], MatricesA33[i],
                -gamma + tau*theta1);
              if (TDatabase::ParamDB->DISCTYPE == SDFEM)
              {
                MatAdd(MatricesM11[i], MatricesS11[i], tau*theta1);
                MatAdd(MatricesM12[i], MatricesS12[i], tau*theta1);
                MatAdd(MatricesM13[i], MatricesS13[i], tau*theta1);
                MatAdd(MatricesM21[i], MatricesS21[i], tau*theta1);
                MatAdd(MatricesM22[i], MatricesS22[i], tau*theta1);
                MatAdd(MatricesM23[i], MatricesS23[i], tau*theta1);
                MatAdd(MatricesM31[i], MatricesS31[i], tau*theta1);
                MatAdd(MatricesM32[i], MatricesS32[i], tau*theta1);
                MatAdd(MatricesM33[i], MatricesS33[i], tau*theta1);
                MatAdd(MatricesM11[i], MatricesK[i], 1);
                MatAdd(MatricesM22[i], MatricesK[i], 1);
                MatAdd(MatricesM33[i], MatricesK[i], 1);
              }
              break;
            case 14:
              MatAdd(MatricesM11[i], MatricesA11[i],
                -gamma + tau*theta1);
              MatAdd(MatricesM12[i], MatricesA12[i],
                -gamma + tau*theta1);
              MatAdd(MatricesM13[i], MatricesA13[i],
                -gamma + tau*theta1);
              MatAdd(MatricesM21[i], MatricesA21[i],
                -gamma + tau*theta1);
              MatAdd(MatricesM22[i], MatricesA22[i],
                -gamma + tau*theta1);
              MatAdd(MatricesM23[i], MatricesA23[i],
                -gamma + tau*theta1);
              MatAdd(MatricesM31[i], MatricesA31[i],
                -gamma + tau*theta1);
              MatAdd(MatricesM32[i], MatricesA32[i],
                -gamma + tau*theta1);
              MatAdd(MatricesM33[i], MatricesA33[i],
                -gamma + tau*theta1);
              MatAdd(MatricesM11[i], MatricesK[i],1);
              MatAdd(MatricesM22[i], MatricesK[i],1);
              MatAdd(MatricesM33[i], MatricesK[i],1);
              break;
          }                                       // endswitch
        }
        // set current factor of steady state matrix
        gamma = tau*theta1;

        //========================================================================
        // end assembling of system matrix
        //========================================================================

        OutPut("CURRENT TIME: ");
        OutPut(TDatabase::TimeDB->CURRENTTIME << endl);
        OutPut("MEMORY: " << setw(10) << GetMemory()/(1048576.0));
        OutPut(" MB" << endl);
        if ((hmin*tau<=5e-5)&&((TDatabase::ParamDB->NSTYPE==1)||(TDatabase::ParamDB->NSTYPE==3)))
        {
          OutPut("ATTENTION: divergence constraint will be multiplied with small time step !!!"<<endl);
          OutPut("remove the exit if this is correct " << endl);
          exit(4711);
        }

        //======================================================================
        // nonlinear loop
        //======================================================================

        N_LinIterCurr = 0;
        solver_time_curr = 0;
        for(j=0;j<=Max_It;j++)                    // solve nonlinear equation
        {
          memcpy(oldsol, sol, SizeOfDouble*N_Unknowns);
          memset(defect, 0, N_Unknowns*SizeOfDouble);

          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
              SQMATRICES[0] = MatricesM[mg_level-1];
              MATRICES[0] = MatricesB1[mg_level-1];
              MATRICES[1] = MatricesB2[mg_level-1];
              MATRICES[2] = MatricesB3[mg_level-1];
              break;
            case 2:
              SQMATRICES[0] = MatricesM[mg_level-1];
              MATRICES[0] = MatricesB1[mg_level-1];
              MATRICES[1] = MatricesB2[mg_level-1];
              MATRICES[2] = MatricesB3[mg_level-1];
              MATRICES[3] = MatricesB1T[mg_level-1];
              MATRICES[4] = MatricesB2T[mg_level-1];
              MATRICES[5] = MatricesB3T[mg_level-1];
              break;
            case 3:
              SQMATRICES[0] = MatricesM11[mg_level-1];
              SQMATRICES[1] = MatricesM12[mg_level-1];
              SQMATRICES[2] = MatricesM13[mg_level-1];
              SQMATRICES[3] = MatricesM21[mg_level-1];
              SQMATRICES[4] = MatricesM22[mg_level-1];
              SQMATRICES[5] = MatricesM23[mg_level-1];
              SQMATRICES[6] = MatricesM31[mg_level-1];
              SQMATRICES[7] = MatricesM32[mg_level-1];
              SQMATRICES[8] = MatricesM33[mg_level-1];
              MATRICES[0] = MatricesB1[mg_level-1];
              MATRICES[1] = MatricesB2[mg_level-1];
              MATRICES[2] = MatricesB3[mg_level-1];
              break;
            case 4:
            case 14:
              SQMATRICES[0] = MatricesM11[mg_level-1];
              SQMATRICES[1] = MatricesM12[mg_level-1];
              SQMATRICES[2] = MatricesM13[mg_level-1];
              SQMATRICES[3] = MatricesM21[mg_level-1];
              SQMATRICES[4] = MatricesM22[mg_level-1];
              SQMATRICES[5] = MatricesM23[mg_level-1];
              SQMATRICES[6] = MatricesM31[mg_level-1];
              SQMATRICES[7] = MatricesM32[mg_level-1];
              SQMATRICES[8] = MatricesM33[mg_level-1];
              MATRICES[0] = MatricesB1[mg_level-1];
              MATRICES[1] = MatricesB2[mg_level-1];
              MATRICES[2] = MatricesB3[mg_level-1];
              MATRICES[3] = MatricesB1T[mg_level-1];
              MATRICES[4] = MatricesB2T[mg_level-1];
              MATRICES[5] = MatricesB3T[mg_level-1];
              if (TDatabase::ParamDB->NSTYPE == 14)
              {
                SQMATRICES[9] = MatricesC[mg_level-1];
              }
              break;
          }

          //#ifdef __CALOTTE__
          //	  SetDirichletNodesFromNeumannNodes(SQMATRICES, MATRICES+3, B, N_U,
          //					    N_neum_to_diri, neum_to_diri);
          //#endif
          // compute defect
          OutPut("B " << Ddot(N_Unknowns,sol,sol) << " "
            <<  Ddot(N_Unknowns,B,B) << endl);
          Defect(sqmatrices,matrices,sol,B,defect);
          // scale divergence constraint appropriately
          //if ((TDatabase::ParamDB->NSTYPE==1) || (TDatabase::ParamDB->NSTYPE==3))
          //    Dscal(N_P, 1.0/tau, sol+3*N_U);

          if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
            IntoL20Vector3D(defect+3*N_U, N_P,pressure_space_code);
          residual =  Ddot(N_Unknowns, defect, defect);
          impuls_residual = Ddot(3*N_U, defect, defect);
          OutPut("nonlinear step " << setw(3) << j);
          OutPut(setw(14) << impuls_residual);
          OutPut(setw(14) << Ddot(N_P,defect+3*N_U,defect+3*N_U));
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
                  Solver(sqmatrixM, matrixB1, matrixB2,  matrixB3, B, sol);
                  break;

                default :
                  OutPut("AMG not implemented" << endl);
                  exit(4711);
                  break;
              }
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
              solver_time += t2-t1;

              // update
              for(k=0;k<N_Unknowns;k++)
              {
                p2 = sol[k]-oldsol[k];
                sol[k] = oldsol[k] + omega * p2;
              }
              break;
          }                                       // endswitch SOLVER_TYPE
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
                MatAdd(MatricesM13[i], MatricesA13[i], -gamma);
                MatAdd(MatricesM21[i], MatricesA21[i], -gamma);
                MatAdd(MatricesM22[i], MatricesA22[i], -gamma);
                MatAdd(MatricesM23[i], MatricesA23[i], -gamma);
                MatAdd(MatricesM31[i], MatricesA31[i], -gamma);
                MatAdd(MatricesM32[i], MatricesA32[i], -gamma);
                MatAdd(MatricesM33[i], MatricesA33[i], -gamma);
                if (TDatabase::ParamDB->DISCTYPE == SDFEM)
                {
                  MatAdd(MatricesM11[i], MatricesS11[i], -tau*theta1);
                  MatAdd(MatricesM12[i], MatricesS12[i], -tau*theta1);
                  MatAdd(MatricesM13[i], MatricesS13[i], -tau*theta1);
                  MatAdd(MatricesM21[i], MatricesS21[i], -tau*theta1);
                  MatAdd(MatricesM22[i], MatricesS22[i], -tau*theta1);
                  MatAdd(MatricesM23[i], MatricesS23[i], -tau*theta1);
                  MatAdd(MatricesM31[i], MatricesS31[i], -tau*theta1);
                  MatAdd(MatricesM32[i], MatricesS32[i], -tau*theta1);
                  MatAdd(MatricesM33[i], MatricesS33[i], -tau*theta1);
                  MatAdd(MatricesM11[i], MatricesK[i], -1);
                  MatAdd(MatricesM22[i], MatricesK[i], -1);
                  MatAdd(MatricesM33[i], MatricesK[i], -1);
                }
                break;
              case 14:
                MatAdd(MatricesM11[i], MatricesA11[i], -gamma);
                MatAdd(MatricesM12[i], MatricesA12[i], -gamma);
                MatAdd(MatricesM13[i], MatricesA13[i], -gamma);
                MatAdd(MatricesM21[i], MatricesA21[i], -gamma);
                MatAdd(MatricesM22[i], MatricesA22[i], -gamma);
                MatAdd(MatricesM23[i], MatricesA23[i], -gamma);
                MatAdd(MatricesM31[i], MatricesA31[i], -gamma);
                MatAdd(MatricesM32[i], MatricesA32[i], -gamma);
                MatAdd(MatricesM33[i], MatricesA33[i], -gamma);
                MatAdd(MatricesM11[i], MatricesK[i], -1);
                MatAdd(MatricesM22[i], MatricesK[i], -1);
                MatAdd(MatricesM33[i], MatricesK[i], -1);
                break;
            }                                     // endswitch
          }                                       // endfor i
          // set current factor of steady state matrix
          gamma = 0;

          if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
            MG->RestrictToAllGrids();
          // subtract term from Newton's method from the rhs array B
          if (nonlinite==1)
          {
            Daxpy(N_Active, -tau*theta1, newton, B);
            Daxpy(N_Active, -tau*theta1, newton+N_U, B+N_U);
            Daxpy(N_Active, -tau*theta1, newton+2*N_U, B+2*N_U);
          }
          //======================================================================
          // assemble new matrix due to nonlinearity
          //======================================================================
          t1 = GetTime();
          TDatabase::ParamDB->INTERNAL_LEVEL = 0;
          for(i=0;i<mg_level;i++)
          {
            if (i==mg_level-1)
              TDatabase::ParamDB->INTERNAL_LEVEL = 1;
            if (((mg_type==1)&&(i<mg_level-1))||((mg_type==2)&&(i<mg_level-2)))
            {
              DiscreteForm = DiscreteFormNLUpwindNC;
              CurrentDiscType =  UPWIND;
            }
            else
              switch(TDatabase::ParamDB->DISCTYPE)
              {
                case GALERKIN:
                  DiscreteForm = DiscreteFormNLGalerkin;
                  CurrentDiscType =  GALERKIN;
                  break;
                case UPWIND:
                  DiscreteForm = DiscreteFormNLUpwind;
                  CurrentDiscType =  UPWIND;
                  break;
                case SMAGORINSKY:
              case VMS_PROJECTION_EXPL:
                DiscreteForm = DiscreteFormNLSmagorinsky;
                CurrentDiscType =  SMAGORINSKY;
                break;
              case  CLASSICAL_LES:
                DiscreteForm = DiscreteFormNLClassicalLES;
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
                DiscreteForm = DiscreteFormNLVMS_Projection;
                CurrentDiscType =  VMS_PROJECTION;
                break;
              case VMS_RFB_EXPL:
              case VMS_RFB_EXPL_COUPLED:
                DiscreteForm = DiscreteFormNLVMS_RFBExplRhs;
                CurrentDiscType = VMS_RFB_EXPL;
                break;
              case SDFEM:
                DiscreteForm = DiscreteFormNLVMS_SUPG;
                CurrentDiscType =  SDFEM;
                break;
              default:
                OutPut("Unknown DISCTYPE " << TDatabase::ParamDB->DISCTYPE << endl);
                exit(1);
            }
            if (DiscreteForm==NULL)
            {
              OutPut("DiscreteForm not implemented !!!"<< endl);
              exit(4711);
            }

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
              case 14:
                N_RectMatrices = 0;
                N_Rhs = 0;
                N_FESpaces = 1;

                if (((TDatabase::ParamDB->LAPLACETYPE==1)
                  &&
                  ((CurrentDiscType == SMAGORINSKY) ||
                  (CurrentDiscType == CLASSICAL_LES) ||
                  (CurrentDiscType == GL00_CONVOLUTION) ||
                  (CurrentDiscType == GL00_AUX_PROBLEM) ||
                  (CurrentDiscType == VMS_PROJECTION)||
                  (CurrentDiscType == VMS_PROJECTION_EXPL)) ||
                  (CurrentDiscType == VMS_RFB_EXPL) ||
                  (CurrentDiscType == VMS_RFB_EXPL_COUPLED)) ||
                  (CurrentDiscType == SDFEM) ||
                  (nonlinite==1))
                {
                  SQMATRICES[0] = MatricesA11[i];
                  SQMATRICES[1] = MatricesA12[i];
                  SQMATRICES[2] = MatricesA13[i];
                  SQMATRICES[3] = MatricesA21[i];
                  SQMATRICES[4] = MatricesA22[i];
                  SQMATRICES[5] = MatricesA23[i];
                  SQMATRICES[6] = MatricesA31[i];
                  SQMATRICES[7] = MatricesA32[i];
                  SQMATRICES[8] = MatricesA33[i];
                  SQMATRICES[0]->Reset();
                  SQMATRICES[1]->Reset();
                  SQMATRICES[2]->Reset();
                  SQMATRICES[3]->Reset();
                  SQMATRICES[4]->Reset();
                  SQMATRICES[5]->Reset();
                  SQMATRICES[6]->Reset();
                  SQMATRICES[7]->Reset();
                  SQMATRICES[8]->Reset();
                  N_SquareMatrices = 9;
                  mid_sq = 4;
                  last_sq = 8;
                  if (CurrentDiscType == VMS_PROJECTION)
                  {
                    switch (TDatabase::ParamDB->VMS_LARGE_VELOCITY_SPACE)
                    {
                      case 17:
                        N_SquareMatrices += 1;
                        SQMATRICES[N_SquareMatrices - 1] = MatricesL[i];
                        SQMATRICES[N_SquareMatrices - 1]->Reset();
                        N_RectMatrices = 6;
                        MATRICES[0] = Matrices_tilde_G11[i];
                        MATRICES[1] = Matrices_tilde_G22[i];
                        MATRICES[2] = Matrices_tilde_G33[i];
                        MATRICES[3] = Matrices_G11[i];
                        MATRICES[4] = Matrices_G22[i];
                        MATRICES[5] = Matrices_G33[i];
                        MATRICES[0]->Reset();
                        MATRICES[1]->Reset();
                        MATRICES[2]->Reset();
                        MATRICES[3]->Reset();
                        MATRICES[4]->Reset();
                        MATRICES[5]->Reset();
                        break;
                      default:
                        N_RectMatrices = 3;
                        MATRICES[0] = Matrices_tilde_G11[i];
                        MATRICES[1] = Matrices_tilde_G22[i];
                        MATRICES[2] = Matrices_tilde_G33[i];
                        MATRICES[0]->Reset();
                        MATRICES[1]->Reset();
                        MATRICES[2]->Reset();
                        break;
                    }
                    N_FESpaces = 3;
                    fesp[2] = ProjectionSpaces[i];
                  }
                  if (CurrentDiscType == SDFEM)
                  {
                    N_FESpaces = 2;
                    fesp[1] = PSpaces[i];
                    N_SquareMatrices = 22;
                    SQMATRICES[N_SquareMatrices-13] =  MatricesM11[i];
                    //SQMATRICES[N_SquareMatrices-5]->Reset();
                    SQMATRICES[N_SquareMatrices-12] =  MatricesM22[i];
                    //SQMATRICES[N_SquareMatrices-4]->Reset();
                    SQMATRICES[N_SquareMatrices-11] =  MatricesM33[i];
                    //SQMATRICES[N_SquareMatrices-4]->Reset();
                    SQMATRICES[N_SquareMatrices-10] =  MatricesK[i];
                    SQMATRICES[N_SquareMatrices-10]->Reset();
                    SQMATRICES[N_SquareMatrices-9] =  MatricesS11[i];
                    SQMATRICES[N_SquareMatrices-9]->Reset();
                    SQMATRICES[N_SquareMatrices-8] =  MatricesS12[i];
                    SQMATRICES[N_SquareMatrices-8]->Reset();
                    SQMATRICES[N_SquareMatrices-7] =  MatricesS13[i];
                    SQMATRICES[N_SquareMatrices-7]->Reset();
                    SQMATRICES[N_SquareMatrices-6] =  MatricesS21[i];
                    SQMATRICES[N_SquareMatrices-6]->Reset();
                    SQMATRICES[N_SquareMatrices-5] =  MatricesS22[i];
                    SQMATRICES[N_SquareMatrices-5]->Reset();
                    SQMATRICES[N_SquareMatrices-4] =  MatricesS23[i];
                    SQMATRICES[N_SquareMatrices-4]->Reset();
                    SQMATRICES[N_SquareMatrices-3] =  MatricesS31[i];
                    SQMATRICES[N_SquareMatrices-3]->Reset();
                    SQMATRICES[N_SquareMatrices-2] =  MatricesS32[i];
                    SQMATRICES[N_SquareMatrices-2]->Reset();
                    SQMATRICES[N_SquareMatrices-1] =  MatricesS33[i];
                    SQMATRICES[N_SquareMatrices-1]->Reset();
                    N_RectMatrices = 0;
                    /*MATRICES[0] = MatricesB1[i];
                    MATRICES[1] = MatricesB2[i];
                    MATRICES[2] = MatricesB3[i];
                    MATRICES[3] = MatricesB1T[i];
                    MATRICES[4] = MatricesB2T[i];
                    MATRICES[5] = MatricesB3T[i];
                    MATRICES[0]->Reset();
                    MATRICES[1]->Reset();
                    MATRICES[2]->Reset();
                    MATRICES[3]->Reset();
                    MATRICES[4]->Reset();
                    MATRICES[5]->Reset();
                    */

                    if (TDatabase::ParamDB->NSTYPE == 14)
                    {
                      // pressure coupling does not need update since it is linear
                      N_SquareMatrices += 1;
                      SQMATRICES[N_SquareMatrices-1] =  MatricesK[i];
                      SQMATRICES[N_SquareMatrices-1]->Reset();
                    }
                  }
                }
                else
                {
                  SQMATRICES[0] = MatricesA11[i];
                  SQMATRICES[1] = MatricesA22[i];
                  SQMATRICES[2] = MatricesA33[i];

                  SQMATRICES[0]->Reset();
                  SQMATRICES[1]->Reset();
                  SQMATRICES[2]->Reset();

                  N_SquareMatrices = 3;
                  mid_sq = 1;
                  last_sq = 2;
                }
                break;
            }

            fesp[0] = USpaces[i];
            fesp[1] = PSpaces[i];

            fefct[0] = U1Array[i];
            fefct[1] = U2Array[i];
            fefct[2] = U3Array[i];

            ferhs[0] = USpaces[i];
            ferhs[1] = USpaces[i];
            ferhs[2] = USpaces[i];

            // Newton's method
            if (nonlinite==1)
            {
              aux =  new TAuxParam3D(TimeNSN_FESpacesVelo_GradVelo,
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
              {
                OutPut("Newton's method in connection with TURBULENT_VISCOSITY_TYPE == 4 not implemented !!!"
                  << endl);
                exit(4711);
              }
            }
            else
            {
              switch(TDatabase::ParamDB->DISCTYPE)
              {
                // turbulent viscosity must be computed
                case SMAGORINSKY:
                case CLASSICAL_LES:
                case GL00_CONVOLUTION:
                case GL00_AUX_PROBLEM:
                  // Smagorinsky matrix will be assembled
                case VMS_PROJECTION_EXPL:
                case SDFEM:
                  // standard turbulent viscosities on coarser grids
                  if ((TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE!=4)
                    || (i<mg_level -1))
                  {
                    aux =  new TAuxParam3D(TimeNSN_FESpacesVelo_GradVelo,
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
                    && (i == mg_level -1))
                  {
                    fesp[1] = uConvSpaces[mg_level-1];
                    N_FESpaces++;
                    fefct[3] = u1ConvArray[mg_level-1];
                    fefct[4] = u2ConvArray[mg_level-1];
                    fefct[5] = u3ConvArray[mg_level-1];
                    //OutPut("aa " << u1ConvArray[mg_level-1]->GetValues()[0] << endl);

                    aux =  new TAuxParam3D(TimeNSN_FESpacesVelo_GradVelo_ConvVelo,
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
                case UPWIND:
                  aux = new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
                  break;
                case VMS_PROJECTION:
                  // standard turbulent viscosities on coarser grids
                  // large scale turbulent viscosity
                  /*  if ((i<mg_level-1)&&(TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE!=5))
                    {
                  OutPut("cc" << endl);
                  N_FESpaces = 5;
                  fesp[4] = label_space_fesp;
                  fefct[9] = label_space_fefct;
                  // turbulent viscosity with u^h
                  ii = TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE;
                  TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE = 1;
                  aux =  new TAuxParam3D(TimeNSN_FESpacesVelo_GradVelo,
                             TimeNSN_FctVelo_GradVelo,
                  TimeNSN_ParamFctVelo_GradVelo,
                  TimeNSN_FEValuesVelo_GradVelo,
                  fesp, fefct,
                  TimeNSFctVelo_GradVelo,
                  TimeNSFEFctIndexVelo_GradVelo,
                  TimeNSFEMultiIndexVelo_GradVelo,
                  TimeNSN_ParamsVelo_GradVelo,
                  TimeNSBeginParamVelo_GradVelo);
                  TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE = ii;
                  }
                  else
                  // finest grid
                  {*/
                  if (i== mg_level -1)
                  {
                    N_FESpaces = 4;
                    fesp[3] = label_space_fesp;
                    TDatabase::ParamDB->INTERNAL_LEVEL = 1;
                    fefct[3] = vms_proj_11;
                    fefct[4] = vms_proj_12;
                    fefct[5] = vms_proj_13;
                    fefct[6] = vms_proj_22;
                    fefct[7] = vms_proj_23;
                    fefct[8] = vms_proj_33;
                    fefct[9] = label_space_fefct;
                    // compute large scales
                    ComputeVMSProjection(Matrices_G11[mg_level-1], Matrices_G22[mg_level-1],
                      Matrices_G33[mg_level-1], MatricesL[mg_level-1],
                      U1Array[mg_level-1], U2Array[mg_level-1],
                      U3Array[mg_level-1], vms_projection_fe);
                  }
                  else
                  {                               // just dummies
                    // turbulent viscosity with u^h
                    ii = TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE;
                    TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE = 1;
                    TDatabase::ParamDB->INTERNAL_LEVEL = 0;
                    fefct[3] = u1;
                    fefct[4] = u1;
                    fefct[5] = u1;
                    fefct[6] = u1;
                    fefct[7] = u1;
                    fefct[8] = u1;
                    fefct[9] = u1;
                    TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE = ii;
                  }

                  // turbulent viscosity for the large scales
                  aux =  new TAuxParam3D(TimeNSN_FESpacesVelo_GradVelo_LargeScale,
                    TimeNSN_FctVelo_GradVelo_LargeScale,
                    TimeNSN_ParamFctVelo_GradVelo_LargeScale,
                    TimeNSN_FEValuesVelo_GradVelo_LargeScale,
                    fesp, fefct,
                    TimeNSFctVelo_GradVelo_LargeScale,
                    TimeNSFEFctIndexVelo_GradVelo_LargeScale,
                    TimeNSFEMultiIndexVelo_GradVelo_LargeScale,
                    TimeNSN_ParamsVelo_GradVelo_LargeScale,
                    TimeNSBeginParamVelo_GradVelo_LargeScale);
                  //}

                  break;

                default:
                  aux =  new TAuxParam3D(TimeNSN_FESpacesVelo, TimeNSN_FctVelo,
                    TimeNSN_ParamFctVelo,
                    TimeNSN_FEValuesVelo,
                    fesp, fefct,
                    TimeNSFctVelo,
                    TimeNSFEFctIndexVelo, TimeNSFEMultiIndexVelo,
                    TimeNSN_ParamsVelo, TimeNSBeginParamVelo);
                  break;
              }
            }
            if (i<mg_level-1)
            {
              // more viscosity on coarser levels, for solver
              cell_measure = TDatabase::ParamDB->CELL_MEASURE;
              TDatabase::ParamDB->CELL_MEASURE = 0;
            }
            //======================================================================
            // assembling of matrices for each level due to nonlinearity
            // A_11, ...
            // no assembling of rhs
            //======================================================================
            Assemble3D(N_FESpaces, fesp,
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

            if (i<mg_level-1)
            {
              // more viscosity on coarser levels, for solver
              TDatabase::ParamDB->CELL_MEASURE = cell_measure;
            }

            if ((DiscreteForm == DiscreteFormNLUpwind)||(DiscreteForm == DiscreteFormNLUpwindNC))
            {
              switch(TDatabase::ParamDB->NSTYPE)
              {
                case 1:
                case 2:
                  // do upwinding with one matrix
                  UpwindForNavierStokes3D(SQMATRICES[0], U1Array[i], U2Array[i], U3Array[i]);
                  //cout << "UPWINDING DONE : level " << i << endl;
                  break;

                case 3:
                case 4:
                case 14:
                  // do upwinding with three matrices
                  UpwindForNavierStokes3D(SQMATRICES[0], U1Array[i], U2Array[i], U3Array[i]);
                  UpwindForNavierStokes3D(SQMATRICES[mid_sq], U1Array[i], U2Array[i], U3Array[i]);
                  UpwindForNavierStokes3D(SQMATRICES[last_sq], U1Array[i], U2Array[i], U3Array[i]);
                  //cout << "UPWINDING DONE(2) : level " << i << endl;
                  //cout << "check correct sqmatrix !!!! " << endl;
                  break;
              }                                   // endswitch
            }

            // slip type bc detected, modify matrices accordingly
            if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION == 1)
            {
              // prepare everything for the assembling of slip with friction bc
              // on level i
              N_FESpaces = 1;
              N_SquareMatrices = 3;
              N_RectMatrices = 0;
              N_Rhs = 3;
              DiscreteForm = NULL;

              SQMATRICES[0] = MatricesA11[i];
              SQMATRICES[1] = MatricesA22[i];
              SQMATRICES[2] = MatricesA33[i];

              fesp[0] = USpaces[i];
              ferhs[0] = USpaces[i];
              ferhs[1] = USpaces[i];
              ferhs[2] = USpaces[i];

              RHSs[0] = RhsArray[i];
              RHSs[1] = RhsArray[i]+N_Uarray[i];
              RHSs[2] = RhsArray[i]+2*N_Uarray[i];

              Assemble3DSlipBC(N_FESpaces, fesp,
                N_SquareMatrices, SQMATRICES,
                N_RectMatrices, MATRICES,
                N_Rhs, RHSs, ferhs,
                DiscreteForm,
                BoundaryConditions,
                BoundValues,
                aux);
            }

            delete aux;

            if (CurrentDiscType == VMS_PROJECTION)
            {
              LumpMassMatrixToDiag(MatricesL[i]);
              SQMATRICES[0] = MatricesA11[i];
              SQMATRICES[1] = MatricesA12[i];
              SQMATRICES[2] = MatricesA13[i];
              SQMATRICES[3] = MatricesA21[i];
              SQMATRICES[4] = MatricesA22[i];
              SQMATRICES[5] = MatricesA23[i];
              SQMATRICES[6] = MatricesA31[i];
              SQMATRICES[7] = MatricesA32[i];
              SQMATRICES[8] = MatricesA33[i];
              SQMATRICES[9] =  MatricesL[i];
              MATRICES[0] = Matrices_tilde_G11[i];
              MATRICES[1] = Matrices_tilde_G22[i];
              MATRICES[2] = Matrices_tilde_G33[i];
              MATRICES[3] = Matrices_G11[i];
              MATRICES[4] = Matrices_G22[i];
              MATRICES[5] = Matrices_G33[i];

              //t5 =  GetTime();
              if ((i==mg_level - 1)||(!TDatabase::ParamDB->VMS_COARSE_MG_SMAGO))
              {
                //OutPut("update ");
                VMS_ProjectionUpdateMatrices(N_Uarray[i], USpaces[i]->GetActiveBound(),
                  ProjectionSpaces[i]->GetN_DegreesOfFreedom(),
                  SQMATRICES,MATRICES);
                OutPut("update done"<<endl);
              }
              //OutPut("update VMS matrices "<<  GetTime() - t5 << endl);
            }

            //======================================================================
            // end of assemble new matrix due to nonlinearity
            //======================================================================

            // build stiffness matrix for next nonlinear iteration step
            // stiffness matrix (left upper block) is stored on
            // M11, (M12, M21, M22)
            // M = M +  tau*theta1 A
            switch(TDatabase::ParamDB->NSTYPE)
            {
              case 1:
              case 2:
                MatAdd(MatricesM[i], MatricesA[i],
                  tau*theta1);
                break;

              case 3:
              case 4:
                MatAdd(MatricesM11[i], MatricesA11[i],
                  tau*theta1);
                MatAdd(MatricesM12[i], MatricesA12[i],
                  tau*theta1);
                MatAdd(MatricesM13[i], MatricesA13[i],
                  tau*theta1);
                MatAdd(MatricesM21[i], MatricesA21[i],
                  tau*theta1);
                MatAdd(MatricesM22[i], MatricesA22[i],
                  tau*theta1);
                MatAdd(MatricesM23[i], MatricesA23[i],
                  tau*theta1);
                MatAdd(MatricesM31[i], MatricesA31[i],
                  tau*theta1);
                MatAdd(MatricesM32[i], MatricesA32[i],
                  tau*theta1);
                MatAdd(MatricesM33[i], MatricesA33[i],
                  tau*theta1);
                if (TDatabase::ParamDB->DISCTYPE == SDFEM)
                {
                  MatAdd(MatricesM11[i], MatricesS11[i], tau*theta1);
                  MatAdd(MatricesM12[i], MatricesS12[i], tau*theta1);
                  MatAdd(MatricesM13[i], MatricesS13[i], tau*theta1);
                  MatAdd(MatricesM21[i], MatricesS21[i], tau*theta1);
                  MatAdd(MatricesM22[i], MatricesS22[i], tau*theta1);
                  MatAdd(MatricesM23[i], MatricesS23[i], tau*theta1);
                  MatAdd(MatricesM31[i], MatricesS31[i], tau*theta1);
                  MatAdd(MatricesM32[i], MatricesS32[i], tau*theta1);
                  MatAdd(MatricesM33[i], MatricesS33[i], tau*theta1);
                  MatAdd(MatricesM11[i], MatricesK[i], 1);
                  MatAdd(MatricesM22[i], MatricesK[i], 1);
                  MatAdd(MatricesM33[i], MatricesK[i], 1);
                }
                break;
              case 14:
                MatAdd(MatricesM11[i], MatricesA11[i],
                  tau*theta1);
                MatAdd(MatricesM12[i], MatricesA12[i],
                  tau*theta1);
                MatAdd(MatricesM13[i], MatricesA13[i],
                  tau*theta1);
                MatAdd(MatricesM21[i], MatricesA21[i],
                  tau*theta1);
                MatAdd(MatricesM22[i], MatricesA22[i],
                  tau*theta1);
                MatAdd(MatricesM23[i], MatricesA23[i],
                  tau*theta1);
                MatAdd(MatricesM31[i], MatricesA31[i],
                  tau*theta1);
                MatAdd(MatricesM32[i], MatricesA32[i],
                  tau*theta1);
                MatAdd(MatricesM33[i], MatricesA33[i],
                  tau*theta1);
                MatAdd(MatricesM11[i], MatricesK[i], 1);
                MatAdd(MatricesM22[i], MatricesK[i], 1);
                MatAdd(MatricesM33[i], MatricesK[i], 1);
                Dscal(MatricesB1T[i]->GetN_Entries(),
                  tau,
                  MatricesB1T[i]->GetEntries());
                Dscal(MatricesB2T[i]->GetN_Entries(),
                  tau,
                  MatricesB2T[i]->GetEntries());
                Dscal(MatricesB3T[i]->GetN_Entries(),
                  tau,
                  MatricesB3T[i]->GetEntries());
                Dscal(MatricesB1[i]->GetN_Entries(),
                  tau,
                  MatricesB1[i]->GetEntries());
                Dscal(MatricesB2[i]->GetN_Entries(),
                  tau,
                  MatricesB2[i]->GetEntries());
                Dscal(MatricesB3[i]->GetN_Entries(),
                  tau,
                  MatricesB3[i]->GetEntries());

                Dscal(MatricesC[i]->GetN_Entries(),
                  tau,
                  MatricesC[i]->GetEntries());
                break;
            }
          }                                       // endfor i
          // Newton's method
          // correction of rhs only on the finest level
          if (nonlinite==1)
          {
            N_FESpaces = 1;
            fesp[0] = USpaces[mg_level-1];
            N_Rhs = 3;
            RHSs[3] = newton;
            RHSs[4] = newton+N_Uarray[mg_level-1];
            RHSs[5] = newton+2*N_Uarray[mg_level-1];
            memset(newton, 0, 3*N_Uarray[mg_level-1]*SizeOfDouble);
            ferhs[0] = USpaces[mg_level-1];
            ferhs[1] = USpaces[mg_level-1];
            ferhs[2] = USpaces[mg_level-1];
            N_SquareMatrices = 0;
            N_RectMatrices = 0;
            fefct[0] = U1Array[mg_level-1];
            fefct[1] = U2Array[mg_level-1];
            fefct[2] = U3Array[mg_level-1];

            DiscreteForm = DiscreteFormRHSNewtonNL;
            aux =  new TAuxParam3D(TimeNSN_FESpacesVelo_GradVelo,
              TimeNSN_FctVelo_GradVelo,
              TimeNSN_ParamFctVelo_GradVelo,
              TimeNSN_FEValuesVelo_GradVelo,
              fesp, fefct,
              TimeNSFctVelo_GradVelo,
              TimeNSFEFctIndexVelo_GradVelo,
              TimeNSFEMultiIndexVelo_GradVelo,
              TimeNSN_ParamsVelo_GradVelo,
              TimeNSBeginParamVelo_GradVelo);

            Assemble3D(N_FESpaces, fesp,
              N_SquareMatrices, SQMATRICES,
              N_RectMatrices, MATRICES,
              N_Rhs, RHSs+3, ferhs,
              DiscreteForm,
              BoundaryConditions,
              BoundValues,
              aux);

            delete aux;

            Daxpy(N_Active, tau*theta1, newton, B);
            Daxpy(N_Active, tau*theta1, newton+N_U, B+N_U);
            Daxpy(N_Active, tau*theta1, newton+2*N_U, B+2*N_U);
          }

          t2 = GetTime();
          ass_time += t2-t1;
          OutPut("time for assembling " << t2-t1 << "s " << ass_time  << "s "<<  endl);
          // set current factor of steady state matrix
          gamma = tau*theta1;
        }                                         // endfor Max_It (solution of nonlinear equation)

        //======================================================================
        // end of nonlinear loop
        //======================================================================

        if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
          IntoL20FEFunction3D(PArray[mg_level-1]->GetValues(),
            N_Parray[mg_level-1], PSpaces[mg_level-1]);

        t22 = GetTime();
        OutPut( "time for timestep: " << t22-t11 << "s"<< endl);

      }                                           // endfor l (sub steps of fractional step theta)
    }                                             // end of the two disc schemes (methods)

    if (time_discs==2)
    {
      // compute difference of solutions
      for (i=0;i<N_Unknowns;i++)
        sol[i]-=frac_step_sol[i];

      // compute norms of difference
      fesp[0] = USpaces[mg_level-1];
      fefct[0] = U1Array[mg_level-1];
      fefct[1] = U2Array[mg_level-1];
      fefct[2] = U3Array[mg_level-1];

      aux =  new TAuxParam3D(TimeNSN_FESpacesVelo, TimeNSN_FctVelo,
        TimeNSN_ParamFctVelo,
        TimeNSN_FEValuesVelo,
        fesp, fefct,
        TimeNSFctVelo,
        TimeNSFEFctIndexVelo, TimeNSFEMultiIndexVelo,
        TimeNSN_ParamsVelo, TimeNSBeginParamVelo);
      // errors
      U1Array[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors);

      U2Array[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors+2);

      U3Array[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors+4);

      PArray[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, PSpaces+mg_level-1, errors+6);
      // compute L^2 error
      errors[8] = sqrt(errors[0]*errors[0]+errors[2]*errors[2]+errors[4]*errors[4]
        +errors[6]*errors[6]);
      if (TDatabase::TimeDB->CURRENTTIME< end_time)
        ComputeNewTimeStep(errors[8]);
      // copy solution of fract. step scheme
      memcpy(sol,frac_step_sol,N_Unknowns*SizeOfDouble);

      delete aux;
    }                                             // adaptive time step control

    /**************************************************************************/
    //
    // the solution in the current discrete time is computed
    //
    /**************************************************************************/
    if (TDatabase::ParamDB->MEASURE_ERRORS)
    {
      fesp[0] = USpaces[mg_level-1];
      fefct[0] = U1Array[mg_level-1];
      fefct[1] = U2Array[mg_level-1];
      fefct[2] = U3Array[mg_level-1];

      aux =  new TAuxParam3D(TimeNSN_FESpacesVelo, TimeNSN_FctVelo,
        TimeNSN_ParamFctVelo,
        TimeNSN_FEValuesVelo,
        fesp, fefct,
        TimeNSFctVelo,
        TimeNSFEFctIndexVelo, TimeNSFEMultiIndexVelo,
        TimeNSN_ParamsVelo, TimeNSBeginParamVelo);

      // errors
      // L2: error[0], H1-semi: error[1]
      U1Array[mg_level-1]->GetErrors(ExactU1, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors);

      // L2: error[2], H1-semi: error[3]
      U2Array[mg_level-1]->GetErrors(ExactU2, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors+2);

      // L2: error[4], H1-semi: error[5]
      U3Array[mg_level-1]->GetErrors(ExactU3, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors+4);

      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "L2(u): " << sqrt(errors[0]*errors[0]+errors[2]*errors[2]
        +errors[4]*errors[4]));
      OutPut( "   H1-semi(u):  " << sqrt(errors[1]*errors[1]+errors[3]*errors[3]
        +errors[5]*errors[5])<<endl);

      // error in L^infty(0,t,L^2)
      if (sqrt(errors[0]*errors[0]+errors[2]*errors[2]+errors[4]*errors[4])
        > l_infty_l_2)
      {
        l_infty_l_2 = sqrt(errors[0]*errors[0]+errors[2]*errors[2]
          +errors[4]*errors[4]);
        l_infty_l_2_time =  TDatabase::TimeDB->CURRENTTIME;
      }
      OutPut( l_infty_l_2_time <<  " l_infty(L2(u)) " << l_infty_l_2 << endl);

      // error in L^2(0,t,L^2)
      l_2_l_2u += (errors[0]*errors[0] + errors[2]*errors[2]
        +errors[4]*errors[4]
        +olderror_l_2_l_2u)*
        TDatabase::TimeDB->TIMESTEPLENGTH/2.0;

      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "L2(0,t,L2)(u) : " <<  sqrt(l_2_l_2u) << endl);

      olderror_l_2_l_2u = errors[0]*errors[0] + errors[2]*errors[2]+errors[4]*errors[4];

      //error in L^2(0,t,H^1)
      l_2_h_1u += (errors[1]*errors[1] + errors[3]*errors[3]+ errors[5]*errors[5]
        +olderror_l_2_h_1u)*
        TDatabase::TimeDB->TIMESTEPLENGTH/2.0;

      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "L2(0,t,H1-semi)(u) : " << sqrt(l_2_h_1u) << endl);
#ifdef  __CHANNEL_CAROLINA__
      OutPut(TDatabase::TimeDB->CURRENTTIME <<" mean value of energy dissipation rate "
        << l_2_h_1u/TDatabase::TimeDB->CURRENTTIME << endl);
#endif
      olderror_l_2_h_1u = errors[1]*errors[1] + errors[3]*errors[3]+ errors[5]*errors[5];

      U1Array[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors);

      U2Array[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors+2);

      U3Array[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors+4);
      OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
      OutPut( "kinetic energy " << (errors[0]*errors[0]+errors[2]*errors[2]
        +errors[4]*errors[4])/2);
      OutPut(endl);

      PArray[mg_level-1]->GetErrors(ExactP, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, PSpaces+mg_level-1, errors);
      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "L2(p): " << errors[0]);
      OutPut( "   H1-semi(p): " << errors[1]<<endl);

      //  if  ((TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE != 4)&&
      //  (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE != 5)
      //  &&(TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE != 100))
      // {
      //  // subgrid dissipation
      //  UArray[mg_level-1]->GetDeformationTensorErrors
      //   (ExactU1, ExactU2, ExactU3,
      //    4, TimeNSAllDerivatives,
      //   1, SubGridDissipation,
      //   NULL, aux, 1, USpaces+mg_level-1, errors);
      //
      //  OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      //  OutPut( "subgrid dissipation : " << errors[0] << endl);
      //}
      //else
      //  OutPut( "subgrid dissipation not implemented " << endl);

#ifdef  __CHANNEL_OBSTACLE__
      // errors
      U1Array[mg_level-1]->GetErrors(ExactU1, 4, TimeNSAllDerivatives,
        2, L2H1InteriorErrors,
        NULL, aux, 1, USpaces+mg_level-1, errors);

      U2Array[mg_level-1]->GetErrors(ExactU2, 4, TimeNSAllDerivatives,
        2, L2H1InteriorErrors,
        NULL, aux, 1, USpaces+mg_level-1, errors+2);

      U3Array[mg_level-1]->GetErrors(ExactU3, 4, TimeNSAllDerivatives,
        2, L2H1InteriorErrors,
        NULL, aux, 1, USpaces+mg_level-1, errors+4);
      OutPut( TDatabase::TimeDB->CURRENTTIME << " xmin " << TDatabase::ParamDB->P6 << " xmax " <<TDatabase::ParamDB->P7);
      OutPut( " L2-interior(u): " << sqrt(errors[0]*errors[0]+errors[2]*errors[2]
        +errors[4]*errors[4]));
      OutPut( "   H1-semi-interior(u):  " << sqrt(errors[1]*errors[1]+errors[3]*errors[3]
        +errors[5]*errors[5])<<endl);
      U1Array[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1InteriorErrors, NULL, aux, 1, USpaces+mg_level-1, errors);

      U2Array[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1InteriorErrors, NULL, aux, 1, USpaces+mg_level-1, errors+2);

      U3Array[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1InteriorErrors, NULL, aux, 1, USpaces+mg_level-1, errors+4);
      OutPut( TDatabase::TimeDB->CURRENTTIME << " xmin " << TDatabase::ParamDB->P6 << " xmax " <<TDatabase::ParamDB->P7);
      OutPut( " kinetic energy-interior " << (errors[0]*errors[0]+errors[2]*errors[2]
        +errors[4]*errors[4])/2 << endl);
#endif

      //OutPut( "L2(p): " << errors[0] << endl);
      //OutPut( "H1(p): " << errors[1] << endl);

      /* l2p[mg_level-1] = errors[0];
         h1p[mg_level-1] = errors[1];

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
      OutPut( "Du : " << errors[0] << " " << sqrt(l_2_l_2Du) << endl);
      */

      UArray[mg_level-1]->GetDeformationTensorErrors
        (ExactNull, ExactNull, ExactNull,
        3, TimeNSAllDerivatives,
        2, DivergenceError,
        NULL, aux, 1, USpaces+mg_level-1, errors);

      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "divergence error (L1/L2) : " << errors[0]*errors[0]  <<
        " " << errors[1] << endl);
      delete aux;
    }                                             // endif MEASURE_ERRORS

    if (TDatabase::ParamDB->COMPUTE_VORTICITY_DIVERGENCE)
    {
      ComputeVorticityDivergence(USpaces[mg_level-1], U1Array[mg_level-1],
        U2Array[mg_level-1], U3Array[mg_level-1],
        vorticity_space, Vorticity->GetComponent(0),
        Vorticity->GetComponent(1),
        Vorticity->GetComponent(2),
        div);
    }
#ifdef __MIXINGLAYERSLIP3D__
    ComputeVorticityThickness3D(Vorticity->GetComponent(2),errors);
    OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
    OutPut( "vorticity thickness: " << errors[0] << " " <<
      errors[0]/vort_zero << endl);
    ComputeMomentumThickness3D(U1Array[mg_level-1],errors);
    OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
    OutPut( "momentum thickness: " << errors[0] << " " <<
      errors[0]/moment_zero << endl);
    comp_vort = 1;

    aux =  new TAuxParam3D(TimeNSN_FESpacesVelo, TimeNSN_FctVelo,
      TimeNSN_ParamFctVelo,
      TimeNSN_FEValuesVelo,
      fesp, fefct,
      TimeNSFctVelo,
      TimeNSFEFctIndexVelo, TimeNSFEMultiIndexVelo,
      TimeNSN_ParamsVelo, TimeNSBeginParamVelo);
    // enstrophy
    Vorticity->GetComponent(0)->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
      2, L2H1Errors,
      NULL, aux, 1, &vorticity_space, errors);

    Vorticity->GetComponent(1)->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
      2, L2H1Errors,
      NULL, aux, 1, &vorticity_space, errors+2);

    Vorticity->GetComponent(2)->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
      2, L2H1Errors,
      NULL, aux, 1, &vorticity_space, errors+4);
    OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
    OutPut( "enstrophy: " << (errors[0]*errors[0]+errors[2]*errors[2]
      +errors[4]*errors[4])/2.0 << endl);
    delete aux;
#endif

#ifdef __CHANNEL_TAU180__
    if ((TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION) ||
      (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION_EXPL))
    {
      ComputeVMSProjection(Matrices_G11[mg_level-1], Matrices_G22[mg_level-1],
        Matrices_G33[mg_level-1], MatricesL[mg_level-1],
        U1Array[mg_level-1], U2Array[mg_level-1],
        U3Array[mg_level-1], vms_projection_fe);
    }

    if (TDatabase::TimeDB->CURRENTTIME>=TDatabase::TimeDB->T0)
    {
      if (TDatabase::ParamDB->INTERNAL_MEAN_COMPUTATION == 0)
      {
        TDatabase::ParamDB->INTERNAL_MEAN_COMPUTATION = 1;
        TDatabase::TimeDB->T0 = TDatabase::TimeDB->CURRENTTIME;
        // reset mean_velocity
        memset(mean_velocity,0,N_z_layers*SizeOfDouble);
      }
    }

    ComputeMeanVelocity(coll, USpaces[mg_level-1], projection_space,
      N_z_layers, coord_z_layers,  rms_velocity,
      rms_velocity2, rms_velocity3, rms_velocity1_type1,
      rms_velocity2_type1, rms_velocity3_type1,
      mean_velocity, mean_velocity_u2, mean_velocity_u3,
      dmean_velocity,
      R_xx, R_xy, R_xz, R_yy, R_yz, R_zz,
      A_M_xx, A_M_xy, A_M_xz, A_M_yy, A_M_yz, A_M_zz,
      U1Array[mg_level-1], U2Array[mg_level-1], U3Array[mg_level-1],
      vms_projection_fe,
      x_dof, y_dof, z_dof, u1x, projection_u1x);

    // increase the damping for the nonlinear iteration
#endif

    if (TDatabase::ParamDB->CONVOLUTE_SOLUTION)
    {
      ConvolveSolution(MG, USpaces,
        U1Array, U2Array, U3Array,
        DiscreteFormRHSAuxProblemU,
        sqmatrixGL00AuxProblem,
        rhsGL00AuxProblem,
        u_uConv,
        mg_level, N_U);

      aux =  new TAuxParam3D(TimeNSN_FESpacesVelo, TimeNSN_FctVelo,
        TimeNSN_ParamFctVelo,
        TimeNSN_FEValuesVelo,
        fesp, fefct,
        TimeNSFctVelo,
        TimeNSFEFctIndexVelo, TimeNSFEMultiIndexVelo,
        TimeNSN_ParamsVelo, TimeNSBeginParamVelo);

      // error in first component
      u1ConvArray[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors);

      // error in second component
      u2ConvArray[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors+2);

      // error in third component
      u3ConvArray[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors+4);

      OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
      OutPut( "conv H1-semi(u): " << sqrt(errors[1]*errors[1]+errors[3]*errors[3]
        +errors[5]*errors[5]) << endl);
      OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
      OutPut( "conv kinetic energy " << (errors[0]*errors[0]+errors[2]*errors[2]
        +errors[4]*errors[4])/2<< endl);

#ifdef  __CHANNEL_OBSTACLE__
      u1ConvArray[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1InteriorErrors,
        NULL, aux, 1, uConvSpaces+mg_level-1, errors);

      u2ConvArray[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1InteriorErrors,
        NULL, aux, 1, uConvSpaces+mg_level-1, errors+2);

      u3ConvArray[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1InteriorErrors,
        NULL, aux, 1, uConvSpaces+mg_level-1, errors+4);
      OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
      OutPut( "conv H1-semi-interior(u): " << sqrt(errors[1]*errors[1]+errors[3]*errors[3]
        +errors[5]*errors[5]) << endl);
      OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
      OutPut( "conv kinetic energy interior " << (errors[0]*errors[0]+errors[2]*errors[2]
        +errors[4]*errors[4])/2<< endl);
      cout << errors[0] << " " << errors[2] << " " << errors[4] << endl;
#endif

#ifdef __MIXINGLAYERSLIP3D__
      // compute vorticity
      ComputeVorticityDivergence(velocity_space, u1ConvArray[mg_level-1],
        u2ConvArray[mg_level-1], u3ConvArray[mg_level-1],
        vorticity_space,Vorticity->GetComponent(3),
        Vorticity->GetComponent(4),
        Vorticity->GetComponent(5),div);

      ComputeVorticityThickness3D(Vorticity->GetComponent(5),errors);
      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "vorticity thickness (uconv): " << errors[0] << " " <<
        errors[0]/vort_zero_conv << endl);
      ComputeMomentumThickness3D(u1ConvArray[mg_level-1],errors);
      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "momentum thickness (uconv): " << errors[0] << " " <<
        errors[0]/moment_zero_conv << endl);
      // enstrophy
      Vorticity->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1Errors,
        NULL, aux, 1,  &vorticity_space, errors);

      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "enstrophy (uconv) " << (errors[0]*errors[0]+errors[2]*errors[2]
        +errors[4]*errors[4])/2.0  << endl);

      // smooth filtered vorticity for output
      // copy vorticity to sol to use mg structure
      if (mixing_layer_galerkin)
      {
        OutPut("smooth solution " << endl);
        memcpy(sol_vort_tmp,sol,N_U*SizeOfDouble);
        memcpy(sol,Vorticity->GetComponent(5)->GetValues(),N_U*SizeOfDouble);
        // restrict
        auxConv = new double[N_U];
        for(ii = mg_level-1 ; ii > mg_level-1-smoothing_depth;ii--)
        {
          RestrictFunction(USpaces[ii-1], USpaces[ii],
            U1Array[ii-1]->GetValues(),
            U1Array[ii]->GetValues(),
            auxConv);
        }
        // prolongate restricted vorticity
        for(ii = mg_level-1-smoothing_depth ; ii < mg_level-1;ii++)
        {
          Prolongate(USpaces[ii], USpaces[ii+1],
            U1Array[ii]->GetValues(),
            U1Array[ii+1]->GetValues(),
            auxConv);

        }
        // copy result to vorticity vector
        memcpy(Vorticity->GetComponent(5)->GetValues(),sol,N_U*SizeOfDouble);
        memcpy(sol,sol_vort_tmp,N_U*SizeOfDouble);
        // restrict solution to all grids
        MG->RestrictToAllGrids();
        delete auxConv;
      }
#endif
      delete aux;
    }                                             // end CONVOLUTE_SOLUTION

#ifdef __BENCH__
    // compute characteristic values (deltaP, Cd, Cl)
    GetCdCl(U1Array[mg_level-1], U2Array[mg_level-1], U3Array[mg_level-1],
      PArray[mg_level-1],  U1old, U2old, Cd, Cl);

    PArray[mg_level-1]->FindGradient(0.45, 0.2, 0.205, dP1);
    PArray[mg_level-1]->FindGradient(0.55, 0.2, 0.205, dP2);

    OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
    OutPut( "C_drag = " << setprecision(16) <<Cd );
    OutPut( " C_lift = " << setprecision(16) << Cl);
    OutPut( " deltaP = " << setprecision(16) << dP1[0] - dP2[0] << endl);
    OutPut( setprecision(7) << endl);
    memcpy(former_sol,sol,2*N_U*SizeOfDouble);
    if ( TDatabase::TimeDB->CURRENTTIME >= TDatabase::TimeDB->T1-1e-5)
    {
	ComputeFrictionVelocities(coll, U1Array[mg_level-1], U2Array[mg_level-1],
				  velo_friction, fric_count);
	PressureAtCylinder(coll, PArray[mg_level-1], press_cyl, N_press_cyl, press_count);
	CenterlineVelocities(coll, U1Array[mg_level-1], U2Array[mg_level-1],
			     center_velo, N_center_velo, center_velo_count);
	VelocityAtCylinder(coll, U1Array[mg_level-1], 
			   cyl_velo, N_cyl_velo, cyl_velo_count);
    }
#endif

#ifdef  __CHANNEL_OBSTACLE__
    // compute characteristic values (deltaP, Cd, Cl)
    GetCdCl(U1Array[mg_level-1], U2Array[mg_level-1], U3Array[mg_level-1],
      PArray[mg_level-1], U1old, U2old, Cd, Cl);

    PArray[mg_level-1]->FindGradient(2.45, 0.2, 0.205, dP1);
    PArray[mg_level-1]->FindGradient(2.55, 0.2, 0.205, dP2);

    OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
    OutPut( "C_drag = " << setprecision(16) <<Cd );
    OutPut( " C_lift = " << setprecision(16) << Cl);
    OutPut( " deltaP = " << setprecision(16) << dP1[0] - dP2[0] << endl);
    if (TDatabase::ParamDB->CONVOLUTE_SOLUTION)
    {
      // !!! time derivative in computation of drag and lift is neglected !!!
      //cout << dP1[0] << " " << dP2[0] << endl;
      GetCdCl(u1ConvArray[mg_level-1], u2ConvArray[mg_level-1],
        u3ConvArray[mg_level-1], pConvArray[mg_level-1],
        u1ConvArray[mg_level-1], u2ConvArray[mg_level-1],
        Cd, Cl);
      pConvArray[mg_level-1]->FindGradient(2.45, 0.2, 0.205, dP1);
      pConvArray[mg_level-1]->FindGradient(2.55, 0.2, 0.205, dP2);

      OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
      OutPut( "conv velo pres: C_drag (wo. time der.) = " << setprecision(16) <<Cd );
      OutPut( " C_lift (wo. time der.) = " << setprecision(16) << Cl);
      OutPut( " deltaP = " << setprecision(16) << dP1[0] - dP2[0] << endl);

      GetCdCl(u1ConvArray[mg_level-1], u2ConvArray[mg_level-1],
        u3ConvArray[mg_level-1], PArray[mg_level-1],
        u1ConvArray[mg_level-1], u2ConvArray[mg_level-1],
        Cd, Cl);
      OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
      OutPut( "conv velo: C_drag (wo. time der.) = " << setprecision(16) <<Cd );
      OutPut( " C_lift (wo. time der.) = " << setprecision(16) << Cl);
      OutPut( " deltaP = " << setprecision(16) << dP1[0] - dP2[0] << endl);
    }
    OutPut( setprecision(7) << endl);
    memcpy(former_sol,sol,2*N_U*SizeOfDouble);
#endif

#ifdef __CHANNELSTEP__
    OutPut("reatttime "<< TDatabase::TimeDB->CURRENTTIME << " " );
    GetReattachmentLine(U1Array[mg_level-1], reatt_pt);
    //OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
    //OutPut( "reattachment: " << reatt_pt<< endl);
#endif

#ifdef __CALOTTE__
    for (j=0;j<N_neum_bdry;j++)
    {
      OutPut("neumann on bdry " << sol[neum_bdry[j]] << " " <<  sol[neum_bdry[j]+N_U]
        << " "  << sol[neum_bdry[j]+2*N_U] << endl);
    }
#endif

    if ((TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION) ||
      (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION_EXPL))
    {
      ComputeVMSProjection(Matrices_G11[mg_level-1], Matrices_G22[mg_level-1],
        Matrices_G33[mg_level-1], MatricesL[mg_level-1],
        U1Array[mg_level-1], U2Array[mg_level-1],
        U3Array[mg_level-1], vms_projection_fe);
      // compute the size of the small scales in VMS methods
      ComputeSizeOfSmallScales(Matrices_G11[mg_level-1], Matrices_G22[mg_level-1],
        Matrices_G33[mg_level-1], MatricesL[mg_level-1],
        U1Array[mg_level-1], U2Array[mg_level-1],
        U3Array[mg_level-1], vms_projection_fe, size_small_scales);
      // time average for mean, and largest_size
      MeanAndLargestSize( projection_space, size_small_scales, &mean, &largest_size);
      mean_time_average = m * (mean_time_average/(m+1))  + mean/(m+1);
      max_time_average =  m * (max_time_average/(m+1))  + largest_size/(m+1);

      OutPut(TDatabase::TimeDB->CURRENTTIME << " small scales: max_time_average  " <<
        max_time_average << " mean_time_average " << mean_time_average << endl);

      // adaptive large scale space
      if ((TDatabase::ParamDB->VMS_LARGE_VELOCITY_SPACE==17)&& (m%TDatabase::ParamDB->VMS_ADAPT_STEPS==0))
      {
        delete sqstructureL;
        delete structure_tilde_G;
        delete structure_G;
        delete sqmatrixL;
        delete matrix_tilde_G11;
        delete matrix_tilde_G22;
        delete matrix_tilde_G33;
        delete matrix_G11;
        delete matrix_G22;
        delete matrix_G33;
        delete vms_projection;
        delete vms_projection_fe;
        delete vms_proj_11;
        delete vms_proj_12;
        delete vms_proj_13;
        delete vms_proj_22;
        delete vms_proj_23;
        delete vms_proj_33;

        // set for the large scale space the fe orders
        // for the individual mesh cells
        AdaptProjectionSpace(projection_space, size_small_scales, fes,
          mean, mean_time_average, largest_size,
          max_time_average, label_space);
        // define new large scale space

        // delete everything which is not longer needed
        // new large scale space
        fes1 = new FE3D[coll->GetN_Cells()];
        memcpy(fes1,fes,coll->GetN_Cells()*sizeof(FE3D));
        // this command deletes fes
        delete projection_space;
        OutPut("done delete "<<endl);
        fes = fes1;
        OutPut("start space "<<endl);
        projection_space = new TFESpace3D(coll, NameString, PString, BoundCondition,
          fes);
        OutPut("done space "<<endl);

        OutPut("projection space adapted, dof " <<
          projection_space->GetN_DegreesOfFreedom()<<endl);
        ProjectionSpaces[mg_level-1] = projection_space;
        sqstructureL = new TSquareStructure3D(projection_space);
        sqstructureL->Sort();
        structure_tilde_G = new TStructure3D(velocity_space, projection_space);
        structure_G = new TStructure3D(projection_space, velocity_space);
        sqmatrixL = new TSquareMatrix3D(sqstructureL);
        MatricesL[mg_level-1] = sqmatrixL;
        LumpMassMatrixToDiag(MatricesL[mg_level-1]);
        matrix_tilde_G11 = new TMatrix3D(structure_tilde_G);
        Matrices_tilde_G11[mg_level-1] = matrix_tilde_G11;
        matrix_tilde_G22 = new TMatrix3D(structure_tilde_G);
        Matrices_tilde_G22[mg_level-1] = matrix_tilde_G22;
        matrix_tilde_G33 = new TMatrix3D(structure_tilde_G);
        Matrices_tilde_G33[mg_level-1] = matrix_tilde_G33;
        matrix_G11 = new TMatrix3D(structure_G);
        Matrices_G11[mg_level-1] = matrix_G11;
        matrix_G22 = new TMatrix3D(structure_G);
        Matrices_G22[mg_level-1] = matrix_G22;
        matrix_G33 = new TMatrix3D(structure_G);
        Matrices_G33[mg_level-1] = matrix_G33;
        N_L = projection_space->GetN_DegreesOfFreedom();
        OutPut(TDatabase::TimeDB->CURRENTTIME << " dof projection : " << setw(10) << N_L << endl);

        vms_projection = new double[6*N_L];
        memset(vms_projection,0,6*N_L*SizeOfDouble);
        vms_projection_fe = new TFEVectFunct3D(projection_space, PString, PString,
          vms_projection, N_L, 6);
        vms_proj_11 = vms_projection_fe->GetComponent(0);
        vms_proj_12 = vms_projection_fe->GetComponent(1);
        vms_proj_13 = vms_projection_fe->GetComponent(2);
        vms_proj_22 = vms_projection_fe->GetComponent(3);
        vms_proj_23 = vms_projection_fe->GetComponent(4);
        vms_proj_33 = vms_projection_fe->GetComponent(5);
      }
    }

    if ((TDatabase::ParamDB->WRITE_GRAPE)||(TDatabase::ParamDB->WRITE_GMV)
      || (TDatabase::ParamDB->WRITE_AMIRA)||(TDatabase::ParamDB->WRITE_VTK)
      || (TDatabase::ParamDB->SAVE_DATA))
    {
      if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      {
        if (TDatabase::ParamDB->COMPUTE_VORTICITY_DIVERGENCE)
        {
          if (!comp_vort)
            ComputeVorticityDivergence(USpaces[mg_level-1], U1Array[mg_level-1],
              U2Array[mg_level-1], U3Array[mg_level-1],
              vorticity_space, Vorticity->GetComponent(0),
              Vorticity->GetComponent(1),
              Vorticity->GetComponent(2),
              div);
        }
        if (TDatabase::ParamDB->WRITE_GMV)
        {
          os.seekp(std::ios::beg);
          os << GMVBaseName << m << ".gmv" << ends;
          Output->WriteGMV(os.str().c_str());
        }
        if (TDatabase::ParamDB->WRITE_VTK)
        {
          os.seekp(std::ios::beg);
          os << VTKBaseName << m << ".vtk" << ends;
          Output->WriteVtk(os.str().c_str());
        }
        if (TDatabase::ParamDB->WRITE_AMIRA)
        {
          os.seekp(std::ios::beg);
          os << GrapeBaseName << m << ".am" << ends;
          Output->WriteAmira(os.str().c_str());
        }
        if (TDatabase::ParamDB->WRITE_GRAPE)
        {
          os.seekp(std::ios::beg);
          os << GrapeBaseName << m << ".dat" << ends;
          Output->WriteGrape(os.str().c_str());
        }
        if (TDatabase::ParamDB->SAVE_DATA)
        {
          save_sol[0] = sol;
          save_N_Unknowns[0] = N_Unknowns;
          SaveData(SaveDataFileName,1,save_sol,save_N_Unknowns);
        }
        N_GRAPE_images++;
      }
    }

    comp_vort =0;
  }                                               // while

  //======================================================================
  // end of time cycle
  //======================================================================

  if ((TDatabase::ParamDB->WRITE_GRAPE)||(TDatabase::ParamDB->WRITE_GMV)
    ||(TDatabase::ParamDB->WRITE_VTK))
  {
    if (TDatabase::ParamDB->COMPUTE_VORTICITY_DIVERGENCE)
    {
      ComputeVorticityDivergence(USpaces[mg_level-1], U1Array[mg_level-1],
        U2Array[mg_level-1],U3Array[mg_level-1],
        vorticity_space, Vorticity->GetComponent(0),
        Vorticity->GetComponent(1),
        Vorticity->GetComponent(2),
        div);
    }
    if (TDatabase::ParamDB->WRITE_GMV)
    {
      os.seekp(std::ios::beg);
      os << GMVBaseName << m << ".gmv" << ends;
      Output->WriteGMV(os.str().c_str());
    }
    if (TDatabase::ParamDB->WRITE_VTK)
    {
      os.seekp(std::ios::beg);
      os << VTKBaseName << m << ".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
    }
    if (TDatabase::ParamDB->WRITE_AMIRA)
    {
      os.seekp(std::ios::beg);
      os << GrapeBaseName << m << ".am" << ends;
      Output->WriteAmira(os.str().c_str());
    }
    if (TDatabase::ParamDB->WRITE_GRAPE)
    {
      os.seekp(std::ios::beg);
      os << GrapeBaseName << m << ".dat" << ends;
      Output->WriteGrape(os.str().c_str());
    }
  }

  t4 =  GetTime();
  total_time += t4 - t3;
  OutPut("total running time: " << total_time
    << " assembling " << ass_time << " " << 100*ass_time/total_time << " "
    << " solver " << solver_time << " " << 100*solver_time/total_time << endl );

  CloseFiles();
  return 0;
}


// start with more than 4300 lines
