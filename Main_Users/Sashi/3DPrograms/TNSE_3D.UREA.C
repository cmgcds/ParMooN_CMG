// =======================================================================
//
// Purpose:     main program
//
// Author:      Volker John, Carina Suciu, Ellen Schmeyer
//
// =======================================================================

#include <Assemble3D.h>
#include <AuxParam3D.h>
#include <BoundFace.h>
#include <Bulk_3d4d.h>
#include <Collection.h>
#include <Convolution.h>
#include <Database.h>
#include <DirectSolver.h>
#include <DiscreteForm3D.h>
#include <Domain.h>
#include <FEDatabase3D.h>
#include <FEM_TVD_FCT.h>
#include <FESpace3D.h>
#include <FgmresIte.h>
#include <FixedPointIte.h>
#include <ItMethod.h>
#include <LinAlg.h>
#include <math.h>
#include <MGLevel3D.h>
#include <MultiGrid3D.h>
#include <MultiGridIte.h>
#include <MultiGridScaIte.h>
#include <NSE3D_ParamRout.h>
#include <NSE_MultiGrid.h>
#include <NSE_MGLevel1.h>
#include <NSE_MGLevel2.h>
#include <NSE_MGLevel3.h>
#include <NSE_MGLevel4.h>
#include <Output3D.h>
#include <QuadAffin.h>
#include <RationalLES.h>
#include <RFB.h>
#include <Solver.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <TCD3D.h>
#include <TNSE3D_ParamRout.h>
#include <Upwind.h>
#include <Upwind3D.h>
#include <Urea_3d4d.h>
#include <VMS.h>
#include <fftw3.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <malloc.h>

#include <MainUtilities.h>
#include <TimeUtilities.h>

double bound = 0;

// ======================================================================
// utilities for main program
// ======================================================================

#define AMG 0
#define GMG 1
#define DIRECT 2

// =======================================================================
// include current example
// =======================================================================

#include "../Examples/TNSE_3D/Urea.h"

// =======================================================================
// start of main programm
// =======================================================================

int main(int argc, char* argv[])
{
  //======================================================================
  // begin of the declaration of the variables
  //======================================================================

  // integer variables
  int i, j, k, l, m, n, ii, ll, i_conc, N3, N4;
  int ansatz_order_reaction;
  int CurrentDiscType, comp_vort;
  int LEVELS, last_sq;
  int Max_It, Max_It_c, methods, mg_level, mg_type, mid_sq, low_scalar, CONC_MAXIT;
  int N, N_, N_Active, N_Active_c, N_Active_temp, N_Active_conc;
  int n_aux, N_Cells, N_Entries, N_FESpaces, N_L, N_LinIter, N_LinIterCurr, N_LinIterCurrIte;
  int Nodes, nonlinite, nonlin_ite, N_P, N_RectMatrices, N_Rhs, N_SquareMatrices, N_SubSteps;
  int N_U, N_UConv, N_Unknowns, N_vort, N_Vort, N_z_layers;
  int N_Unknowns_c, N_Unknowns_temp, N_Unknowns_conc, N_Unknowns_Integral_Space;
  int N_V, Nx, Ny, Nz, Na;
  int order, only_first_time, pressure_space_code;
  int ret, smoothing_depth, solver_type_reaction;
  int time_discs, velocity_space_code, zerostart;
        //N3=(Nx+1)*(Ny+1)*(Nz+1);
        //N4=N3*(Na+1);

  // initialised integer variables
  int first_psd_fem = 1;
  int mixing_layer_galerkin = 0;
  int N_neum_to_diri_c = 0, N_neum_to_diri_temp = 0, N_neum_to_diri_conc = 0;
  int N_neum_to_diri_psd = 0;
  int N_Paramters = 1, very_first_time = 0, substance = 1;
  int average_step[1];
  int steady = 0;
  int INTERNAL_STEADY_STATE_MATRICES_OR_RHS;
  // integer pointers
  int *col_ptr, *correspond_3dgrid, *N_Parray, *N_Uarray, *RowPtr;
  int *N_Array, *N_Array_c, *N_Array_temp, *N_Array_conc, *row_ptr;
  int *neum_to_diri, N_neum_to_diri;
  int *neum_to_diri_temp;
  int *neum_to_diri_conc, *neum_to_diri_c;
  int *neum_to_diri_psd;
  int *index_test_ansatz;

  int **downwind;

  // double variables
  double end_time, gamma, gamma_c, h;
  double hmin, hmax, impuls_residual;
  double limit, limit_c;
  double moment_zero, moment_zero_conv;
  double oldresidual, oldtau, omega, p2;
  double real_time, reatt_pt, res, residual, residual_conc;
  double solver_time, solver_time_curr;
  double t, t1, t11, t2, t22, t3, t4, tau, tau1, tau2, theta1, theta2, theta3, theta4, total_time;
  double vort_zero, vort_zero_conv;
  double x, y;
  double CONC_TOL;
  // initialised double variables
  double l_infty_l_2 = 0, l_infty_l_2_time = -4711.0, l_2_h_1u = 0, l_2_l_2u = 0;
  double max_conc = -1, nonlin_resid = 1e8, olderror_l_2_l_2u = 0, olderror_l_2_h_1u = 0;
  double average_median[1];
  //double *average_q3;
  double max_c[6];

  // fixing of the interval limits for the coordinates
  double x_min, x_max;
  double y_min, y_max;
  double z_min, z_max;
  double a_min = 0, a_max = 1;
  // coordinate for the 3d cut -> visulization in paraview
  double cut_coord = 0;
  // declaration and initialisation of the coordinate vectors
  double *x_coord, *y_coord, *z_coord, *a_coord;
  // declaration and initialisation of the vector with the differences between the z layers
  double *a_layers_coord;

  // double pointers
  double *A_M_xx, *A_M_yy, *A_M_zz, *A_M_xy, *A_M_xz, *A_M_yz, *auxConv;
  double *B, *B_c;
  double *coord_z_layers, *concent_C_array,*Temp_array, *current_B, *current_sol, *current_sol_c;
  double *defect, *defect_c, *div, *dmean_velocity;
  double *frac_step_sol, *sol_psd_old,*sol_psd_help, *sol_psd, *h1p, *rhs_psd,*rhs_psd_old ;
  double *integral_val, *itmethod_rhs, *itmethod_sol;
  double *itmethod_rhs_c, *itmethod_sol_c, *itmethod_sol_temp, *itmethod_sol_conc;
  double *l_inf, *l2p, *LESModelRhs;
  double *lump_mass_PSD, *lump_mass_temp, *lump_mass_conc, *lump_mass_c;
  double *matrix_D_Entries_PSD, *matrix_D_Entries_temp, *matrix_D_Entries_conc, *matrix_D_Entries_c;
  double *mean_velocity, *mean_velocity_u2, *mean_velocity_u3;
  double *neum_to_diri_param_temp, *neum_to_diri_param_conc, *neum_to_diri_param_c;
  double *neum_to_diri_temp_x, *neum_to_diri_temp_y, *neum_to_diri_temp_z;
  double *neum_to_diri_conc_x, *neum_to_diri_conc_y, *neum_to_diri_conc_z;
  double *neum_to_diri_x, *neum_to_diri_y, *neum_to_diri_z;
  double *neum_to_diri_psd_x, *neum_to_diri_psd_y, *neum_to_diri_psd_z, *neum_to_diri_psd_a;
  double *newton;
  double *oldrhs, *oldrhs_c, *oldrhs_temp, *oldrhs_conc;
  double *oldrhs_fem_fct0_temp, *oldrhs_fem_fct0_conc;
  double *oldrhs_fem_fct1_temp, *oldrhs_fem_fct1_conc;
  double *oldrhs_fem_fct0_c, *oldrhs_fem_fct1_c;
  double *oldsol, *old_sol, *oldsol_c, *oldsol_temp, *oldsol_conc;
  double *projection_u1x, *psi;
  double *R_xx, *R_yy, *R_zz, *R_xy, *R_xz, *R_yz, *rhsGL00AuxProblem, *rhs_vms_expl;
  double *rhs, *rhs_c, *rhs_c_complete, *rhs_c_complete_A, *rhs_c_complete_B;
  double *rms_velocity, *rms_velocity2, *rms_velocity3, *rms_velocity1_type1, *rms_velocity2_type1, *rms_velocity3_type1;
  double *sol, *solGL00AuxProblem, *sol_timestep_m1, *sol_vort_tmp, *startsol;
  double *sol_c, *sol_temp, *sol_conc, *startsol_c;
  double *tilde_u_c, *tilde_u_temp, *tilde_u_conc;
  double *u1x, *u_conv, *u_uConv;
  double *velo1, *velo2, *velo3, *vms_projection, *vorticity;
  double *x_dof, *y_dof, *z_dof;
  double *size_small_scales,*label_space;

  double **RhsArray, **RhsArray_c, **RhsArray_temp, **RhsArray_conc;

  // double pointer on an array
  double *RHSs[6];

  // double arrays
  double alpha[2], errors[9], Parameters[4];
  double *save_sol[5];
  int save_N_Unknowns[5];

  // variables of MooNMD classees
  BoundCondFunct3D *BoundaryConditions[6], *BoundaryConditionsAuxProblem[6];
  BoundCondFunct3D *BoundaryConditionsPressureSeparation[1];
  BoundCondFunct3D *BoundaryConditions_Scalar[3]; // [3] oder [6] ? Michael

  BoundValueFunct3D *BoundValues[6], *BoundValuesAuxProblem[6];
  BoundValueFunct3D *BoundaryValuesPressureSeparation[1];
  BoundValueFunct3D *BoundValues_Scalar[3];       // [3] oder [6] ? Michael

  CoeffFct3D *Coefficients[4], *Coeff_c;

  DefectProc *Defect, *DefectScalar;

  MatVecProc *MatVect, *MatVectScalar;

  TAuxParam3D *aux, *auxn;

  TBaseCell *cell;

  TCollection *coll, *mortarcoll = NULL;

  TDatabase *Database = new TDatabase();

  TDiscreteForm3D *DiscreteForm;
  TDiscreteForm3D *DiscreteFormC;
  TDiscreteForm3D *DiscreteFormClassicalLES;
  TDiscreteForm3D *DiscreteFormGalerkin;
  TDiscreteForm3D *DiscreteFormGL00AuxProblem;
  TDiscreteForm3D *DiscreteFormGL00AuxProblemRHS;
  TDiscreteForm3D *DiscreteFormGL00Convolution;
  TDiscreteForm3D *DiscreteFormJ;
  TDiscreteForm3D *DiscreteFormMatricesA_Urea;
  TDiscreteForm3D *DiscreteFormMatricesA_Urea_conc;
  TDiscreteForm3D *DiscreteFormMatricesA_Galerkin_Urea;
  TDiscreteForm3D *DiscreteFormMatricesA_Galerkin_Urea_conc;
  TDiscreteForm3D *DiscreteFormMatricesA_SUPG_Urea;
  TDiscreteForm3D *DiscreteFormMatricesA_SUPG_Urea_conc;
  TDiscreteForm3D *DiscreteFormMatrixAuxProblemU;
  TDiscreteForm3D *DiscreteFormMatrixGL00AuxProblem;
  TDiscreteForm3D *DiscreteFormMatrixMUrea;
  TDiscreteForm3D *DiscreteFormNLClassicalLES;
  TDiscreteForm3D *DiscreteFormNLGalerkin;
  TDiscreteForm3D *DiscreteFormNLGL00Convolution;
  TDiscreteForm3D *DiscreteFormNLGL00AuxProblem;
  TDiscreteForm3D *DiscreteFormNLSmagorinsky;
  TDiscreteForm3D *DiscreteFormNLUpwind;
  TDiscreteForm3D *DiscreteFormNLUpwindNC;
  TDiscreteForm3D *DiscreteFormNLVMS_Projection;
  TDiscreteForm3D *DiscreteFormNLVMS_ProjectionExpl;
  TDiscreteForm3D *DiscreteFormNLVMS_RFBExplRhs;
  TDiscreteForm3D *DiscreteFormRHS;
  TDiscreteForm3D *DiscreteFormRhs_Urea;
  TDiscreteForm3D *DiscreteFormRhs_Urea_conc;
  TDiscreteForm3D *DiscreteFormRhs_Galerkin_Urea;
  TDiscreteForm3D *DiscreteFormRhs_Galerkin_Urea_conc;
  TDiscreteForm3D *DiscreteFormRHSAuxProblemU;
  TDiscreteForm3D *DiscreteFormRHSClassicalLES;
  TDiscreteForm3D *DiscreteFormRHSLESModel;
  TDiscreteForm3D *DiscreteFormRHSNewton;
  TDiscreteForm3D *DiscreteFormRHSNewtonNL;
  TDiscreteForm3D *DiscreteFormRhs_SUPG_Urea;
  TDiscreteForm3D *DiscreteFormRhs_SUPG_Urea_conc;
  TDiscreteForm3D *DiscreteFormSmagorinsky;
  TDiscreteForm3D *DiscreteFormUpwind;
  TDiscreteForm3D *DiscreteFormUpwindNC;
  TDiscreteForm3D *DiscreteFormVMS_Projection;
  TDiscreteForm3D *DiscreteFormVMS_SUPG;
  TDiscreteForm3D *DiscreteFormNLVMS_SUPG;
  TDiscreteForm3D *DiscreteFormRHSSUPG;

  TDomain *Domain = new TDomain();

  TFEDatabase3D *FEDatabase = new TFEDatabase3D();

  TFEFunction3D *Approx, *AuxPArray, *temp, *conc, *Divergence;
  TFEFunction3D *du11Conv, *du12Conv, *du13Conv, *du22Conv, *du23Conv, *du33Conv;
  TFEFunction3D *GL00AuxProblemSol11, *GL00AuxProblemSol12;
  TFEFunction3D *GL00AuxProblemSol13, *GL00AuxProblemSol22;
  TFEFunction3D *GL00AuxProblemSol23, *GL00AuxProblemSol33;
  TFEFunction3D *integral_space_conc_fct;
  TFEFunction3D *old_p, *old_u, *p, *pConv;
  TFEFunction3D *separated_pressure_fe_funct, *separated_pressure_rhs_fe_funct;
  TFEFunction3D *soldiff_fe1,*soldiff_fe2, *StreamFct;
  TFEFunction3D *u1, *u2, *u3, *u1Conv, *u2Conv, *u3Conv, *u4Conv, *u5Conv, *u6Conv;
  TFEFunction3D *vort1, *vort2, *vort3, *Vort_x, *Vort_y, *Vort_z;
  TFEFunction3D *vms_proj_11, *vms_proj_12, *vms_proj_13, *vms_proj_22;
  TFEFunction3D *vms_proj_23, *vms_proj_33, *size_small_scales_fefct, *label_space_fefct;

  TFEFunction3D *fefct[12], *fefct_conc[10];

  TFEFunction3D **AuxFEFunctArray;
  TFEFunction3D **du11ConvArray, **du12ConvArray, **du13ConvArray;
  TFEFunction3D **du22ConvArray, **du23ConvArray, **du33ConvArray;
  TFEFunction3D **GL00AuxProblemSol11Array, **GL00AuxProblemSol12Array;
  TFEFunction3D **GL00AuxProblemSol13Array, **GL00AuxProblemSol22Array;
  TFEFunction3D **GL00AuxProblemSol23Array, **GL00AuxProblemSol33Array;
  TFEFunction3D **IntegralSpaces_conc_fct, **PArray, **pConvArray;
  TFEFunction3D **SolArray_temp, **SolArray_conc;
  TFEFunction3D **SolArray_temp_old, **SolArray_conc_old;
  TFEFunction3D **SolArray_old, **SolArray;
  TFEFunction3D **u1ConvArray, **u2ConvArray, **u3ConvArray;
  TFEFunction3D **u4ConvArray, **u5ConvArray, **u6ConvArray;
  TFEFunction3D **U1Array, **U2Array, **U3Array;

  TFEVectFunct3D *u, **UArray, *uconf, *duConv, **duConvArray;
  TFEVectFunct3D *uConv, **uConvArray, **AuxFEVectFunctArray;
  TFEVectFunct3D *Vorticity, *vms_projection_fe;
  TFEVectFunct3D *GL00AuxProblemSol, **GL00AuxProblemSolArray;

  TFESpace3D *concentration_space,  *concentration_space_c;
  TFESpace3D *convolution_space, *integral_space_conc;
  TFESpace3D *old_p_space, *old_u_space;
  TFESpace3D *pressure_separation_space, *pressure_space, *projection_space;
  TFESpace3D *streamfunction_space;
  TFESpace3D *velocity_space, *vorticity_space, *size_small_scales_fesp, *label_space_fesp;

  TFESpace3D *ferhs[6], *fesp[4], *fesp_conc[6];

  TFESpace3D **ConcentrationSpaces;
  TFESpace3D **ConcentrationSpaces_temp, **ConcentrationSpaces_conc;
  TFESpace3D **duConvSpaces, **IntegralSpaces_conc;
  TFESpace3D **pConvSpaces, **ProjectionSpaces, **PsiSpaces, **PSpaces;
  TFESpace3D **uConvSpaces, **USpaces;
  FE3D *fes, *fes1;

  TItMethod *itmethod, *itmethod_c, *itmethod_temp, *itmethod_conc;
  TItMethod *prec, *prec_temp, *prec_conc;

  TMatrix3D *MATRICES[12];

  TMatrix **matrices = (TMatrix **)MATRICES;

  TMatrix3D *matrixB1, *matrixB2, *matrixB3;
  TMatrix3D *matrixB1T, *matrixB2T, *matrixB3T;
  TMatrix3D *matrix_G11, *matrix_G22, *matrix_G33;
  TMatrix3D *matrix_tilde_G11, *matrix_tilde_G22, *matrix_tilde_G33;

  TMatrix3D **MatricesB1, **MatricesB2, **MatricesB3;
  TMatrix3D **MatricesB1T, **MatricesB2T, **MatricesB3T;
  TMatrix3D **Matrices_G11, **Matrices_G22, **Matrices_G33;
  TMatrix3D **Matrices_tilde_G11, **Matrices_tilde_G22, **Matrices_tilde_G33;

  TMGLevel3D *MGLevel_temp, *MGLevel_conc, *MGLevelGL00AuxProblem;

  TMultiGrid3D *MG_temp, *MG_conc, *MGGL00AuxProblem;

  TNSE_MGLevel *MGLevel;

  TNSE_MultiGrid *MG;

  TOutput3D *Output;

  TSquareMatrix3D *SQMATRICES[18];

  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;

  TSquareMatrix3D *mat;
  TSquareMatrix3D *matM, *matM_cons;
  TSquareMatrix3D *sqmatrixA;
  TSquareMatrix3D *sqmatrixA11, *sqmatrixA12, *sqmatrixA13;
  TSquareMatrix3D *sqmatrixA21, *sqmatrixA22, *sqmatrixA23;
  TSquareMatrix3D *sqmatrixA31, *sqmatrixA32, *sqmatrixA33;
  TSquareMatrix3D *sqmatrixGL00AuxProblem, *sqmatrixK, *sqmatrixL, *sqmatrixM;
  TSquareMatrix3D *sqmatrixM11, *sqmatrixM12, *sqmatrixM13;
  TSquareMatrix3D *sqmatrixM21, *sqmatrixM22, *sqmatrixM23;
  TSquareMatrix3D *sqmatrixM31, *sqmatrixM32, *sqmatrixM33;
  TSquareMatrix3D *sqmatrixPressSep;
  TSquareMatrix3D *sqmatrixA_temp, *sqmatrixM_temp, *sqmatrixK_temp;
  TSquareMatrix3D *sqmatrixA_conc, *sqmatrixM_conc, *sqmatrixK_conc;

  TSquareMatrix3D **MatricesGL00AuxProblem;
  TSquareMatrix3D **MatricesM11, **MatricesM12, **MatricesM13;
  TSquareMatrix3D **MatricesM21, **MatricesM22, **MatricesM23;
  TSquareMatrix3D **MatricesM31, **MatricesM32, **MatricesM33;
  TSquareMatrix3D **MatricesA, **MatricesK, **MatricesL, **MatricesM;
  TSquareMatrix3D **MatricesA11, **MatricesA12, **MatricesA13;
  TSquareMatrix3D **MatricesA21, **MatricesA22, **MatricesA23;
  TSquareMatrix3D **MatricesA31, **MatricesA32, **MatricesA33;
  TSquareMatrix3D **MatricesA_c, **MatricesK_c, **MatricesM_c;
  TSquareMatrix3D **MatricesA_temp, **MatricesM_temp, **MatricesK_temp;
  TSquareMatrix3D **MatricesA_conc, **MatricesM_conc, **MatricesK_conc;
  TSquareMatrix3D **MatricesS_c, **MatricesS_temp, **MatricesS_conc;

  TSquareStructure3D *matrix_structure, *sqstructureA, *sqstructureC, *sqstructure_temp;
  TSquareStructure3D *sqstructure_conc, *sqstructureL, *sqstructurePressSep;

  TStructure3D *structureB, *structureBT;
  TStructure3D *structure_G, *structure_tilde_G;

  // strings
  char AuxProblemString[] = "AuxProblem";
  char temp_String[] = "temp";
  char conc_String[] = "conc";
  char DivergenceString[] = "divergence";
  char DString[] = "d";
  char Mass[] = "Mass";
  char MassMatrix[] = "Mass matrix";
  char Name[] = "name";
  char NameString[] = "name";
  char PConvString[] = "p_conv";
  char PsepString[] = "psep";
  char PsiString[] = "psi";
  char PString[] = "p";
  char ReadinDat[] = "readin.dat";
  char UConvString[] = "u_conv";
  char UString[] = "u";
  char VorticityString[] = "vorticity";
  char VortString[] = "vort";
  char SmallScaleString[] = "smallscales";
  char LabelSpaceString[] = "space";

  // file name for f_old
  // char file_name_f_old = "save_f_old";

  // char pointers
  char *GmvBaseName, *GnuBaseName, *GrapeBaseName, *MatlabBaseName;
  char *PsBaseName, *ReadGrapeBaseName, *VtkBaseName;
  char *SaveDataFileName,  *ReadDataFileName;
  char *MAP, *PRM, *GEO;

  struct mallinfo MALLINFO;

  std::ostringstream os;
  os << " ";

  //======================================================================
  // end of the declaration of the variables
  //======================================================================

  //======================================================================
  // read parameter file
  //======================================================================

  total_time = GetTime();
  if(argc>=2)
    ret=Domain->ReadParam(argv[1]);
  else
    ret=Domain->ReadParam(ReadinDat);

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

  if( TDatabase::ParamDB->DISCTYPE==2 )
  {
    OutPut("SDFEM does not work!" << endl);
    Error("SDFEM does not work!" << endl);
    exit(4711);
  }
  if( TDatabase::ParamDB->DISCTYPE==5 )
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

  //======================================================================
  // definitions for Navier-Stokes equations
  //======================================================================

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

  //======================================================================
  // allocation of arrays
  //======================================================================

  // parameter for the size of the arrays
  LEVELS = TDatabase::ParamDB->LEVELS;

  l2p = new double[LEVELS+1];
  h1p = new double[LEVELS+1];
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
      break;
  }                                               // endswitch NSTYPE

  // matrices for VMS_PROJECTION
  if ( (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION) ||
    (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION_EXPL) ||
    (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION_SD) )
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
  if ( (TDatabase::ParamDB->DISCTYPE==GL00_AUX_PROBLEM)
    ||(TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE==4) )
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
  if ( TDatabase::ParamDB->DISCTYPE==VMS_PROJECTION_SD )
  {
    TDatabase::ParamDB->DISCTYPE=VMS_PROJECTION;
    OutPut("Change internally DISCTYPE to " << TDatabase::ParamDB->DISCTYPE
      << " since correct DiscreteForms have been assigned"<< endl);
  }

  //======================================================================
  // read boundary parameterization and initialize coarse grid
  //======================================================================
  Domain->Init(PRM, GEO);
  Domain->ConvertQuadToTri(TDatabase::ParamDB->CONVERT_QUAD_TO_TRI);

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

  // Navier-Stokes
  Coefficients[0] = LinCoeffs;
  // temp
  Coefficients[1] = BilinearCoeffs_temp;          //NOTE
  // conc
  Coefficients[2] = BilinearCoeffs_conc;
  Coefficients[3] = NoCoeffs;

  // create coarsest grid for multigrid method
  for ( i=0 ; i<TDatabase::ParamDB->SC_COARSEST_LEVEL_SADDLE ; i++ )
  {
    Domain->RegRefineAll();
  }

  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
  omega = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SADDLE;
  nonlinite = TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE;

  Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
  Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
  Parameters[2] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;

  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;

  if ( TDatabase::ParamDB->SOLVER_TYPE == GMG )
  {
    i=1;
    MG = new TNSE_MultiGrid(i, N_Paramters, Parameters);
  }

  // multigrid for Galdi/Layton model with auxiliary problem
  if ( TDatabase::ParamDB->DISCTYPE==GL00_AUX_PROBLEM )
  {
    MGGL00AuxProblem = new TMultiGrid3D(i, N_Paramters, Parameters);
  }

  mg_level = LEVELS+mg_level;

  // refine grids
  for ( i=0 ; i<mg_level ; i++ )
  {
    if ( i && (i<LEVELS) )
    {
      Domain->RegRefineAll();
    }
    if(TDatabase::ParamDB->WRITE_PS)
    {
      // write grid into an Postscript file
      os.seekp(std::ios::beg);
      os << PsBaseName << i << ".ps" << ends;
      Domain->PS(os.str().c_str(),It_EQ,i);
    }
  }

  //======================================================================
  // definitions for convection-reactions equations
  //======================================================================
  // force direct solver for cd equations
  solver_type_reaction = 2;
  ansatz_order_reaction = 1;

  // array for pointers to the solutions on the
  // different levels of the multigrid
  SolArray_temp = new TFEFunction3D*[LEVELS+1];
  SolArray_conc = new TFEFunction3D*[LEVELS+1];
  SolArray_temp_old = new TFEFunction3D*[LEVELS+1];
  SolArray_conc_old = new TFEFunction3D*[LEVELS+1];
  SolArray = new TFEFunction3D*[LEVELS+1];
  SolArray_old = new TFEFunction3D*[LEVELS+1];
  IntegralSpaces_conc_fct = new TFEFunction3D*[LEVELS+1];

  // array for pointers to right hand sides on the
  // different levels of the multigrid
  RhsArray_c = new double* [LEVELS+1];
  RhsArray_temp = new double* [LEVELS+1];
  RhsArray_conc = new double* [LEVELS+1];
  N_Array_temp = new int[LEVELS+1];
  N_Array_conc = new int[LEVELS+1];
  N_Array_c = new int[LEVELS+1];

  // array which points to the finite element spaces on the
  // different levels of the multigrid
  ConcentrationSpaces = new TFESpace3D*[LEVELS+1];
  ConcentrationSpaces_temp = new TFESpace3D*[LEVELS+1];
  ConcentrationSpaces_conc = new TFESpace3D*[LEVELS+1];
  IntegralSpaces_conc = new TFESpace3D*[LEVELS+1];

  // array which points to the system matrices on  the
  // different levels of the multigrid
  MatricesA_c = new TSquareMatrix3D*[LEVELS+1];
  MatricesA_temp = new TSquareMatrix3D*[LEVELS+1];
  MatricesA_conc = new TSquareMatrix3D*[LEVELS+1];

  // array which points to the mass matrices on  the
  // different levels of the multigrid
  MatricesM_c = new TSquareMatrix3D*[LEVELS+1];
  MatricesM_temp = new TSquareMatrix3D*[LEVELS+1];
  MatricesM_conc = new TSquareMatrix3D*[LEVELS+1];

  if (TDatabase::ParamDB->UREA_REACTION_DISC == SDFEM)
  {
    // array which points to the stabilization  matrices (sdfem) on the
    // different levels of the multigrid
    MatricesK_c = new TSquareMatrix3D*[LEVELS+1];
    MatricesK_temp = new TSquareMatrix3D*[LEVELS+1];
    MatricesK_conc = new TSquareMatrix3D*[LEVELS+1];

    if (TDatabase::ParamDB->SOLD_TYPE)
    {
      // array which points to the stabilization  matrices (sold) on the
      // different levels of the multigrid
      MatricesS_c = new TSquareMatrix3D*[LEVELS+1];
      MatricesS_temp = new TSquareMatrix3D*[LEVELS+1];
      MatricesS_conc = new TSquareMatrix3D*[LEVELS+1];
    }
  }
  if ((TDatabase::ParamDB->UREA_REACTION_DISC == GALERKIN) &&
    (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
  {
    // array which points to the stabilization  matrices (sdfem) on the
    // different levels of the multigrid
    MatricesK_c = new TSquareMatrix3D*[LEVELS+1];
    MatricesK_temp = new TSquareMatrix3D*[LEVELS+1];
    MatricesK_conc = new TSquareMatrix3D*[LEVELS+1];
    TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME = 1;
    // this is only for the PSD equation
    //TDatabase::ParamDB->INTERNAL_SORT_AMG = 2;
  }

  // pointers to the routines which compute matrix-vector
  // products and the defect
  MatVectScalar = MatVect_Scalar;
  DefectScalar = Defect_Scalar;

  //======================================================================
  // initialize discrete forms
  //======================================================================

  // discrete form for assembling mass matrix and rhs (Galerkin FEM)
  DiscreteFormMatrixMUrea = new TDiscreteForm3D
    (MassMatrix, Mass, N_Terms_MatrixM_Urea, Derivatives_MatrixM_Urea,
    SpacesNumbers_MatrixM_Urea, N_Matrices_MatrixM_Urea, N_Rhs_MatrixM_Urea,
    RowSpace_MatrixM_Urea, ColumnSpace_MatrixM_Urea, RhsSpace_MatrixM_Urea,
    MatrixMAssemble_Bulk3D, NoCoeffs, NULL);

  // discrete form for assembling stiffness matrix, stabilization matrix and rhs (SDFEM)
  DiscreteFormMatricesA_SUPG_Urea = new TDiscreteForm3D
    (MassMatrix, Mass, N_Terms_MatricesA_SUPG_Urea, Derivatives_MatricesA_SUPG_Urea,
    SpacesNumbers_MatricesA_SUPG_Urea, N_Matrices_MatricesA_SUPG_Urea, N_Rhs_MatricesA_SUPG_Urea,
    RowSpace_MatricesA_SUPG_Urea, ColumnSpace_MatricesA_SUPG_Urea, RhsSpace_MatricesA_SUPG_Urea,
    MatricesA_Assemble_SUPG_Bulk3D, BilinearCoeffs_temp, NULL);

  DiscreteFormMatricesA_SUPG_Urea_conc = new TDiscreteForm3D
    (MassMatrix, Mass, N_Terms_MatricesA_SUPG_Urea, Derivatives_MatricesA_SUPG_Urea,
    SpacesNumbers_MatricesA_SUPG_Urea, N_Matrices_MatricesA_SUPG_Urea, N_Rhs_MatricesA_SUPG_Urea,
    RowSpace_MatricesA_SUPG_Urea, ColumnSpace_MatricesA_SUPG_Urea, RhsSpace_MatricesA_SUPG_Urea,
    MatricesA_Assemble_SUPG_Bulk3D, BilinearCoeffs_conc, NULL);

  DiscreteFormRhs_SUPG_Urea = new TDiscreteForm3D
    (MassMatrix, Mass, N_Terms_Rhs_SUPG_Urea, Derivatives_Rhs_SUPG_Urea,
    SpacesNumbers_Rhs_SUPG_Urea, N_Matrices_Rhs_SUPG_Urea, N_Rhs_Rhs_SUPG_Urea,
    RowSpace_Rhs_SUPG_Urea, ColumnSpace_Rhs_SUPG_Urea, RhsSpace_Rhs_SUPG_Urea,
    Rhs_Assemble_SUPG_Bulk3D, BilinearCoeffs_temp, NULL);

  DiscreteFormRhs_SUPG_Urea_conc = new TDiscreteForm3D
    (MassMatrix, Mass, N_Terms_Rhs_SUPG_Urea, Derivatives_Rhs_SUPG_Urea,
    SpacesNumbers_Rhs_SUPG_Urea, N_Matrices_Rhs_SUPG_Urea, N_Rhs_Rhs_SUPG_Urea,
    RowSpace_Rhs_SUPG_Urea, ColumnSpace_Rhs_SUPG_Urea, RhsSpace_Rhs_SUPG_Urea,
    Rhs_Assemble_SUPG_Bulk3D, BilinearCoeffs_conc, NULL);

  // discrete form for assembling stiffness matrix, upwinding
  DiscreteFormMatricesA_Urea = new TDiscreteForm3D
    (MassMatrix, Mass, N_Terms_MatricesA_SUPG_Urea, Derivatives_MatricesA_SUPG_Urea,
    SpacesNumbers_MatricesA_SUPG_Urea, N_Matrices_MatricesA_Galerkin_Urea, N_Rhs_MatricesA_SUPG_Urea,
    RowSpace_MatricesA_Galerkin_Urea, ColumnSpace_MatricesA_Galerkin_Urea, RhsSpace_MatricesA_SUPG_Urea,
    MatricesA_Assemble_Bulk3D, BilinearCoeffs_temp, NULL);

  DiscreteFormMatricesA_Urea_conc = new TDiscreteForm3D
    (MassMatrix, Mass, N_Terms_MatricesA_SUPG_Urea, Derivatives_MatricesA_SUPG_Urea,
    SpacesNumbers_MatricesA_SUPG_Urea, N_Matrices_MatricesA_Galerkin_Urea, N_Rhs_MatricesA_SUPG_Urea,
    RowSpace_MatricesA_Galerkin_Urea, ColumnSpace_MatricesA_Galerkin_Urea, RhsSpace_MatricesA_SUPG_Urea,
    MatricesA_Assemble_Bulk3D, BilinearCoeffs_conc, NULL);

  DiscreteFormRhs_Urea =  new TDiscreteForm3D
    (MassMatrix, Mass, N_Terms_Rhs_Galerkin_Urea, Derivatives_Rhs_Galerkin_Urea,
    SpacesNumbers_Rhs_Galerkin_Urea, N_Matrices_Rhs_SUPG_Urea, N_Rhs_Rhs_SUPG_Urea,
    RowSpace_Rhs_SUPG_Urea, ColumnSpace_Rhs_SUPG_Urea, RhsSpace_Rhs_SUPG_Urea,
    Rhs_Assemble_Bulk3D, BilinearCoeffs_temp, NULL);

  DiscreteFormRhs_Urea_conc =  new TDiscreteForm3D
    (MassMatrix, Mass, N_Terms_Rhs_Galerkin_Urea, Derivatives_Rhs_Galerkin_Urea,
    SpacesNumbers_Rhs_Galerkin_Urea, N_Matrices_Rhs_SUPG_Urea, N_Rhs_Rhs_SUPG_Urea,
    RowSpace_Rhs_SUPG_Urea, ColumnSpace_Rhs_SUPG_Urea, RhsSpace_Rhs_SUPG_Urea,
    Rhs_Assemble_Bulk3D, BilinearCoeffs_conc, NULL);

  // discrete form for assembling stiffness matrix, Galerkin (FEM-FCT)
  DiscreteFormMatricesA_Galerkin_Urea = new TDiscreteForm3D
    (MassMatrix, Mass, N_Terms_MatricesA_SUPG_Urea, Derivatives_MatricesA_SUPG_Urea,
    SpacesNumbers_MatricesA_SUPG_Urea, N_Matrices_MatricesA_Galerkin_Urea, N_Rhs_MatricesA_SUPG_Urea,
    RowSpace_MatricesA_Galerkin_Urea, ColumnSpace_MatricesA_Galerkin_Urea, RhsSpace_MatricesA_SUPG_Urea,
    MatricesA_Assemble_Galerkin_Bulk3D, BilinearCoeffs_temp, NULL);

  DiscreteFormMatricesA_Galerkin_Urea_conc = new TDiscreteForm3D
    (MassMatrix, Mass, N_Terms_MatricesA_SUPG_Urea, Derivatives_MatricesA_SUPG_Urea,
    SpacesNumbers_MatricesA_SUPG_Urea, N_Matrices_MatricesA_Galerkin_Urea, N_Rhs_MatricesA_SUPG_Urea,
    RowSpace_MatricesA_Galerkin_Urea, ColumnSpace_MatricesA_Galerkin_Urea, RhsSpace_MatricesA_SUPG_Urea,
    MatricesA_Assemble_Galerkin_Bulk3D, BilinearCoeffs_conc, NULL);

  DiscreteFormRhs_Galerkin_Urea =  new TDiscreteForm3D
    (MassMatrix, Mass, N_Terms_Rhs_Galerkin_Urea, Derivatives_Rhs_Galerkin_Urea,
    SpacesNumbers_Rhs_Galerkin_Urea, N_Matrices_Rhs_SUPG_Urea, N_Rhs_Rhs_SUPG_Urea,
    RowSpace_Rhs_SUPG_Urea, ColumnSpace_Rhs_SUPG_Urea, RhsSpace_Rhs_SUPG_Urea,
    Rhs_Assemble_Bulk3D, BilinearCoeffs_temp, NULL);

  DiscreteFormRhs_Galerkin_Urea_conc =  new TDiscreteForm3D
    (MassMatrix, Mass, N_Terms_Rhs_Galerkin_Urea, Derivatives_Rhs_Galerkin_Urea,
    SpacesNumbers_Rhs_Galerkin_Urea, N_Matrices_Rhs_SUPG_Urea, N_Rhs_Rhs_SUPG_Urea,
    RowSpace_Rhs_SUPG_Urea, ColumnSpace_Rhs_SUPG_Urea, RhsSpace_Rhs_SUPG_Urea,
    Rhs_Assemble_Bulk3D, BilinearCoeffs_conc, NULL);

  // initialize multigrid
  if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
    Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
    i=1;
    MG_temp = new TMultiGrid3D(i, N_Paramters, Parameters);
    MG_conc = new TMultiGrid3D(i, N_Paramters, Parameters);
  }

  BoundaryConditions_Scalar[0] =  BoundCondition_temp;
  BoundaryConditions_Scalar[2] =  BoundCondition_conc;

  BoundValues_Scalar[0] = BoundValue_temp;
  BoundValues_Scalar[2] = BoundValue_conc;

  //======================================================================
  // definitions for population balance equation
  //======================================================================
  a_min = TDatabase::ParamDB->UREA_D_P_0/TDatabase::ParamDB->UREA_D_P_MAX;
  solver_type_reaction = 2;
  ansatz_order_reaction = 1;

  t3 = GetTime();
  total_time = t3 - total_time;

  //======================================================================
  // loop over all levels
  //======================================================================

  for ( i=0 ; i<mg_level ; i++ )
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

    if (i==mg_level-1)
      coll=Domain->GetCollection(It_Finest, 0);
    else
      coll=Domain->GetCollection(It_EQ, i+TDatabase::ParamDB->SC_COARSEST_LEVEL_SADDLE);
    if ((i==mg_level-2)&& (mg_type==2))
      coll=Domain->GetCollection(It_Finest, 0);

    OutPut("cells " << coll->GetN_Cells()<< endl);

    TDatabase::ParamDB->INTERNAL_LEVEL = 0;
    if (i== mg_level-1)
      TDatabase::ParamDB->INTERNAL_LEVEL = 1;

    // get spaces for low order disc on finest geo grid
    if (((mg_type==1)&&(i<mg_level-1))||((mg_type==2)&&(i<mg_level-2)))
    {
      velocity_space = new TFESpace3D(coll, NameString, UString, BoundCondition, Non_USpace, 1);
      pressure_space = new TFESpace3D(coll, NameString, PString, BoundCondition, DiscP_PSpace, 0);
      velocity_space_code = -1;
      pressure_space_code = 0;
      order = -1;
      convolution_space = new TFESpace3D(coll, NameString, UString, BoundConditionAuxProblem, Non_USpace, 1);
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

    // build fespace hierarchy
    // set values and pointers for fe space
    USpaces[i] = velocity_space;
    PSpaces[i] = pressure_space;
    N_U = velocity_space->GetN_DegreesOfFreedom();
    N_P = pressure_space->GetN_DegreesOfFreedom();
    N_Uarray[i] = N_U;
    N_Parray[i] = N_P;
    N_Active = velocity_space->GetActiveBound();

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

    // matrices for VMS_PROJECTION
    if ( (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION)||
      (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION_EXPL) )
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

    if (i==mg_level-1)
      sol_timestep_m1 =  new double[N_Unknowns];

    B = new double [N_Unknowns];

    // ( A B' )
    // ( B 0  )

    if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
    {
      // determine number of auxiliary arrays
      if ( (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE)
        || (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE) )
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

      for ( j=0 ; j<coll->GetN_Cells() ; j++ )
      {
        downwind[i][j] = j;
      }

#ifdef __DOWNWIND__
      DownwindNumberingCells(coll, downwind[i]);
#endif

      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          MGLevel = new TNSE_MGLevel1(i, sqmatrixM, matrixB1, matrixB2, matrixB3, structureBT, B, sol, n_aux,
            alpha, velocity_space_code, pressure_space_code,NULL,downwind[i]);
          break;
        case 2:
          MGLevel = new TNSE_MGLevel2(i, sqmatrixM, matrixB1, matrixB2, matrixB3,  matrixB1T, matrixB2T, matrixB3T,
            B, sol, n_aux, alpha, velocity_space_code, pressure_space_code,NULL,downwind[i]);
          break;
        case 3:
          MGLevel = new TNSE_MGLevel3(i, sqmatrixM11, sqmatrixM12, sqmatrixM13, sqmatrixM21, sqmatrixM22, sqmatrixM23,
            sqmatrixM31, sqmatrixM32, sqmatrixM33, matrixB1, matrixB2, matrixB3,
            structureBT, B, sol, n_aux, alpha, velocity_space_code,  pressure_space_code,NULL,downwind[i]);
          break;
        case 4:
          MGLevel = new TNSE_MGLevel4(i, sqmatrixM11, sqmatrixM12, sqmatrixM13, sqmatrixM21, sqmatrixM22, sqmatrixM23, sqmatrixM31,
            sqmatrixM32, sqmatrixM33, matrixB1, matrixB2, matrixB3, matrixB1T, matrixB2T, matrixB3T,
            B, sol, n_aux, alpha, velocity_space_code, pressure_space_code,NULL,downwind[i]);
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
    if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE==4)

    {
      uConvSpaces[i] = convolution_space;
      // define vector fe function for convolution
      // allocate memory for values of convoluted function
      u_uConv = new double[3*N_U];
      // initialize u_conv to 0
      memset(u_uConv,0,3*N_U*SizeOfDouble);
      // allocate fe vector function for convolution
      uConv = new TFEVectFunct3D(convolution_space, UConvString, UConvString, u_uConv, N_U, 3);
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

     // read initial solution of finest level from grape file
    if ((i==mg_level-1)&&(TDatabase::ParamDB->READ_DATA>=1))
    {
      save_sol[0] = sol;
      save_N_Unknowns[0] = N_Unknowns;
      //ReadData(ReadDataFileName,1,save_sol,save_N_Unknowns);
      ReadData(ReadGrapeBaseName,1,save_sol,save_N_Unknowns);
    }

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
      default:
        OutPut("Unknown DISCTYPE " << TDatabase::ParamDB->DISCTYPE << endl);
        exit(1);
    }

    // parameters which are the same for all NSTYPEs
    N_Rhs = 3;
    N_FESpaces = 3;

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
          N_SquareMatrices = 13;
          SQMATRICES[12] = MatricesGL00AuxProblem[i];
          SQMATRICES[12]->Reset();
        }
        if ((CurrentDiscType == VMS_PROJECTION)|| (CurrentDiscType == VMS_PROJECTION_EXPL))
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
        break;
    }

    t1 = GetTime();

    // set rhs
    fesp[0] = velocity_space;
    fesp[1] = pressure_space;
    fesp[2] = convolution_space;

    fefct[0] = u1;
    fefct[1] = u2;
    fefct[2] = u3;
    fefct[3] = du11ConvArray[i];
    fefct[4] = du12ConvArray[i];
    fefct[5] = du13ConvArray[i];
    fefct[6] = du22ConvArray[i];
    fefct[7] = du23ConvArray[i];
    fefct[8] = du33ConvArray[i];

    ferhs[0] = velocity_space;
    ferhs[1] = velocity_space;
    ferhs[2] = velocity_space;

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
        case GL00_AUX_PROBLEM:
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
          // no parameter necessary
        case UPWIND:
          aux = new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
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
    OutPut("time for assembling " << t2-t1 << "s" << endl);
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
    // update matrices
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
          OutPut("update ");
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

    // declarations for reaction equations
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

    // if multiple discretization multilevel method is used
    // get space for low order disc on finest geo grid
    if ((mg_type==1)&&(i<mg_level-1))
    {
      // nonconforming velocity space
      concentration_space = new TFESpace3D(coll,NameString,UString,BoundCondition_temp,-1);
      concentration_space_c = new TFESpace3D(coll,NameString,UString,BoundCondition_conc,-1);
    }
    // standard multigrid or finest level
    // get fe space of high order disc on finest geo grid
    else
    {
      concentration_space = new TFESpace3D(coll, NameString,
        UString, BoundCondition_temp, ansatz_order_reaction);
      concentration_space_c = new TFESpace3D(coll, NameString,
        UString, BoundCondition_conc, ansatz_order_reaction);
    }
    // array of the fe spaces
    ConcentrationSpaces_temp[i] = concentration_space;
    ConcentrationSpaces_conc[i] = concentration_space_c;
    N_Unknowns_temp = concentration_space->GetN_DegreesOfFreedom();
    N_Unknowns_conc = concentration_space_c->GetN_DegreesOfFreedom();
    N_Array_temp[i] = N_Unknowns_temp;
    N_Array_conc[i] = N_Unknowns_conc;

    //N_Unknowns_conc = N_Unknowns_temp;

    // active dof (i.e. dof without Dirichlet dofs)
    N_Active_temp = concentration_space->GetActiveBound();
    N_Active_conc = concentration_space_c->GetActiveBound();

    OutPut("degrees of freedom (temp): "<< setw(10) << N_Unknowns_temp << " active: " <<
      N_Active_temp << endl);
    OutPut("degrees of freedom (conc): "<< setw(10) << N_Unknowns_conc << " active: " <<
      N_Active_conc << endl);

    // **********************************************************************************
    // substance A
    // **********************************************************************************
    // build matrices
    // first build matrix structure
    sqstructure_temp = new TSquareStructure3D(ConcentrationSpaces_temp[i]);
    sqstructure_temp->Sort();

    // two matrices used
    // A contains the non time dependent part of the discretization
    sqmatrixA = new TSquareMatrix3D(sqstructure_temp);
    MatricesA_temp[i] = sqmatrixA;

    // M is the mass matrix
    // the iterative solver uses M
    sqmatrixM_temp = new TSquareMatrix3D(sqstructure_temp);
    MatricesM_temp[i] = sqmatrixM_temp;

    if (TDatabase::ParamDB->UREA_REACTION_DISC == SDFEM)
    {
      // stabilisation matrix K
      sqmatrixK_temp = new TSquareMatrix3D(sqstructure_temp);
      MatricesK_temp[i] = sqmatrixK_temp;
    }

    // allocate array for solution
    sol_temp = new double[N_Unknowns_temp];
    memset(sol_temp, 0, N_Unknowns_temp*SizeOfDouble);
    current_sol = sol_temp;
    oldsol_temp = new double[N_Unknowns_temp];

    // allocate fe function on the current level
    temp = new TFEFunction3D(concentration_space, temp_String, temp_String, sol_temp, N_Unknowns_temp);
    SolArray_temp[i] = temp;

    // allocate array on which the time scheme works
    B_c = new double [N_Unknowns_temp];
    current_B = B_c;

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
      MGLevel_temp = new TMGLevel3D(i, sqmatrixM_temp, current_B, current_sol,  n_aux, NULL);
      MG_temp->AddLevel(MGLevel_temp);
    }
    // interpolate initial condition
    temp->Interpolate(InitialCondition_temp);

    // **********************************************************************************
    // concentration
    // **********************************************************************************

    // build matrices
    // first build matrix structure
    sqstructure_conc = new TSquareStructure3D(ConcentrationSpaces_conc[i]);
    sqstructure_conc->Sort();

    // two matrices used
    // A contains the non time dependent part of the discretization
    sqmatrixA = new TSquareMatrix3D(sqstructure_conc);
    MatricesA_conc[i] = sqmatrixA;

    // M is the mass matrix
    // the iterative solver uses M
    sqmatrixM_conc = new TSquareMatrix3D(sqstructure_conc);
    MatricesM_conc[i] = sqmatrixM_conc;

    if (TDatabase::ParamDB->UREA_REACTION_DISC == SDFEM)
    {
      // stabilisation matrix K
      sqmatrixK_conc = new TSquareMatrix3D(sqstructure_conc);
      MatricesK_conc[i] = sqmatrixK_conc;
    }

    // allocate array for solution
    sol_conc = new double[N_Unknowns_conc];
    memset(sol_conc, 0, N_Unknowns_conc*SizeOfDouble);
    current_sol = sol_conc;
    oldsol_conc = new double[N_Unknowns_conc];

    // allocate fe function on the current level
    conc = new TFEFunction3D(concentration_space_c, conc_String, conc_String, sol_conc, N_Unknowns_conc);
    SolArray_conc[i] = conc;

    // allocate array on which the time scheme works
    B_c = new double [N_Unknowns_conc];
    current_B = B_c;

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
      MGLevel_conc = new TMGLevel3D(i, sqmatrixM_conc, current_B, current_sol, n_aux, NULL);
      MG_conc->AddLevel(MGLevel_conc);
    }
    // interpolate initial condition
    conc->Interpolate(InitialCondition_conc);

    // define bilinear fe space which contains the integrals of the decrease of f
    integral_space_conc = new TFESpace3D(coll, NameString, UString, BoundCondition_conc, 1);
    IntegralSpaces_conc[i] = integral_space_conc;
    N_Unknowns_Integral_Space = integral_space_conc->GetN_DegreesOfFreedom();
    integral_val = new double[N_Unknowns_Integral_Space];
    memset(integral_val,0, N_Unknowns_Integral_Space*SizeOfDouble);
    integral_space_conc_fct = new TFEFunction3D(integral_space_conc, temp_String, temp_String,
      integral_val,N_Unknowns_Integral_Space);
    IntegralSpaces_conc_fct[i] = integral_space_conc_fct;

    // **********************************************************************************
    // general definitions
    // **********************************************************************************

    // allocate rhs
    rhs_c = new double[N_Unknowns_temp];
    memset(rhs_c, 0, N_Unknowns_temp*SizeOfDouble);
    RhsArray_temp[i] = rhs_c;
    rhs_c = new double[N_Unknowns_conc];
    memset(rhs_c, 0, N_Unknowns_conc*SizeOfDouble);
    RhsArray_conc[i] = rhs_c;

    // prepare output
    if (i==mg_level-1)
    {
      if ((TDatabase::ParamDB->WRITE_GRAPE)||(TDatabase::ParamDB->WRITE_GMV)||(TDatabase::ParamDB->WRITE_VTK))
      {
        // output velocity, pressure, temp, conc
        Output = new TOutput3D(5, 4, 1, 2, Domain);
        Output->AddFEVectFunct(u);
        Output->AddFEFunction(p);
        Output->AddFEFunction(temp);
        Output->AddFEFunction(conc);
        os.seekp(std::ios::beg);
        Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());
        os.seekp(std::ios::beg);
        Output->AddParameter(real_time,os.str().c_str());
      }
    }                                             // end of preparing output

    // read initial solution of finest level from grape file
    if ((i==mg_level-1)&&(TDatabase::ParamDB->READ_GRAPE_FILE))
    {
      // only velocity and pressure is read
      AuxFEFunctArray = new TFEFunction3D*[4];
      AuxFEFunctArray[0] = PArray[mg_level-1];
      AuxFEFunctArray[1] = temp;
      AuxFEFunctArray[3] = conc;
      AuxFEVectFunctArray = new TFEVectFunct3D*[1];
      AuxFEVectFunctArray[0] = UArray[mg_level-1];
      ReadGrapeFile3D(ReadGrapeBaseName, 4 , 1, AuxFEFunctArray,AuxFEVectFunctArray);
      if (TDatabase::TimeDB->RESET_CURRENTTIME > 0)
      {
        TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->RESET_CURRENTTIME_STARTTIME;
        OutPut("start time reset to " << TDatabase::TimeDB->CURRENTTIME << endl);
        memset(sol_temp, 0, N_Unknowns_temp*SizeOfDouble);
        OutPut("temperature set to zero" << endl);
        memset(sol_conc, 0, N_Unknowns_conc*SizeOfDouble);
        OutPut("concentrations set to zero" << endl);
      }
    }
    // read initial solution of finest level
    // only temperature and concentration
    // flow field is just dummy
    if ((i==mg_level-1)&&(TDatabase::ParamDB->READ_DATA>=2))
    {
      save_sol[0] = oldsol;
      save_sol[1] = sol_temp;
      save_sol[2] = sol_conc;
      save_N_Unknowns[0] = N_Unknowns;
      save_N_Unknowns[1] = N_Unknowns_temp;
      save_N_Unknowns[2] = N_Unknowns_conc;
      ReadData(ReadDataFileName,3,save_sol,save_N_Unknowns);
      memcpy(oldsol, sol, N_Unknowns*SizeOfDouble);
      //OutPut("setting concentration to zero " << endl);
      //memset(sol_conc,0,N_Unknowns_conc*SizeOfDouble);
    }

    //======================================================================
    // assembling of mass matrices
    //======================================================================

    // set parameters
    N_Rhs = 0;
    N_FESpaces = 1;
    N_SquareMatrices = 1;
    DiscreteForm = DiscreteFormMatrixMUrea;
    aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

    // substance A
    fesp[0] = ConcentrationSpaces_temp[i];
    SQMATRICES[0] = MatricesM_temp[i];
    SQMATRICES[0]->Reset();

    if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
      (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
      BoundValues_Scalar[0] = BoundValue_FEM_FCT;

    Assemble3D(N_FESpaces, fesp,
      N_SquareMatrices, SQMATRICES,
      0, NULL,
      0, NULL, NULL,
      DiscreteForm,
      BoundaryConditions_Scalar,
      BoundValues_Scalar,
      aux);

    if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
      (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
      BoundValues_Scalar[0] = BoundValue_temp;

    if ((TDatabase::ParamDB->UREA_REACTION_DISC == GALERKIN) &&
      (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)
      && (i==mg_level-1))
    {
      CheckWrongNeumannNodes_temp(coll, ConcentrationSpaces_temp[i], N_neum_to_diri_temp, neum_to_diri_temp,
        neum_to_diri_temp_x, neum_to_diri_temp_y, neum_to_diri_temp_z);
    }
    if ((TDatabase::ParamDB->UREA_REACTION_DISC == GALERKIN) &&
      (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)
      && (i==mg_level-1))
    {
      oldrhs_fem_fct0_temp = new double[N_Unknowns_temp];
      memcpy(oldrhs_fem_fct0_temp,  RhsArray_temp[i], N_Unknowns_temp*SizeOfDouble);
      oldrhs_fem_fct1_temp = new double[N_Unknowns_temp];
      lump_mass_temp = new double [N_Unknowns_temp];
      LumpMassMatrixToVector(MatricesM_temp[i], lump_mass_temp);
      matrix_D_Entries_temp = new double[MatricesA_temp[i]->GetN_Entries()];
      // matrix K for copy of mass matrix
      sqmatrixK_temp = new TSquareMatrix3D(sqstructure_temp);
      MatricesK_temp[i] = sqmatrixK_temp;
      tilde_u_temp = new double [N_Unknowns_temp];
      // save mass matrix in matricesK
      memcpy(MatricesK_temp[i]->GetEntries(), MatricesM_temp[i]->GetEntries(),
        MatricesM_temp[i]->GetN_Entries() * SizeOfDouble);
    }

    // if (TDatabase::ParamDB->UREA_REACTION_MASS_LUMPING)
    // {
    //   LumpMassMatrixToDiag_Urea(MatricesM_temp[i]);
    // }

    // substance C
    fesp[0] = ConcentrationSpaces_conc[i];
    SQMATRICES[0] = MatricesM_conc[i];
    SQMATRICES[0]->Reset();

    if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
      (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
      BoundValues_Scalar[2] = BoundValue_FEM_FCT;

    Assemble3D(N_FESpaces, fesp,
      N_SquareMatrices, SQMATRICES,
      0, NULL,
      0, NULL, NULL,
      DiscreteForm,
      BoundaryConditions_Scalar+2,
      BoundValues_Scalar+2,
      aux);

    if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
      (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
      BoundValues_Scalar[2] = BoundValue_conc;

    if ((TDatabase::ParamDB->UREA_REACTION_DISC == GALERKIN) &&
      (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)
      && (i==mg_level-1))
    {
      CheckWrongNeumannNodes_conc(coll, ConcentrationSpaces_conc[i], N_neum_to_diri_conc, neum_to_diri_conc,
        neum_to_diri_conc_x, neum_to_diri_conc_y, neum_to_diri_conc_z);
    }
    if ((TDatabase::ParamDB->UREA_REACTION_DISC == GALERKIN) &&
      (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)
      && (i==mg_level-1))
    {
      oldrhs_fem_fct0_conc = new double[N_Unknowns_conc];
      memcpy(oldrhs_fem_fct0_conc,  RhsArray_conc[i], N_Unknowns_conc*SizeOfDouble);
      oldrhs_fem_fct1_conc = new double[N_Unknowns_conc];
      lump_mass_conc = new double [N_Unknowns_conc];
      LumpMassMatrixToVector(MatricesM_conc[i], lump_mass_conc);
      matrix_D_Entries_conc = new double[MatricesA_conc[i]->GetN_Entries()];
      // matrix K for copy of mass matrix
      sqmatrixK_conc = new TSquareMatrix3D(sqstructure_conc);
      MatricesK_conc[i] = sqmatrixK_conc;
      tilde_u_conc= new double [N_Unknowns_conc];
      memcpy(MatricesK_conc[i]->GetEntries(), MatricesM_conc[i]->GetEntries(),
        MatricesM_conc[i]->GetN_Entries() * SizeOfDouble);
    }

    delete aux;
  }                                               // endfor i  loop over all cells (line 1044)

  t4 =  GetTime();
  total_time += t4 - t3;
  t3 = t4;

  //======================================================================
  // end of space cycle, finest grid reached
  // everything happens on the same grid
  //======================================================================
  // definitions for Navier-Stokes equations
  // copy sol for extrapolation after time step
  memcpy(sol_timestep_m1,sol,N_Unknowns*SizeOfDouble);
  comp_vort =0;

  if ((TDatabase::ParamDB->WRITE_GRAPE)||(TDatabase::ParamDB->WRITE_GMV)||
    (TDatabase::ParamDB->WRITE_VTK))
  {
    if (TDatabase::ParamDB->WRITE_GMV)
    {
      os.seekp(std::ios::beg);
      os << GmvBaseName << 0 << ".gmv" << ends;
      Output->WriteGMV(os.str().c_str());
    }
    if (TDatabase::ParamDB->WRITE_VTK)
    {
      os.seekp(std::ios::beg);
      os << VtkBaseName << 0 << ".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      //os.seekp(std::ios::beg);
      //os << VtkBaseName << "end.psd." << m << ".vtk" << ends;
      //write_vtk_file(Nx, Ny, Nz, Na, cut_coord, x_coord, y_coord, z_coord, a_coord, sol_psd,os.str().c_str());
    }
    if (TDatabase::ParamDB->WRITE_GRAPE)
    {
      os.seekp(std::ios::beg);
      os << GrapeBaseName << 0 << ".dat" << ends;
      Output->WriteGrape(os.str().c_str());
    }
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
  {
    newton = new double[3*N_Uarray[mg_level-1]];
  }

  N_Active = velocity_space->GetActiveBound();

  solver_time = 0.0;
  N_LinIter = 0;

  gamma = 0;
  gamma_c = 0;
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

  // definitions for convection-reaction equations and
  // population balance equation

  // allocate arrays for solver
  //sol_c = new double[N_Unknowns_temp];
  defect_c = new double[N_Unknowns_temp];
  startsol_c = new double[N_Unknowns_temp];
  rhs_c_complete_A  = new double[N_Unknowns_temp];
  oldrhs_temp =  new double[N_Unknowns_temp];
  oldrhs_conc =  new double[N_Unknowns_conc];
  memset(oldrhs_temp, 0, N_Unknowns_temp*SizeOfDouble);
  memset(oldrhs_conc, 0, N_Unknowns_conc*SizeOfDouble);

  memcpy(oldsol_temp,sol_temp,N_Unknowns_temp*SizeOfDouble);
  memcpy(oldsol_conc,sol_conc,N_Unknowns_conc*SizeOfDouble);

  // number of active d.o.f.
  N_Active_temp = ConcentrationSpaces_temp[mg_level-1]->GetActiveBound();
  N_Unknowns_temp = ConcentrationSpaces_temp[mg_level-1]->GetN_DegreesOfFreedom();

  // initialize solver
  if (solver_type_reaction==1)
  {
    low_scalar = 0;
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
        prec_temp = new TMultiGridScaIte(MatVectScalar, DefectScalar, NULL,
          0, N_Unknowns_temp, MG_temp, zerostart);
        prec_conc = new TMultiGridScaIte(MatVectScalar, DefectScalar, NULL,
          0, N_Unknowns_conc, MG_conc, zerostart);
        break;
      default:
        OutPut("Unknown preconditioner !!!" << endl);
        exit(4711);
    }
    switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
    {
      // fixed point iteration
      case 11:
        itmethod_temp = new TFixedPointIte(MatVectScalar, DefectScalar, prec_temp,
          0, N_Unknowns_temp, 0);
        itmethod_conc = new TFixedPointIte(MatVectScalar, DefectScalar, prec_conc,
          0, N_Unknowns_conc, 0);
        if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
        {
          itmethod_sol_temp = new double[N_Unknowns_temp];
          itmethod_sol_conc = itmethod_sol_temp;
          itmethod_rhs_c = new double[N_Unknowns_temp];
        }
        else
        {
          itmethod_sol_temp = sol_temp;
          itmethod_sol_conc = sol_conc;
          itmethod_rhs_c = rhs;
        }
        break;
      case 16:
        // FGMRES
        itmethod_temp = new TFgmresIte(MatVectScalar, DefectScalar, prec_temp,
          0, N_Unknowns_temp, 1);
        itmethod_conc = new TFgmresIte(MatVectScalar, DefectScalar, prec_conc,
          0, N_Unknowns_conc, 1);
        if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
        {
          itmethod_sol_temp = new double[N_Unknowns_temp];
          itmethod_sol_conc = itmethod_sol_temp;
          itmethod_rhs_c = new double[N_Unknowns_temp];
        }
        else
        {
          itmethod_sol_temp = sol_temp;
          itmethod_sol_conc = sol_conc;
          itmethod_rhs_c = rhs;
        }
        break;
      default:
        OutPut("Unknown solver !!!" << endl);
        exit(4711);
    }
  }
  else
  {
    low_scalar = mg_level - 1;
  }

  // definitions for turbulence models
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

  // declarations for the population balance equation
  //number of layers in a-direction
  Na = TDatabase::ParamDB->N_CELL_LAYERS_PSD;
  a_layers_coord = new double[Na+1];
  memset(a_layers_coord,0,(Na+1)*SizeOfDouble);
  // generation of the grid
  // outputs are Nx, Ny, Nz, x_coord, y_coord, z_coord, z_layers_coord
  OutPut("grid generation for population balance equation " << endl);

  grid_generator_4d(coll, Nx, Ny, Nz,
    x_min, x_max, y_min, y_max, z_min, z_max,
    a_min, a_max, Na,
    x_coord, y_coord, z_coord, a_coord,
    a_layers_coord);

  OutPut("grid for population balance: " << Nx << " x " << Ny << " x "
	 << Nz << " x " << Na << " -- nodes: " << (Nx+1)*(Ny+1)*(Nz+1)*(Na+1)<< endl);
  OutPut("size of 3D domain: ["<<x_min<<","<<x_max<<"] x ["
    <<y_min<<","<<y_max<<"] x [" <<z_min<<","<<z_max<<"]" <<endl);

  // total number of nodes or grid points
  Nodes = (Nx+1)*(Ny+1)*(Nz+1)*(Na+1);
  if (TDatabase::ParamDB->UREA_PB_DISC== UREA_FEM_FCT)
  {
      sol_psd_help = new double[Nodes];
      memset(sol_psd_help, 0, Nodes*SizeOfDouble);
  }
  sol_psd = new double[Nodes];
  memset(sol_psd, 0, Nodes*SizeOfDouble);
  sol_psd_old = new double[Nodes];
  memset(sol_psd_old, 0, Nodes*SizeOfDouble);
  rhs_psd = new double[Nodes];
  memset(rhs_psd, 0, Nodes*SizeOfDouble);
  rhs_psd_old = new double[Nodes];
  memset(rhs_psd_old, 0, Nodes*SizeOfDouble);
  if (TDatabase::ParamDB->READ_DATA>=3)
  {
      save_sol[0] = sol;
      save_sol[1] = sol_temp;
      save_sol[2] = sol_conc;
      save_sol[3] = sol_psd;
      save_N_Unknowns[0] = N_Unknowns;
      save_N_Unknowns[1] = N_Unknowns_temp;
      save_N_Unknowns[2] = N_Unknowns_conc;
      save_N_Unknowns[3] = Nodes;
      ReadData(ReadDataFileName,4,save_sol,save_N_Unknowns);
  }
  OutPut("psd " << Ddot(Nodes,sol_psd,sol_psd) << endl);

  velo1 = new double[(Nx+1)*(Ny+1)*(Nz+1)];
  memset(velo1, 0, (Nx+1)*(Ny+1)*(Nz+1)*SizeOfDouble);
  velo1[0] = -4711;
  velo2 = new double[(Nx+1)*(Ny+1)*(Nz+1)];
  memset(velo2, 0, (Nx+1)*(Ny+1)*(Nz+1)*SizeOfDouble);
  velo3 = new double[(Nx+1)*(Ny+1)*(Nz+1)];
  memset(velo3, 0, (Nx+1)*(Ny+1)*(Nz+1)*SizeOfDouble);
  concent_C_array = new double[(Nx+1)*(Ny+1)*(Nz+1)];
  memset(concent_C_array, 0, (Nx+1)*(Ny+1)*(Nz+1)*SizeOfDouble);
  Temp_array = new double[(Nx+1)*(Ny+1)*(Nz+1)];
  memset(Temp_array, 0, (Nx+1)*(Ny+1)*(Nz+1)*SizeOfDouble);

  limit_c = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR;
  Max_It_c = TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR;

  if ((TDatabase::ParamDB->UREA_PB_DISC==UREA_BWE_FDM_UPWIND) ||
    (TDatabase::ParamDB->UREA_PB_DISC==UREA_FEM_FCT))
  {
    // declaration and initialisation of the column and the row pointer
    col_ptr = new int[81*Nodes];
    memset(col_ptr,0,(81*Nodes)*SizeOfInt);
    row_ptr = new int[Nodes+1];
    memset(row_ptr,0,(Nodes+1)*SizeOfInt);

    // filling of the column and the row pointer
    filling_row_and_col_ptr(&N_Entries, Nodes, Nx, Ny, Nz, x_max, x_min, y_max,
      y_min, z_max, z_min, a_max, a_min, x_coord, y_coord, z_coord, a_coord,
      row_ptr, col_ptr);

    OutPut("numerate collection" << endl);

    // number of mesh cells -> necessary for the corresponding 2d grid
    N_Cells = coll->GetN_Cells();

    // declaration, initialisation and computing of the corresponding 2d grid
    correspond_3dgrid = new int[N_Cells];
    memset(correspond_3dgrid,0,N_Cells*SizeOfInt);

    generate_correspond_3d_grid(Nx, Ny, Nz, x_coord, y_coord, z_coord, coll, correspond_3dgrid);
    OutPut("matrix structure " << Nodes <<  " " << N_Entries << endl);
    matrix_structure = new TSquareStructure3D(Nodes, N_Entries, col_ptr, row_ptr);
    OutPut("matrix structure done" << endl);
    mat = new TSquareMatrix3D(matrix_structure);

    if((TDatabase::ParamDB->UREA_PB_DISC_STAB == GALERKIN) &&
      (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
    {
      matM = new TSquareMatrix3D(matrix_structure);
      matM_cons = new TSquareMatrix3D(matrix_structure);
      lump_mass_PSD = new double[Nodes];
      matrix_D_Entries_PSD = new double[mat->GetN_Entries()];
      OutPut("lump_mass_PSD " << Nodes << endl);
    }
    else
    {
      matM = NULL;
      matM_cons = NULL;
      lump_mass_PSD = NULL;
      matrix_D_Entries_PSD = NULL;
    }
  }

  if (TDatabase::ParamDB->UREA_PB_DISC==UREA_FWE_FDM_UPWIND)
  {
    // number of mesh cells -> necessary for the corresponding 2d grid
    N_Cells = coll->GetN_Cells();

    // declaration, initialisation and computing of the corresponding 2d grid
    correspond_3dgrid = new int[N_Cells];
    memset(correspond_3dgrid,0,N_Cells*SizeOfInt);

    generate_correspond_3d_grid(Nx, Ny, Nz, x_coord, y_coord, z_coord, coll, correspond_3dgrid);
  }

  /*
    // write data for pictures
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

  if(TDatabase::ParamDB->WRITE_VTK)
  {
  os.seekp(std::ios::beg);
  os << VtkBaseName << 0 << ".vtk" << ends;
  Output->WriteVtk(os.str().c_str());
  os.seekp(std::ios::beg);
  os << VtkBaseName << "end.psd." << m << ".vtk" << ends;
  write_vtk_file(Nx, Ny, Nz, Na, cut_coord, x_coord, y_coord, z_coord, a_coord, sol_psd,os.str().c_str());
  }
  */
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
  total_time = GetTime();
  while(TDatabase::TimeDB->CURRENTTIME< end_time)
  {                                               // time cycle
    OutPut("time cycle memory : " << setw(10) << GetMemory() << endl);
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

      // sub steps of fractional step theta
      for(l=0;l<N_SubSteps;l++)
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
       if ((TDatabase::TimeDB->CURRENTTIME > 6 - 1e-3) && (TDatabase::TimeDB->CURRENTTIME < 6 + 1e-3))
       {  
          TDatabase::TimeDB->TIMESTEPLENGTH *=2.0;
          TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = TDatabase::TimeDB->TIMESTEPLENGTH;
       }
       
         if (TDatabase::TimeDB->CURRENTTIME > 10)
        {
          TDatabase::TimeDB->TIMESTEPLENGTH = 0.1;
          TDatabase::TimeDB->CURRENTTIMESTEPLENGTH = TDatabase::TimeDB->TIMESTEPLENGTH;
       }
 
	 tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
        if (!very_first_time)
          TDatabase::TimeDB->CURRENTTIME += tau;
	//if (!steady)
	//{
        // working array for rhs is B, initialize B
        memset(B, 0, N_Unknowns*SizeOfDouble);

        // old rhs array (f) multiplied with current subtime step and theta3 on B
        Daxpy(N_Active, tau*theta3, rhs, B);
        Daxpy(N_Active, tau*theta3, rhs+N_U, B+N_U);
        Daxpy(N_Active, tau*theta3, rhs+2*N_U, B+2*N_U);

        //======================================================================
        // prepare input data (parameters) for several discretizations

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
            N_FESpaces = 2;
            fesp[0] = USpaces[mg_level-1];
            fesp[1] = ProjectionSpaces[mg_level-1];

            fefct[0] = U1Array[mg_level-1];
            fefct[1] = U2Array[mg_level-1];
            fefct[2] = U3Array[mg_level-1];

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

            MATRICES[0] = Matrices_tilde_G11[mg_level-1];
            MATRICES[1] = Matrices_tilde_G22[mg_level-1];
            MATRICES[2] = Matrices_tilde_G33[mg_level-1];
            MATRICES[0]->Reset();
            MATRICES[1]->Reset();
            MATRICES[2]->Reset();

            Assemble3D(N_FESpaces, fesp,
              0, NULL,
              3, MATRICES,
              0, NULL, ferhs,
              DiscreteFormNLVMS_ProjectionExpl,
              BoundaryConditions,
              BoundValues,
              aux);
            delete aux;

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
        if (TDatabase::ParamDB->DISCTYPE==GL00_AUX_PROBLEM)
          Daxpy(3*N_U, tau, LESModelRhs, B);

        // update rhs in explicit projection-based VMS
        if (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION_EXPL)
          Daxpy(3*N_U, tau, rhs_vms_expl, B);

        // update rhs in explicit bubble VMS
        if ((TDatabase::ParamDB->DISCTYPE == VMS_RFB_EXPL)||
          (TDatabase::ParamDB->DISCTYPE == VMS_RFB_EXPL_COUPLED))
        {
          //  	  ApproximateTimeRFBSolutionQuad_Q2_NSE3D(coll, U1Array[mg_level-1], U2Array[mg_level-1],
          // 					U3Array[mg_level-1], PArray[mg_level-1],
          //                                         Coefficients[0], B);
          // 	  ApproximateTimeRFBSolutionQuadNSE3D(coll, U1Array[mg_level-1], U2Array[mg_level-1],
          // 					U3Array[mg_level-1], PArray[mg_level-1],
          //                                         Coefficients[0], B);

          ApproximateTimeRFB_coupled_SolutionQuad_Q2_NSE3D(coll, U1Array[mg_level-1], U2Array[mg_level-1],
            U3Array[mg_level-1], PArray[mg_level-1],
            Coefficients[0], B);
          OutPut(" RFB DONE"<<endl);
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
        // the matrices A_12, A_21, M_11, M_12, M_21, M_22, B1T, B2T
        //     stay unchanged during the complete solution process
        // the matrices A_11, A_22 are manipulated after their new
        //     assembling during the nonlinear iteration

        if ((m==1)&& (l==0) &&
          (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION == 1))
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
              break;
          }                                       // endswitch
        }                                         // end i
        // set current factor of steady state matrix
        gamma = -tau*theta2;

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

        // set Dirichlet values
        // RHSs[0] still available from assembling
        //         memcpy(B+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
        //         memcpy(B+N_Active+N_U, RHSs[1]+N_Active, (N_U-N_Active)*SizeOfDouble);
        //         memcpy(B+N_Active+2*N_U, RHSs[2]+N_Active, (N_U-N_Active)*SizeOfDouble);

        // copy Dirichlet values from rhs into sol
        //         memcpy(sol+N_Active, RHSs[0]+N_Active, (N_U-N_Active)*SizeOfDouble);
        //         memcpy(sol+N_Active+N_U, RHSs[1]+N_Active, (N_U-N_Active)*SizeOfDouble);
        //         memcpy(sol+N_Active+2*N_U, RHSs[2]+N_Active, (N_U-N_Active)*SizeOfDouble);

        //========================================================================
        // end assembling of rhs
        //========================================================================

        // extrapolate solution to get starting value at next time step
        // equidistant time steps assumed
        // current solution sol
        // solution of last time step sol_timestep_m1

        // save sol
        memcpy(oldsol, sol, SizeOfDouble*N_Unknowns);
        // sol := tau1 *sol - tau2 * sol_timestep_m1
        // at first time step: sol = sol_timestep_m1 -> result is sol
        tau2 = tau/oldtau;
        tau1 = 1 + tau2;
        for (k=0;k<3*N_U;k++)
          sol[k] = tau1*sol[k] - tau2*sol_timestep_m1[k];
        // save current solution
        memcpy(sol_timestep_m1, oldsol, SizeOfDouble*N_Unknowns);

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
              break;
          }                                       // endswitch
        }
        // set current factor of steady state matrix
        gamma = tau*theta1;

        //========================================================================
        // end assembling of system matrix
        //========================================================================

        real_time = TDatabase::TimeDB->CURRENTTIME * TDatabase::ParamDB->UREA_l_infty/ TDatabase::ParamDB->UREA_u_infty; 
        OutPut(endl << "CURRENT TIME: ");
        OutPut(TDatabase::TimeDB->CURRENTTIME);
        OutPut(" (real time: " << real_time << " s)" << endl);

        OutPut("CURRENT TIME: ");
        OutPut(TDatabase::TimeDB->CURRENTTIME);
        OutPut(" MEMORY: " << setw(10) << GetMemory()/(1048576.0));
        OutPut(" MB" << endl);

        //======================================================================
        // nonlinear loop
        //======================================================================
        INTERNAL_STEADY_STATE_MATRICES_OR_RHS=TDatabase::ParamDB->INTERNAL_STEADY_STATE_MATRICES_OR_RHS;
        N_LinIterCurr = 0;
        solver_time_curr = 0;

        // solve nonlinear equation
     if(!INTERNAL_STEADY_STATE_MATRICES_OR_RHS)
      {  
        for(j=0;j<Max_It;j++)
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
              break;
          }

          // compute defect
          Defect(sqmatrices,matrices,sol,B,defect);

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
          
          if ((((sqrt(residual)<=limit)||(j==Max_It-1)))
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
            //OutPut(TDatabase::TimeDB->CURRENTTIME << " total elapsed running time " <<    total_time << endl);
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
              default:
                OutPut("Unknown DISCTYPE " << TDatabase::ParamDB->DISCTYPE << endl);
                exit(1);
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
                N_RectMatrices = 0;
                N_Rhs = 0;
                N_FESpaces = 1;

                if (((TDatabase::ParamDB->LAPLACETYPE==1)
                  &&
                  ((CurrentDiscType == SMAGORINSKY) ||
                  (CurrentDiscType == GL00_AUX_PROBLEM) ||
                  (CurrentDiscType == VMS_PROJECTION)||
                  (CurrentDiscType == VMS_PROJECTION_EXPL)) ||
                  (CurrentDiscType == VMS_RFB_EXPL) ||
                  (CurrentDiscType == VMS_RFB_EXPL_COUPLED)) ||
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
                case GL00_AUX_PROBLEM:
                case VMS_PROJECTION_EXPL:
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

              if ((i==mg_level - 1)||(!TDatabase::ParamDB->VMS_COARSE_MG_SMAGO))
              {
                VMS_ProjectionUpdateMatrices(N_Uarray[i], USpaces[i]->GetActiveBound(),
                  ProjectionSpaces[i]->GetN_DegreesOfFreedom(),
                  SQMATRICES,MATRICES);
                OutPut("update done"<<endl);
              }
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
          OutPut("time for assembling " << t2-t1 << "s" << endl);
          // set current factor of steady state matrix
          gamma = tau*theta1;
        }                                         // endfor j / Max_It (solution of nonlinear equation)
     } 
        //======================================================================
        // end of nonlinear loop
        //======================================================================
        if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
        {
          IntoL20FEFunction3D(PArray[mg_level-1]->GetValues(),
            N_Parray[mg_level-1], PSpaces[mg_level-1]);
        }

        //======================================================================
        // computation of temperature
        //======================================================================
        CONC_TOL = TDatabase::ParamDB->UREA_CONC_TOL;
        CONC_MAXIT = TDatabase::ParamDB->UREA_CONC_MAXIT;
   
        // BEGIN LOOP
        // step 1: compute temp
    for(i_conc=0;i_conc<CONC_MAXIT;i_conc++) 
      // BEGIN LOOP
   {
       // OutPut("********* Computing temperature **********"<< endl);
        N_Unknowns_c = N_Unknowns_temp;
        N_Active_c = N_Active_temp;
        sol_c = sol_temp;
        oldsol_c = oldsol_temp;
        current_sol_c = sol_temp;
        itmethod_c = itmethod_temp;
        itmethod_sol_c = itmethod_sol_temp;
        RhsArray_c = RhsArray_temp;
        oldrhs_c = oldrhs_temp;
        rhs_c_complete = rhs_c_complete_A;

        Coeff_c = Coefficients[1];
        N_neum_to_diri_c = N_neum_to_diri_temp;
        neum_to_diri_c = neum_to_diri_temp;
        neum_to_diri_x = neum_to_diri_temp_x;
        neum_to_diri_y = neum_to_diri_temp_y;
        neum_to_diri_z = neum_to_diri_temp_z;
        oldrhs_fem_fct0_c = oldrhs_fem_fct0_temp;
        oldrhs_fem_fct1_c = oldrhs_fem_fct1_temp;
        lump_mass_c = lump_mass_temp;
        matrix_D_Entries_c = matrix_D_Entries_temp;
        tilde_u_c = tilde_u_temp;

        BoundaryConditions_Scalar[0] =  BoundCondition_temp;
        BoundValues_Scalar[0] = BoundValue_temp;

        for (i=0;i<mg_level;i++)
        {
          MatricesA_c[i] = MatricesA_temp[i];

          if(TDatabase::ParamDB->UREA_REACTION_DISC  == SDFEM)
          {
            MatricesK_c[i] = MatricesK_temp[i];
            if (TDatabase::ParamDB->SOLD_TYPE)
            {
              MatricesS_c[i] = MatricesS_temp[i];
            }
          }

          MatricesM_c[i] = MatricesM_temp[i];
          N_Array_c[i] = N_Array_temp[i];
          ConcentrationSpaces[i] = ConcentrationSpaces_temp[i];

          if ((TDatabase::ParamDB->UREA_REACTION_DISC  == SDFEM) &&
            (TDatabase::ParamDB->SOLD_TYPE))
          {
            SolArray[i] = SolArray_temp[i];
            SolArray_old[i] = SolArray_temp_old[i];
          }
          if ((TDatabase::ParamDB->UREA_REACTION_DISC == GALERKIN) &&
            (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
          {
            MatricesK_c[i] = MatricesK_temp[i];
          }
        }
        OutPut("********* Computing temperature **********"<< endl);
        // rhs
        // working array is B_c
        rhs_c = RhsArray_c[mg_level-1];
        // working array for rhs is B, initialize B
        memset(B_c, 0, N_Unknowns_c*SizeOfDouble);
        // old rhs multiplied with current subtime step and theta3 on B
        // old rhs not needed any longer now
        Daxpy(N_Active_c, tau*TDatabase::TimeDB->THETA3, oldrhs_c, B_c);

        // set parameter for compute term with new rhs
        N_Rhs = 1;
        RHSs[0] = RhsArray_c[mg_level-1];
        memset(RHSs[0],0,N_Array_c[mg_level-1]*SizeOfDouble);

       /* ferhs[0] = ConcentrationSpaces[mg_level-1];
        N_FESpaces = 2;
        fesp[0] = ConcentrationSpaces[mg_level-1];
        fesp[1] = USpaces[mg_level-1];
        fefct[0] = U1Array[mg_level-1];
        fefct[1] = U2Array[mg_level-1];
        fefct[2] = U3Array[mg_level-1];

        aux =  new TAuxParam3D(TimeCDParamsUreaN_FESpaces,
          TimeCDParamsUreaN_Fct,
          TimeCDParamsUreaN_ParamFct,
          TimeCDParamsUreaN_FEValues,
          fesp+1, fefct,
          TimeCDParamsUreaFct,
          TimeCDParamsUreaFEFctIndex,
          TimeCDParamsUreaFEMultiIndex,
          TimeCDParamsUreaN_Params,
          TimeCDParamsUreaBeginParam); */
        ferhs[0] = ConcentrationSpaces[mg_level-1];
        N_FESpaces = 4;
        //for conc
        fesp_conc[0] = ConcentrationSpaces[mg_level-1];
        //for the velocity
        fesp_conc[1] = USpaces[mg_level-1];
        //for conc
        fesp_conc[2] = ConcentrationSpaces_temp[mg_level-1];
        fesp_conc[3] = IntegralSpaces_conc[mg_level-1];

        fefct_conc[0] = U1Array[mg_level-1];
        fefct_conc[1] = U2Array[mg_level-1];
        fefct_conc[2] = U3Array[mg_level-1];
        fefct_conc[3] = SolArray_conc[mg_level-1];
        fefct_conc[4] = SolArray_temp[mg_level-1];
        fefct_conc[5] = IntegralSpaces_conc_fct[mg_level-1];

        aux =  new TAuxParam3D(TimeCDParamsUrea_tempN_FESpaces,
          TimeCDParamsUrea_tempN_Fct,
          TimeCDParamsUrea_tempN_ParamFct,
          TimeCDParamsUrea_tempN_FEValues,
          fesp_conc, fefct_conc,
          TimeCDParamsUrea_tempFct,
          TimeCDParamsUrea_tempFEFctIndex,
          TimeCDParamsUrea_tempFEMultiIndex,
          TimeCDParamsUrea_tempN_Params,
          TimeCDParamsUrea_tempBeginParam);

        switch(TDatabase::ParamDB->UREA_REACTION_DISC)
        {
          case UPWIND:
            OutPut("UREA_REACTION_DISC " << TDatabase::ParamDB->UREA_REACTION_DISC
              <<" not available !!!"<<endl);
            exit(1);
            //DiscreteForm = DiscreteFormRhs_Urea;
            break;
          case SDFEM:
            DiscreteForm = DiscreteFormRhs_SUPG_Urea;
            break;
          case GALERKIN:
            DiscreteForm = DiscreteFormRhs_Galerkin_Urea;
            break;
          default:
            OutPut("UREA_REACTION_DISC " << TDatabase::ParamDB->UREA_REACTION_DISC
              <<" not available !!!"<<endl);
            exit(4711);
        }

        if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
          (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
          BoundValues_Scalar[0] = BoundValue_FEM_FCT;

        Assemble3D(N_FESpaces, fesp_conc,
          0, NULL,
          0, NULL,
          N_Rhs, RHSs, ferhs,
          DiscreteForm,
          BoundaryConditions_Scalar,
          BoundValues_Scalar,
          aux);
        delete aux;

        if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
          (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
          BoundValues_Scalar[0] = BoundValue_temp;

        if ((TDatabase::ParamDB->UREA_REACTION_DISC == GALERKIN) &&
          (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
        {
          memcpy(oldrhs_fem_fct1_c, rhs_c, N_Unknowns_c*SizeOfDouble);
        }

        // add rhs from current sub time step to rhs array B_c
        Daxpy(N_Active_c, tau*TDatabase::TimeDB->THETA4, rhs_c, B_c);
        // save rhs for next time step
        memcpy(oldrhs_c, rhs_c, N_Unknowns_c * SizeOfDouble);

        // assembling of matrices on all levels
        for(i=low_scalar;i<mg_level;i++)
        {
          // set parameters
          N_FESpaces = 2;
          fesp[0] = ConcentrationSpaces[i];
          fesp[1] = USpaces[i];
          fefct[0] = U1Array[i];
          fefct[1] = U2Array[i];
          fefct[2] = U3Array[i];

          aux =  new TAuxParam3D(TimeCDParamsUreaN_FESpaces,
            TimeCDParamsUreaN_Fct,
            TimeCDParamsUreaN_ParamFct,
            TimeCDParamsUreaN_FEValues,
            fesp+1, fefct,
            TimeCDParamsUreaFct,
            TimeCDParamsUreaFEFctIndex,
            TimeCDParamsUreaFEMultiIndex,
            TimeCDParamsUreaN_Params,
            TimeCDParamsUreaBeginParam);

          N_SquareMatrices = 1;
          SQMATRICES[0] = MatricesA_c[i];
          SQMATRICES[0]->Reset();
          switch(TDatabase::ParamDB->UREA_REACTION_DISC)
          {
            case UPWIND:
              exit(1);
              break;
            case SDFEM:
              DiscreteForm = DiscreteFormMatricesA_SUPG_Urea;
              SQMATRICES[1] = MatricesK_c[i];
              SQMATRICES[1]->Reset();
              N_SquareMatrices = 2;
              if (TDatabase::ParamDB->SOLD_TYPE)
              {
                SQMATRICES[2] = MatricesS_c[i];
                SQMATRICES[2]->Reset();
                N_SquareMatrices = 3;
              }
              break;
            case GALERKIN:
              DiscreteForm = DiscreteFormMatricesA_Galerkin_Urea;
              break;
          }

          Assemble3D(N_FESpaces, fesp,
            N_SquareMatrices, SQMATRICES,
            0, NULL,
            0, NULL, NULL,
            DiscreteForm,
            BoundaryConditions_Scalar,
            BoundValues_Scalar,
            aux);

          //if(TDatabase::ParamDB->UREA_REACTION_DISC  == UPWIND)
          //{
          // UpwindForConvDiff(SQMATRICES[0],RHSs[0],
          //  fesp[0],DiscreteForm,
          //  U1Array[i], U2Array[i], 1);
          //}
          delete aux;
        }                                         // endfor i

        if (very_first_time==1)
        {
          very_first_time=0;
          l--;
          continue;
        }

        //======================================================================
        // manipulation of matrices due to current time discretization
        // the stiffness matrix is stored on M
        //======================================================================

        oldtau = tau;

        if (!((TDatabase::ParamDB->UREA_REACTION_DISC == GALERKIN) &&
          (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)))
        {

          // update rhs by Laplacian and convective term from previous
          // time step
          // scaled by current sub time step length and theta2
          // currently : M := M + gamma A
          // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2) A

          // complete rhs with data from the old time step
          // including mass term
          // no contribution from SOLD method

          for(i=low_scalar;i<mg_level;i++)
          {
            MatAdd(MatricesM_c[i], MatricesA_c[i],-gamma_c - tau*TDatabase::TimeDB->THETA2);
          }

          // set current factor of steady state matrix
          gamma_c = -tau*TDatabase::TimeDB->THETA2;

          MatVectActive(MatricesM_c[mg_level-1], oldsol_c, defect_c);
          Daxpy(N_Active_c, 1, defect_c, B_c);

          // contribution of SUPG term from time derivative
          if(TDatabase::ParamDB->UREA_REACTION_DISC  == SDFEM)
          {
            MatVectActive(MatricesK_c[mg_level-1], oldsol_c, defect_c);
            Daxpy(N_Active_c, 1, defect_c, B_c);
          }

          // set Dirichlet values
          memcpy(B_c+N_Active_c, RHSs[0]+N_Active_c, (N_Unknowns_c-N_Active_c)*SizeOfDouble);

          // save rhs B_c for the following iterations
          memcpy(rhs_c_complete, B_c, N_Unknowns_c*SizeOfDouble);

          // copy Dirichlet values from rhs into sol
          memcpy(sol_c+N_Active_c, B_c + N_Active_c, (N_Unknowns_c-N_Active_c)*SizeOfDouble);
          // M = M + (-gamma + tau*TDatabase::TimeDB->THETA1) A
          for(i=low_scalar;i<mg_level;i++)
          {
            MatAdd(MatricesM_c[i], MatricesA_c[i],-gamma_c + tau*TDatabase::TimeDB->THETA1);
          }
          // set current factor of steady state matrix
          gamma_c = tau*TDatabase::TimeDB->THETA1;
          // contribution of SUPG term from time derivative
          if(TDatabase::ParamDB->UREA_REACTION_DISC  == SDFEM)
          {
            for(i=low_scalar;i<mg_level;i++)
            {
              MatAdd(MatricesM_c[i], MatricesK_c[i], 1);
              // contribution from SOLD method only to left hand side
              if (TDatabase::ParamDB->SOLD_TYPE)
              {
                MatAdd(MatricesM_c[i], MatricesS_c[i],1.0);
              }
            }
          }
        }
        else
        {
          // FEM-FCT
          only_first_time = 1;
          // if(j>0)
          // only_first_time = 0;

          FEM_FCT_ForConvDiff(MatricesK_c[mg_level-1], MatricesA_c[mg_level-1],
            N_Unknowns_c, N_Active_c,
            lump_mass_c, matrix_D_Entries_c,
            sol_c, oldsol_c,
            B_c, RhsArray_c[mg_level-1], oldrhs_fem_fct0_c, tilde_u_c,
            N_neum_to_diri_c, neum_to_diri_c,
            neum_to_diri_x, neum_to_diri_y,  neum_to_diri_z,
            only_first_time,
            BoundValues_Scalar[0],NULL);

          SQMATRICES[0] = MatricesM_c[mg_level-1];
          // only first iteration
          // if (j==0)
          // {
          MatricesM_c[mg_level-1]->Reset();
          // system matrix for FEM-FCT   M_lump + theta1*tau*A
          // A = Galerkin + D

          FEM_FCT_SystemMatrix(MatricesM_c[mg_level-1], MatricesA_c[mg_level-1],
            lump_mass_c, N_Unknowns_c);
          // }

          // set Diriclet nodes
          if (N_neum_to_diri_c)
          {
            SetDirichletNodesFromNeumannNodes(SQMATRICES, B_c, sol_c,
              N_neum_to_diri_c, neum_to_diri_c,
              neum_to_diri_x, neum_to_diri_y, neum_to_diri_z,
              BoundValues_Scalar[0]);
          }
        }

        //======================================================================
        // solution of linear system
        //======================================================================

        memset(defect_c, 0, N_Unknowns_c*SizeOfDouble);
        SQMATRICES[0] = MatricesM_c[mg_level-1];
        // compute defect
        DefectScalar(sqmatrices,NULL,sol_c,B_c,defect_c);
        residual =  Ddot(N_Unknowns_c, defect_c, defect_c);
        residual_conc = residual*residual;
        OutPut("initial residual ");
        OutPut(setw(14) << residual << endl);
        if ((isnan(residual))||(isinf(residual)))
          {
            OutPut("Iteration diverged !!!" << endl);
            exit(4711);
          }
        //======================================================================
        // solve linear system
        //======================================================================
        switch(solver_type_reaction)
        {
          case AMG:
            TDatabase::ParamDB->SC_VERBOSE=0;
            t1 = GetTime();
            Solver(SQMATRICES[0], B, sol_c);
            t2 = GetTime();
            solver_time_curr = t2-t1;
            solver_time += solver_time_curr;
            break;

          case DIRECT:
            t1 = GetTime();
            DirectSolver(SQMATRICES[0], B_c, sol_c);
            t2 = GetTime();
            solver_time_curr = t2-t1;
            solver_time += solver_time_curr;
            break;

          case GMG:
            t1 = GetTime();
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
            {
              memcpy(itmethod_sol_c, sol_c, N_Unknowns_c*SizeOfDouble);
              memcpy(itmethod_rhs_c, B_c, N_Unknowns_c*SizeOfDouble);
            }

            N_LinIterCurrIte = itmethod_c->Iterate(sqmatrices,NULL,itmethod_sol_c,itmethod_rhs_c);
            N_LinIter+=N_LinIterCurrIte;
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
            {
              memcpy(sol_c, itmethod_sol_c, N_Unknowns_c*SizeOfDouble);
              memcpy(B_c, itmethod_rhs_c, N_Unknowns_c*SizeOfDouble);
            }
            t2 = GetTime();
            solver_time_curr += t2-t1;
            solver_time += solver_time_curr;
            break;
        }                                         // endswitch SOLVER_TYPE

        //======================================================================
        // end solve linear system
        //======================================================================

        // restore mass matrices by subtracting the K-matrices
        if(TDatabase::ParamDB->UREA_REACTION_DISC  == SDFEM)
        {
          for(i=low_scalar;i<mg_level;i++)
          {
            MatAdd(MatricesM_c[i], MatricesK_c[i], -1);
          }

        }
        // restore mass matrices by subtracting the A-matrices
        for(i=low_scalar;i<mg_level;i++)
        {
          MatAdd(MatricesM_c[i], MatricesA_c[i], -gamma_c);
        }
        // set current factor of steady state matrix
        gamma_c = 0;

        if ((TDatabase::ParamDB->UREA_REACTION_DISC == GALERKIN) &&
          (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
        {
          memcpy(oldrhs_fem_fct0_temp, oldrhs_fem_fct1_temp, N_Unknowns_temp*SizeOfDouble);
        }
        MG_temp->RestrictToAllGrids();
        //  steady = 0;
	//steady = 1;
//	} 
// matches if (!steady), line 2826
        // *********************************************************
        // temperature computed, compute now concentration
        // *********************************************************
  //       CONC_TOL = TDatabase::ParamDB->UREA_CONC_TOL;
   //      CONC_MAXIT = TDatabase::ParamDB->UREA_CONC_MAXIT;
   
        // BEGIN LOOP
        // step 1: compute conc
   // for(i_conc=0;i_conc<CONC_MAXIT;i_conc++) 
//// BEGIN LOOP
   //{
        N_Unknowns_c = N_Unknowns_conc;
        N_Active_c = N_Active_conc;
        sol_c = sol_conc;
        oldsol_c = oldsol_conc;
        current_sol_c = sol_conc;
        itmethod_c = itmethod_conc;
        itmethod_sol_c = itmethod_sol_conc;
        RhsArray_c = RhsArray_conc;
        oldrhs_c = oldrhs_conc;

        Coeff_c = Coefficients[2];
        N_neum_to_diri_c = N_neum_to_diri_conc;
        neum_to_diri_c = neum_to_diri_conc;
        neum_to_diri_x = neum_to_diri_conc_x;
        neum_to_diri_y = neum_to_diri_conc_y;
        neum_to_diri_z = neum_to_diri_conc_z;
        oldrhs_fem_fct0_c = oldrhs_fem_fct0_conc;
        oldrhs_fem_fct1_c = oldrhs_fem_fct1_conc;
        lump_mass_c = lump_mass_conc;
        matrix_D_Entries_c = matrix_D_Entries_conc;
        tilde_u_c = tilde_u_conc;

        BoundaryConditions_Scalar[0] =  BoundCondition_conc;
        BoundValues_Scalar[0] = BoundValue_conc;

        for (i=low_scalar;i<mg_level;i++)
        {
          MatricesA_c[i] = MatricesA_conc[i];
          if(TDatabase::ParamDB->UREA_REACTION_DISC  == SDFEM)
          {
            MatricesK_c[i] = MatricesK_conc[i];
            if (TDatabase::ParamDB->SOLD_TYPE)
              MatricesS_c[i] = MatricesS_conc[i];
          }
          MatricesM_c[i] = MatricesM_conc[i];
          N_Array_c[i] = N_Array_conc[i];
          ConcentrationSpaces[i] = ConcentrationSpaces_conc[i];
          if ((TDatabase::ParamDB->UREA_REACTION_DISC == GALERKIN) &&
            (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
            MatricesK_c[i] = MatricesK_conc[i];
        }
        OutPut("******** Computing concentration  ********"<< endl);

        gamma_c = 0;
        // assemble rhs only on finest level
        rhs_c = RhsArray_c[mg_level-1];
        // working array for rhs is B, initialize B
        memset(B_c, 0, N_Unknowns_c*SizeOfDouble);
        // old rhs multiplied with current subtime step and theta3 on B
        Daxpy(N_Active_c, tau*TDatabase::TimeDB->THETA3, oldrhs_c, B_c);

        // set parameter for compute term with new rhs
        N_Rhs = 1;
        RHSs[0] = RhsArray_c[mg_level-1];
        memset(RHSs[0],0,N_Array_c[mg_level-1]*SizeOfDouble);

        ferhs[0] = ConcentrationSpaces[mg_level-1];
        N_FESpaces = 4;
        //for conc
        fesp_conc[0] = ConcentrationSpaces[mg_level-1];
        //for the velocity
        fesp_conc[1] = USpaces[mg_level-1];
        //for conc
        fesp_conc[2] = ConcentrationSpaces_temp[mg_level-1];
        fesp_conc[3] = IntegralSpaces_conc[mg_level-1];

        fefct_conc[0] = U1Array[mg_level-1];
        fefct_conc[1] = U2Array[mg_level-1];
        fefct_conc[2] = U3Array[mg_level-1];
        fefct_conc[3] = SolArray_conc[mg_level-1];
        fefct_conc[4] = SolArray_temp[mg_level-1];
        fefct_conc[5] = IntegralSpaces_conc_fct[mg_level-1];

        aux =  new TAuxParam3D(TimeCDParamsUrea_concN_FESpaces,
          TimeCDParamsUrea_concN_Fct,
          TimeCDParamsUrea_concN_ParamFct,
          TimeCDParamsUrea_concN_FEValues,
          fesp_conc, fefct_conc,
          TimeCDParamsUrea_concFct,
          TimeCDParamsUrea_concFEFctIndex,
          TimeCDParamsUrea_concFEMultiIndex,
          TimeCDParamsUrea_concN_Params,
          TimeCDParamsUrea_concBeginParam);

        switch(TDatabase::ParamDB->UREA_REACTION_DISC)
        {
          case UPWIND:
            DiscreteForm = DiscreteFormRhs_Urea_conc;
            break;
          case SDFEM:
            DiscreteForm = DiscreteFormRhs_SUPG_Urea_conc;
            break;
          case GALERKIN:
            DiscreteForm = DiscreteFormRhs_Galerkin_Urea_conc;
            break;
        }

        if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
          (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
          BoundValues_Scalar[0] = BoundValue_FEM_FCT;

        Assemble3D(N_FESpaces, fesp_conc,
          0, NULL,
          0, NULL,
          N_Rhs, RHSs, ferhs,
          DiscreteForm,
          BoundaryConditions_Scalar,
          BoundValues_Scalar,
          aux);

        delete aux;

        if ((TDatabase::ParamDB->UREA_REACTION_DISC == GALERKIN) &&
          (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
        {
          BoundValues_Scalar[0] = BoundValue_conc;
          memcpy(oldrhs_fem_fct1_c, rhs_c, N_Unknowns_c*SizeOfDouble);
        }

        // add rhs from current sub time step to rhs array B_c
        Daxpy(N_Active_c, tau*TDatabase::TimeDB->THETA4, rhs_c, B_c);
        // save rhs for next time step
        memcpy(oldrhs_c, rhs_c, N_Unknowns_c * SizeOfDouble);

        // assembling of A and rhs
        // for SDFEM: in addition stabilisation matrix K
        for(i=low_scalar;i<mg_level;i++)
        {
          N_FESpaces = 4;
          //for conc
          fesp_conc[0] = ConcentrationSpaces[i];
          //for the velocity
          fesp_conc[1] = USpaces[mg_level-1];
          //for conc
          fesp_conc[2] = ConcentrationSpaces_temp[i];
          fesp_conc[3] = IntegralSpaces_conc[i];

          fefct_conc[0] = U1Array[i];
          fefct_conc[1] = U2Array[i];
          fefct_conc[2] = U3Array[i];
          fefct_conc[3] = SolArray_conc[i];
          fefct_conc[4] = SolArray_temp[i];
          fefct_conc[5] = IntegralSpaces_conc_fct[i];

          aux =  new TAuxParam3D(TimeCDParamsUrea_concN_FESpaces,
            TimeCDParamsUrea_concN_Fct,
            TimeCDParamsUrea_concN_ParamFct,
            TimeCDParamsUrea_concN_FEValues,
            fesp_conc, fefct_conc,
            TimeCDParamsUrea_concFct,
            TimeCDParamsUrea_concFEFctIndex,
            TimeCDParamsUrea_concFEMultiIndex,
            TimeCDParamsUrea_concN_Params,
            TimeCDParamsUrea_concBeginParam);

          //======================================================================
          // assembling of mass matrix and rhs
          //======================================================================

          N_SquareMatrices = 1;
          SQMATRICES[0] = MatricesA_c[i];
          SQMATRICES[0]->Reset();
          switch(TDatabase::ParamDB->UREA_REACTION_DISC)
          {
            case UPWIND:
              DiscreteForm = DiscreteFormMatricesA_Urea_conc;
              SQMATRICES[1] =  MatricesA_c[i];
              break;
            case SDFEM:
              DiscreteForm = DiscreteFormMatricesA_SUPG_Urea_conc;
              SQMATRICES[1] = MatricesK_c[i];
              SQMATRICES[1]->Reset();
              N_SquareMatrices = 2;
              if (TDatabase::ParamDB->SOLD_TYPE)
              {
                SQMATRICES[2] = MatricesS_c[i];
                SQMATRICES[2]->Reset();
                N_SquareMatrices = 3;
              }
              break;
            case GALERKIN:
              DiscreteForm = DiscreteFormMatricesA_Galerkin_Urea_conc;
              break;
          }

          Assemble3D(N_FESpaces, fesp_conc,
            N_SquareMatrices, SQMATRICES,
            0, NULL,
            0, NULL, NULL,
            DiscreteForm,
            BoundaryConditions_Scalar,
            BoundValues_Scalar,
            aux);

          //if(TDatabase::ParamDB->UREA_REACTION_DISC  == UPWIND)
          //{
          //  UpwindForConvDiff(SQMATRICES[0],RHSs[0],
          //    fesp_conc[0],DiscreteForm,
          //    U1Array[i], U2Array[i], 1);
          //}
          delete aux;

        }                                         // endfor i

        if (very_first_time==1)
        {
          very_first_time=0;
          l--;
          continue;
        }

        //======================================================================
        // manipulation of matrices due to current time discretization
        // the stiffness matrix is stored on M
        //======================================================================

        oldtau = tau;

        if (!((TDatabase::ParamDB->UREA_REACTION_DISC == GALERKIN) &&
          (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)))
        {
          // update rhs by Laplacian and convective term from previous
          // time step
          // scaled by current sub time step length and theta2
          // currently : M := M + gamma A
          // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2) A
          for(i=low_scalar;i<mg_level;i++)
          {
            MatAdd(MatricesM_c[i], MatricesA_c[i], -gamma_c - tau*TDatabase::TimeDB->THETA2);
          }
          // set current factor of steady state matrix
          gamma_c = -tau*TDatabase::TimeDB->THETA2;

          MatVectActive(MatricesM_c[mg_level-1], oldsol_c, defect_c);
          Daxpy(N_Active_c, 1, defect_c, B_c);
          // contribution of SUPG term from time derivative
          if(TDatabase::ParamDB->UREA_REACTION_DISC  == SDFEM)
          {
            MatVectActive(MatricesK_c[mg_level-1], oldsol_c, defect_c);
            Daxpy(N_Active_c, 1, defect_c, B_c);
          }

          // set Dirichlet values
          // RHSs[0] still available from assembling
          memcpy(B_c+N_Active_c, RHSs[0]+N_Active_c, (N_Unknowns_c-N_Active_c)*SizeOfDouble);
          // copy Dirichlet values from rhs into sol
          memcpy(sol_c+N_Active_c, RHSs[0]+N_Active_c, (N_Unknowns_c-N_Active_c)*SizeOfDouble);

          // M = M + (-gamma + tau*TDatabase::TimeDB->THETA1) A
          for(i=low_scalar;i<mg_level;i++)
          {
            MatAdd(MatricesM_c[i], MatricesA_c[i], -gamma_c + tau*TDatabase::TimeDB->THETA1);
          }
          // set current factor of steady state matrix
          gamma_c = tau*TDatabase::TimeDB->THETA1;
          // contribution of SUPG term from time derivative
          if(TDatabase::ParamDB->UREA_REACTION_DISC  == SDFEM)
          {
            for(i=low_scalar;i<mg_level;i++)
            {
              MatAdd(MatricesM_c[i], MatricesK_c[i], 1);
              if (TDatabase::ParamDB->SOLD_TYPE)
              {
                MatAdd(MatricesM_c[i], MatricesS_c[i],1.0);
              }
            }
          }
        }
        else
        {
          // FEM-FCT
          only_first_time = 1;
          // if(j>0)
          // only_first_time = 0;

          FEM_FCT_ForConvDiff(MatricesK_c[mg_level-1], MatricesA_c[mg_level-1],
            N_Unknowns_c, N_Active_c,
            lump_mass_c, matrix_D_Entries_c,
            sol_c, oldsol_c,
            B_c, RhsArray_c[mg_level-1], oldrhs_fem_fct0_c, tilde_u_c,
            N_neum_to_diri_c, neum_to_diri_c,
            neum_to_diri_x, neum_to_diri_y, neum_to_diri_z,
            only_first_time,
            BoundValues_Scalar[0],NULL);

          SQMATRICES[0] = MatricesM_c[mg_level-1];
          // only first iteration
          // if (j==0)
          // {
          MatricesM_c[mg_level-1]->Reset();
          // system matrix for FEM-FCT   M_lump + theta1*tau*A
          // A = Galerkin + D
          FEM_FCT_SystemMatrix(MatricesM_c[mg_level-1], MatricesA_c[mg_level-1],
            lump_mass_c, N_Unknowns_c);
          // }
          // set Diriclet nodes
          if (N_neum_to_diri_c)
            SetDirichletNodesFromNeumannNodes(SQMATRICES, B_c, sol_c,
              N_neum_to_diri_c, neum_to_diri_c,
              neum_to_diri_x, neum_to_diri_y, neum_to_diri_z,
              BoundValues_Scalar[0]);
        }

        //======================================================================
        // solution of linear system
        //======================================================================

        memset(defect_c, 0, N_Unknowns_c*SizeOfDouble);
        SQMATRICES[0] = MatricesM_c[mg_level-1];

        // compute defect
        DefectScalar(sqmatrices,NULL,sol_c,B_c,defect_c);
	//OutPut("sol_c " << Ddot(N_Unknowns_c, sol_c, sol_c) << " rhs_c " << Ddot(N_Unknowns_c, B_c, B_c) << endl)
        residual_conc =  Ddot(N_Unknowns_c, defect_c, defect_c);
        residual_conc += residual*residual;

        //nonlin_resid += residual;
        OutPut("initial residual ");
        OutPut(setw(14) << residual_conc << endl);
 
        //======================================================================
        // solve linear system
        //======================================================================
// nonlinear system solved(temp, concentration and psd)
 if ((sqrt(residual_conc) < CONC_TOL)&&(i_conc>0))
{
// text ausgeben
         OutPut("nonlinear system after " << i_conc <<" iterations " << residual_conc << endl);
        // restore mass matrices by subtracting the K-matrices
        if(TDatabase::ParamDB->UREA_REACTION_DISC  == SDFEM)
        {
          for(i=low_scalar;i<mg_level;i++)
          {
            MatAdd(MatricesM_c[i], MatricesK_c[i], -1);
          }
        }
        // restore mass matrices by subtracting the A-matrices
        for(i=low_scalar;i<mg_level;i++)
        {
          MatAdd(MatricesM_c[i], MatricesA_c[i], -gamma_c);
        }
        // set current factor of steady state matrix
        gamma_c = 0;

        //OutPut(TDatabase::TimeDB->CURRENTTIME << " NONLIN: " << N_LinIter << endl);

        max_c[0] = 0;
        max_c[1] = -1;
        for (i=0;i<N_Unknowns_c;i++)
        {
          if (sol_conc[i] < 0)
          {
            if (-sol_conc[i] >  max_c[0])
              max_c[0] = -sol_conc[i];
            // cut off negative values
            if (!((TDatabase::ParamDB->UREA_REACTION_DISC == GALERKIN) &&
              (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)))
              sol_conc[i] = 0;
          }
          if (sol_conc[i] > max_c[1])
            max_c[1] = sol_conc[i];
        }
        OutPut(TDatabase::TimeDB->CURRENTTIME <<
          " -conc " << max_c[0] <<  " +conc " << max_c[1] << endl);
         if (max_c[1] == -1)
        {
          OutPut("indefinite values for conc computed " << endl);
          exit(4711);
        }
        
	//OutPut("sol_c " << Ddot(N_Unknowns_c, sol_c, sol_c) << endl)

        // store solutions of this time step
      //  memcpy(oldsol_temp,sol_temp,N_Unknowns_temp*SizeOfDouble);
        //memcpy(oldsol_conc,sol_conc,N_Unknowns_conc*SizeOfDouble);
        MG_conc->RestrictToAllGrids();

      break;
 }
 
   switch(solver_type_reaction)
        {
          case AMG:
            TDatabase::ParamDB->SC_VERBOSE=0;
            t1 = GetTime();
            Solver(SQMATRICES[0], B_c, sol_c);
            t2 = GetTime();
            solver_time_curr = t2-t1;
            solver_time += solver_time_curr;
            break;

          case DIRECT:
            TDatabase::ParamDB->SC_VERBOSE=0;
            t1 = GetTime();
            DirectSolver(SQMATRICES[0], B_c, sol_c);
            t2 = GetTime();
            solver_time_curr = t2-t1;
            solver_time += solver_time_curr;
            break;

          case GMG:
            t1 = GetTime();
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
            {
              memcpy(itmethod_sol_c, sol_c, N_Unknowns_c*SizeOfDouble);
              memcpy(itmethod_rhs_c, B_c, N_Unknowns_c*SizeOfDouble);
            }
            N_LinIterCurrIte = itmethod_c->Iterate(sqmatrices,NULL,itmethod_sol_c,itmethod_rhs_c);
            N_LinIter+=N_LinIterCurrIte;
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
            {
              memcpy(sol_c, itmethod_sol_c, N_Unknowns_c*SizeOfDouble);
              memcpy(B_c, itmethod_rhs_c, N_Unknowns_c*SizeOfDouble);
            }
            t2 = GetTime();
            solver_time_curr += t2-t1;
            solver_time += solver_time_curr;
            break;
        }                                         // endswitch SOLVER_TYPE

        //======================================================================
        // end solve linear system
        //======================================================================

        // restore mass matrices by subtracting the K-matrices
        if(TDatabase::ParamDB->UREA_REACTION_DISC  == SDFEM)
        {
          for(i=low_scalar;i<mg_level;i++)
          {
            MatAdd(MatricesM_c[i], MatricesK_c[i], -1);
          }
        }
        // restore mass matrices by subtracting the A-matrices
        for(i=low_scalar;i<mg_level;i++)
        {
          MatAdd(MatricesM_c[i], MatricesA_c[i], -gamma_c);
        }
        // set current factor of steady state matrix
        gamma_c = 0;

        OutPut(TDatabase::TimeDB->CURRENTTIME << " LIN ITE: " << N_LinIter << endl);

        max_c[0] = 0;
        max_c[1] = -1;
        for (i=0;i<N_Unknowns_c;i++)
        {
          if (sol_conc[i] < 0)
          {
            if (-sol_conc[i] >  max_c[0])
              max_c[0] = -sol_conc[i];
            // cut off negative values
            if (!((TDatabase::ParamDB->UREA_REACTION_DISC == GALERKIN) &&
              (TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)))
              sol_conc[i] = 0;
          }
          if (sol_conc[i] > max_c[1])
            max_c[1] = sol_conc[i];
        }
         if (max_c[1] == -1)
        {
          OutPut("indefinite values for conc computed " << endl);
          exit(4711);
        }
        
	//OutPut("sol_c " << Ddot(N_Unknowns_c, sol_c, sol_c) << endl)

        // store solutions of this time step
      //  memcpy(oldsol_temp,sol_temp,N_Unknowns_temp*SizeOfDouble);
        //memcpy(oldsol_conc,sol_conc,N_Unknowns_conc*SizeOfDouble);
        MG_conc->RestrictToAllGrids();

        // **************************************************
        // concentration is computed
        // **************************************************

        // **************************************************
        // computing population size distribution
        // **************************************************

        OutPut("******** Computing f  ********"<< endl);
	/* compute rhs */
        if (TDatabase::ParamDB->UREA_MODEL != 1)
        {
	  OutPut("===== Aggregation and breakage ======="<< endl);
          //with aggregation     
          PrepareAgglomerationBreakage(coll, u1, u2, u2, temp, 
			     Nx, Ny, Nz, Na,
		     x_coord, y_coord, z_coord,
		     a_layers_coord, sol_psd, rhs_psd);
          
           //aggregation_conv_urea(coll, x_coord, y_coord, z_coord, Nx, Ny, Nz, Na,  sol_psd, rhs_psd, a_layers_coord, temp);
        }
        else 
       //without aggregation
        memset(rhs_psd, 0, Nodes*SizeOfDouble);
        
        if (TDatabase::ParamDB->UREA_PB_DISC== UREA_FWE_FDM_UPWIND)
        {
          //first time step
	  if (( m==1 )&&(TDatabase::ParamDB->READ_DATA<3))
          {
            //initial condition, since explicit scheme
            memset(sol_psd, 0, Nodes*SizeOfDouble);
          }
          else
          {
            OutPut("===== Begin Forward Euler FD Upwind Method ======="<< endl);

            Urea_FWE_FDM_Upwind_4D(coll, u1, u2, u3, conc, temp, sol_psd, sol_psd_old, rhs_psd,
                 Nx, Ny, Nz, Na,
                 x_coord, y_coord, z_coord, a_coord,
                 x_min, x_max, y_min, y_max, z_min, z_max, a_min, a_max,
                 velo1, velo2, velo3, concent_C_array, Temp_array, correspond_3dgrid);

                OutPut("===== End Forward Euler FD Upwind Method ======="<< endl);
          }
        }


        if ((TDatabase::ParamDB->UREA_PB_DISC== UREA_BWE_FDM_UPWIND)&&(TDatabase::ParamDB->READ_DATA<3))
        {
          OutPut("===== Begin Backward Euler FD Upwind Method ======="<< endl);

          Urea_BWE_FDM_Upwind_4D(coll, u1, u2, u3, conc,temp, sol_psd, rhs_psd,
            correspond_3dgrid, Nx, Ny, Nz, Na,
            x_coord, y_coord, z_coord, a_coord, mat);

          OutPut("===== End Backward Euler FD Upwind Method ======="<< endl);
        }

        // save the actual array of sol_psd for later computations
        // save_f_old_in_txt_file( Nodes, sol_psd, file_name_f_old);

        if (TDatabase::ParamDB->UREA_PB_DISC== UREA_FEM_FCT)
        {
          OutPut("===== Begin Crank-Nicolson FEM-FCT Method ======="<< endl);
          if (first_psd_fem)
          {
              Build_4D_FEM_FCT_MassMatrix_Q1(coll, Nx, Ny, Nz, Na, x_coord, y_coord,
              z_coord, a_coord, index_test_ansatz, matM_cons,
              lump_mass_PSD);

            Urea_Compute_Neum_To_Diri_FEM_FCT(Nx, Ny, Nz, Na,
					      x_coord, y_coord, z_coord, a_coord,
					      N_neum_to_diri_psd,
					      neum_to_diri_psd,
					      neum_to_diri_psd_x,
					      neum_to_diri_psd_y,
					      neum_to_diri_psd_z,
					      neum_to_diri_psd_a);
	    OutPut("Neum_to_diri " << N_neum_to_diri_psd << endl);
          first_psd_fem = 0;
          }

	  Urea_Build_4D_FEM_FCT_Matrix_Q1(coll, u1, u2, u3, conc, temp,
					  sol_psd_help, sol_psd, rhs_psd,
					  lump_mass_PSD, matrix_D_Entries_PSD,
					  correspond_3dgrid, Nx, Ny, Nz, Na,
					  x_coord, y_coord, z_coord, a_coord, mat, matM_cons, matM,
					  index_test_ansatz, N_neum_to_diri_psd, neum_to_diri_psd,
					  neum_to_diri_psd_x, neum_to_diri_psd_y,
					  neum_to_diri_psd_z, neum_to_diri_psd_a);
	  
          OutPut("===== End Crank-Nicolson FEM-FCT Method ======="<< endl);
       
         Dcopy(Nodes, rhs_psd, rhs_psd_old);
        }



        // compute f at the outflow boundary
       //  Evaluate_f_at_outflow1(Nx, Ny, Nz, Na, x_coord, y_coord, z_coord,
                        //              a_layers_coord, sol_psd);

        // compute integral values which are needed in the next time step

       
	 Integral_For_Particle_Increase_Term(integral_space_conc, integral_space_conc_fct,
              Nx, Ny, Nz, Na,
              x_coord, y_coord, z_coord, a_layers_coord,
              concent_C_array, sol_psd);
      
  }      // END LOOP
       //memset(sol_psd_old, 0, Nodes*SizeOfDouble);
        // store solutions of this time step
        //average_q3 = new double[Na];
        int step[1];
	step[0]=0;
        memcpy(oldsol_temp,sol_temp,N_Unknowns_temp*SizeOfDouble);
        memcpy(oldsol_conc,sol_conc,N_Unknowns_conc*SizeOfDouble);
        memcpy(sol_psd_old,sol_psd,N_Unknowns_temp*SizeOfDouble);
        //memset(average_q3,0,Na*SizeOfDouble);
         //memset(step,0,Na*SizeOfInt);

       
        
 //      OutPut("psd bis hier" << Ddot((Nx+1)*(Ny+1)*(Nz+1)*(Na+1),sol_psd,sol_psd) << endl);
       
	//Calculate_PSD_on_node(Nx, Ny, Nz, Na,
	//		      x_coord, y_coord, z_coord,
	//	      a_layers_coord, sol_psd, 0.0, 0.5, 0.5);
        //Calculate_PSD_on_node(Nx, Ny, Nz, Na,
          //                    x_coord, y_coord, z_coord,
            //                 a_layers_coord, sol_psd, 49., 0.5, 0.5);
	//Calculate_PSD_on_node(Nx, Ny, Nz, Na,
	//		      x_coord, y_coord, z_coord,
	//		      a_layers_coord, sol_psd, 100., 0.5, 0.5);
	//Calculate_PSD_on_node(Nx, Ny, Nz, Na,
	//		      x_coord, y_coord, z_coord,
	//		      a_layers_coord, sol_psd, 200.0, 0.5, 0.5);
       Calculate_PSD_outflow(coll, Nx, Ny, Nz, Na,
			      x_coord, y_coord, z_coord,
			      a_layers_coord, sol_psd,step, 200.0);
	    
	//Calculate_Volume_Distribution(coll, Nx, Ny, Nz, Na, 
	//			       x_coord, y_coord, z_coord,
	//			       a_layers_coord, sol_psd);
        N3=(Nx+1)*(Ny+1)*(Nz+1);
        N4=N3*(Na+1);
        memset(sol_psd + (N4-N3),0,N3*SizeOfDouble); 
        t22 = GetTime();
        OutPut( "time for timestep: " << t22-t11 << "s"<< endl);

                 
      }                                           // endfor l (sub steps of fractional step theta)
    }                                             // end of the two disc schemes (methods)

    if (time_discs==2)
    {
      // compute difference of solutions
      for ( i=0 ; i<N_Unknowns ; i++ )
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
      errors[8] = sqrt(errors[0]*errors[0]+errors[2]*errors[2]+errors[4]*errors[4]+errors[6]*errors[6]);

      if (TDatabase::TimeDB->CURRENTTIME< end_time)
        ComputeNewTimeStep(errors[8]);

      // copy solution of fract. step scheme
      memcpy(sol,frac_step_sol,N_Unknowns*SizeOfDouble);

      delete aux;
    }   

                                          // adaptive time step control


    // **************************************************************************
    //
    // the solution in the current discrete time is computed
    //
    // **************************************************************************

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
      OutPut( "L2(u): " << sqrt(errors[0]*errors[0]+errors[2]*errors[2]+errors[4]*errors[4]));
      OutPut( "   H1-semi(u):  " << sqrt(errors[1]*errors[1]+errors[3]*errors[3]+errors[5]*errors[5])<<endl);

      // error in L^infty(0,t,L^2)
      if (sqrt(errors[0]*errors[0]+errors[2]*errors[2]+errors[4]*errors[4]) > l_infty_l_2)
      {
        l_infty_l_2 = sqrt(errors[0]*errors[0] + errors[2]*errors[2] + errors[4]*errors[4]);
        l_infty_l_2_time =  TDatabase::TimeDB->CURRENTTIME;
      }
      OutPut( l_infty_l_2_time <<  " l_infty(L2(u)) " << l_infty_l_2 << endl);

      // error in L^2(0,t,L^2)
      l_2_l_2u += (errors[0]*errors[0] + errors[2]*errors[2] + errors[4]*errors[4]
        +olderror_l_2_l_2u) * TDatabase::TimeDB->TIMESTEPLENGTH/2.0;

      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "L2(0,t,L2)(u) : " <<  sqrt(l_2_l_2u) << endl);

      olderror_l_2_l_2u = errors[0]*errors[0] + errors[2]*errors[2]+errors[4]*errors[4];

      //error in L^2(0,t,H^1)
      l_2_h_1u += (errors[1]*errors[1] + errors[3]*errors[3]+ errors[5]*errors[5]
        +olderror_l_2_h_1u) * TDatabase::TimeDB->TIMESTEPLENGTH/2.0;

      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "L2(0,t,H1-semi)(u) : " << sqrt(l_2_h_1u) << endl);

      olderror_l_2_h_1u = errors[1]*errors[1] + errors[3]*errors[3]+ errors[5]*errors[5];

      U1Array[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors);
      U2Array[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors+2);
      U3Array[mg_level-1]->GetErrors(ExactNull, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, USpaces+mg_level-1, errors+4);
      OutPut( TDatabase::TimeDB->CURRENTTIME << "  " );
      OutPut( "kinetic energy " << (errors[0]*errors[0]+errors[2]*errors[2]+errors[4]*errors[4])/2 << endl);

      PArray[mg_level-1]->GetErrors(ExactP, 4, TimeNSAllDerivatives,
        2, L2H1Errors, NULL, aux, 1, PSpaces+mg_level-1, errors);
      OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
      OutPut( "L2(p): " << errors[0]);
      OutPut( "   H1-semi(p): " << errors[1]<<endl);

      if  ((TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE != 4)
        &&(TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE != 100))
      {
        // subgrid dissipation
        UArray[mg_level-1]->GetDeformationTensorErrors(ExactU1, ExactU2, ExactU3,
          4, TimeNSAllDerivatives,
          1, SubGridDissipation,
          NULL, aux, 1, USpaces+mg_level-1, errors);
        OutPut( TDatabase::TimeDB->CURRENTTIME << " " );
        OutPut( "subgrid dissipation : " << errors[0] << endl);
      }
      else
        OutPut( "subgrid dissipation not implemented " << endl);

      delete aux;
    }                                             // endif MEASURE_ERRORS

    if ((TDatabase::ParamDB->WRITE_GRAPE)||(TDatabase::ParamDB->WRITE_GMV)||(TDatabase::ParamDB->WRITE_VTK)
      || (TDatabase::ParamDB->SAVE_DATA))
    {
      if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      {
        if (TDatabase::ParamDB->WRITE_GMV)
        {
          os.seekp(std::ios::beg);
          os << GmvBaseName << m << ".gmv" << ends;
          Output->WriteGMV(os.str().c_str());
        }
        if (TDatabase::ParamDB->WRITE_VTK)
        {
          os.seekp(std::ios::beg);
          os << VtkBaseName << m << ".vtk" << ends;
          Output->WriteVtk(os.str().c_str());
          //os.seekp(std::ios::beg);
          //os << VtkBaseName << "end.psd." << m << ".vtk" << ends;
          //write_vtk_file(Nx, Ny, Nz, Na, cut_coord, x_coord, y_coord, z_coord, a_coord, sol_psd,os.str().c_str());
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
          save_sol[1] = sol_temp;
          save_sol[2] = sol_conc;
          save_sol[3] = sol_psd;
          save_N_Unknowns[0] = N_Unknowns;
          save_N_Unknowns[1] = N_Unknowns_temp;
          save_N_Unknowns[2] = N_Unknowns_conc;
          save_N_Unknowns[3] = Nodes;
          SaveData(SaveDataFileName,4,save_sol,save_N_Unknowns);
        }
      }
    }
    comp_vort =0;
    OutPut(TDatabase::TimeDB->CURRENTTIME << " total elapsed running time " <<    GetTime()-total_time << endl);
 
  }                                               // while

  //======================================================================
  // end of time cycle
  //======================================================================

  if ((TDatabase::ParamDB->WRITE_GRAPE)||(TDatabase::ParamDB->WRITE_GMV)||(TDatabase::ParamDB->WRITE_VTK))
  {
    if (TDatabase::ParamDB->WRITE_GMV)
    {
      os.seekp(std::ios::beg);
      os << GmvBaseName << m << ".gmv" << ends;
      Output->WriteGMV(os.str().c_str());
    }
    if (TDatabase::ParamDB->WRITE_VTK)
    {
      os.seekp(std::ios::beg);
      os << VtkBaseName << "end."<< m << ".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      //os.seekp(std::ios::beg);
      //  os << VtkBaseName << "end.psd." << m << ".vtk" << ends;
      //write_vtk_file(Nx, Ny, Nz, Na, cut_coord, x_coord, y_coord, z_coord, a_coord, sol_psd,os.str().c_str());
    }
    if (TDatabase::ParamDB->WRITE_GRAPE)
    {
      os.seekp(std::ios::beg);
      os << GrapeBaseName << m << ".dat" << ends;
      Output->WriteGrape(os.str().c_str());
    }
  }

  t4 = GetTime();
  total_time += t4 - t3;
  OutPut("total running time: " << total_time << endl);

  CloseFiles();

  return 0;


}
// =======================================================================
// end of main programm
// =======================================================================
