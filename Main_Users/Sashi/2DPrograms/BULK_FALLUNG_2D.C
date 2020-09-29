// =======================================================================
//
// Purpose:     main program
//
// Author:      Volker John
//
// =======================================================================

// =======================================================================
//   Header Files
// =======================================================================
#include <Assemble2D.h>
#include <AuxParam2D.h>
#include <Bulk_2d3d.h>
#include <Collection.h>
#include <Convolution.h>
#include <Database.h>
#include <DirectSolver.h>
#include <DiscreteForm2D.h>
#include <Domain.h>
#include <FEDatabase2D.h>
#include <FEM_TVD_FCT.h>
#include <FESpace2D.h>
#include <FgmresIte.h>
#include <FixedPointIte.h>
#include <ItMethod.h>
#include <LinAlg.h>
#include <math.h>
#include <MGLevel2D.h>
#include <MultiGrid2D.h>
#include <MultiGridIte.h>
#include <MultiGridScaIte.h>
#include <NSE2D_ParamRout.h>
#include <NSE_MultiGrid.h>
#include <NSE_MGLevel1.h>
#include <NSE_MGLevel2.h>
#include <NSE_MGLevel3.h>
#include <NSE_MGLevel4.h>
#include <Output2D.h>
#include <QuadAffin.h>
#include <Solver.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <TCD2D.h>
#include <TNSE2D_ParamRout.h>
#include <Upwind.h>
#include <VMS.h>

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

#include "Examples/TNSE_2D/Bulk_Fallung_Driven_Cavity.h"
//#include "Examples/TNSE_2D/Bulk_Fallung_Driven_Cavity_Academic.h"

int main(int argc, char* argv[])
{
  //======================================================================
  // begin of the declaration of the variables
  //======================================================================
  // integer variables
  int i, j, k, l, m, n, ii, ll;
  int ansatz_order_reaction, BASELEVEL, comp_vort, CurrentDiscType;
  int FirstSolve;
  int last_digit_ite, last_sq, Len, level_down, LEVELS, low;
  int Max_It, Max_It_c, methods, mg_level, mg_type;
  int mom, N_mom, N_Unknowns_mom, N_Active_mom;
  int N, N_, N_Active, N_Active_c, N_Active_c_A, N_Active_c_B, N_Active_c_C, n_aux, N_Cells, N_CConv, N_Columns, n_dof;
  int N_Entries, N_FESpaces, N_L, N_LinIter, N_LinIterCurr, N_LinIterCurrIte, N_NonActive, Nodes, nonlin_ite;
  int N_P, N_P_low, N_P_sep, N_RectMatrices, N_Rhs, N_Rows, N_SquareMatrices, N_SubSteps;
  int N_U, N_UConv, N_U_low, number_pressep;
  int N_Unknowns, N_Unknowns_c, N_Unknowns_c_A, N_Unknowns_c_B, N_Unknowns_c_C, N_Unknowns_Integral_Space, N_Unknowns_low;
  int N_V, N_Vort, Nx, Ny, Nz, low_scalar;
  int ORDER, order, pde, pressure_space_code, ret;
  int slow_conv, solver_type_reaction, step_length, time_discs;
  int velocity_space_code, which, zerostart, only_first_time;

  // initialised integer variables
  int calculations = 1, fully_implicit = 0, mixing_layer = 0, N_GNU_images = 0;
  int N_neum_to_diri_c = 0, N_neum_to_diri_c_A = 0, N_neum_to_diri_c_B = 0, N_neum_to_diri_c_C = 0;
  int N_neum_to_diri_psd = 0;
  int N_Paramters = 1, N_Parameters = 2, pre_calculation = 1, pressure_separation = 0;
  int very_first_time=0, substance = 1, first_psd_fem = 1;
  int average_step[1];

  // integer pointers
  int *N_Array, *N_Array_c, *N_Array_c_A, *N_Array_c_B, *N_Array_c_C, *N_Uarray, *N_Parray;
  int *correspond_2dgrid, *col_ptr, *permutation, *row_ptr,*RowPtr, *N_Array_mom;
  int *neum_to_diri_c_A, *neum_to_diri_bdry_c_A, *neum_to_diri_c_B, *neum_to_diri_bdry_c_B;
  int *neum_to_diri_c_C, *neum_to_diri_bdry_c_C, *neum_to_diri_c, *neum_to_diri_bdry_c;
  int *neum_to_diri_psd;
  int *index_test_ansatz, *index_test_ansatz_diag;

  int **downwind;

  // double variables
  double average;
  double cd, cl, convergence_speed;
  double delta, delta0, delta1, delta_time_step, DiffL2, DiffH1, divergence, dp;
  double end_time;
  double firsterror, firsterrorl2;
  double gamma, gamma_c;
  double h, hmin, hmax;
  double impuls_residual;
  double lasterror, lasterrorl2, limit, limit_c, linredfac;
  double max, min;
  double negPower;
  double oldres, oldresidual, oldtau, omega;
  double p1, p2, p3, p4;
  double real_time, reatt_pt, RE_NR, res, res2, residual, residual_Cc;
  double solver_time, solver_time_curr, sum;
  double t, t1, t2, t3, t4, tau, tau1, tau2, tol, tolerance, total_time;
  double velo_l2, vort_zero, vort_zero_conv;
  double x, y, residual_mom;
  double theta1, theta2, theta3, theta4;

  // initialised double variables
  double l_2_h_1u = 0.0, l_infty_l_2 = 0.0, l_infty_l_2_time = -4711.0;
  double l_2_l_2Du = 0.0, l_2_l_2u = 0.0, nonlin_resid = 1e8;
  double olderror = 0.0, olderror_l_2_h_1u = 0.0, olderror_l_2_l_2u = 0.0;
  double average_median[1];
  double max_c[6];

  // fixing of the interval limits for the coordinates
  double x_min = 0, x_max = 1;
  double y_min = 0, y_max = 1;
  double z_min, z_max = 1;

  // double pointers
  double *app, *auxConv, *Auxitmethod_sol, *Auxitmethod_rhs;
  double *B, *B_c, *concent_C_array, *conv_vort, *current_B, *current_sol, *current_sol_c;
  double *defect, *defect_c, *div, *du_tensor;
  double *fesol, *frac_step_sol, *frac_step_sol_c, *f_time_space_approx_old, *f_time_space_approx;
  double *h1p, *h1u1, *h1u2;
  double *integral_val, *itmethod_rhs, *itmethod_rhs_c;
  double *itmethod_sol, *itmethod_sol_c, *itmethod_sol_c_A, *itmethod_sol_c_B, *itmethod_sol_c_C;
  double *l_inf, *l2p, *l2u1, *l2u2, *LESModelRhs, *nosep_p;
  double *lump_mass_PSD, *lump_mass_c_A, *lump_mass_c_B, *lump_mass_c_C, *lump_mass_c;
  double *oldrhs, *oldrhs_c, *oldrhs_c_A, *oldrhs_c_B, *oldrhs_c_C;
  double *oldsol, *oldsol_c, *oldsol_c_A, *oldsol_c_B, *oldsol_c_C;
  double *oldrhs_fem_fct0_c_A, *oldrhs_fem_fct0_c_B, *oldrhs_fem_fct0_c_C;
  double *oldrhs_fem_fct1_c_A, *oldrhs_fem_fct1_c_B, *oldrhs_fem_fct1_c_C;
  double *oldrhs_fem_fct0_c, *oldrhs_fem_fct1_c;
  double *rhs_c_complete_A, *rhs_c_complete_B, * rhs_c_complete;
  double *pressure_aux_array, *psi;
  double *rhs, *rhs_high, *rhs_low, *rhsPressSep, *rhs_c;
  double *sd, *separated_pressure_array, *separated_pressure_aux_array;
  double *startsol, *startsol_c, *sol, *sol_c, *sol_c_A, *sol_c_B, *sol_c_C, *soldiff;
  double *sol_low, *sol_timestep_m1;
  double *u_conv, *u_uConv;
  double *val, *velo1, *velo2, *vorticity;
  double *sol_mom, *oldsol_mom, *B_mom, *itmethod_sol_mom, *itmethod_rhs_mom, *rhs_mom;
  double *oldrhs_mom, *defect_mom, *rhs_single_mom;
  double *neum_to_diri_param_c_A, *neum_to_diri_param_c_B, *neum_to_diri_param_c_C, *neum_to_diri_param_c;
  double *matrix_D_Entries_PSD, *matrix_D_Entries_c_A, *matrix_D_Entries_c_B, *matrix_D_Entries_c_C, *matrix_D_Entries_c;
  double *tilde_u_c, *tilde_u_c_A, *tilde_u_c_B, *tilde_u_c_C;
  double *neum_to_diri_psd_x, *neum_to_diri_psd_y, *neum_to_diri_psd_z;

  // declaration and initialisation of the coordinate vectors
  double *x_coord, *y_coord, *z_coord;
  // declaration and initialisation of the vector with the differences between the z layers
  double *z_layers_coord;

  double **RhsArray, **RhsArray_c, **RhsArray_c_A, **RhsArray_c_B, **RhsArray_c_C;
  double **RhsArray_mom;

  // double pointer on an array
  double *RHSs[4];

  // double arrays
  double alpha[2], alpha_fine[2], errors[7], errors_mg[4], exactvect[3], Parameters[4];
  double reatt_point[3], residuals[10], vect[3];
  double *save_sol[5];
  int save_N_Unknowns[5];

  // variables of MooNMD classes
  BoundCondFunct2D *BoundaryConditions[2], *BoundaryConditions_Scalar[4];
  BoundCondFunct2D *BoundaryConditionsAuxProblem[3], *BoundaryConditionsPressureSeparation[1];

  BoundValueFunct2D *BoundValues[2], *BoundValues_Scalar[4], *BoundValuesAuxProblem[3];
  BoundValueFunct2D *BoundaryValuesPressureSeparation[1];

  CoeffFct2D *Coefficients[4], *Coeff_c;
  DefectProc *Defect, *DefectScalar;

  MatVecProc *MatVect, *MatVectScalar;

  TAuxParam2D *aux;

  TBaseCell *cell;

  TCollection *coll, *mortarcoll = NULL;

  TDatabase *Database = new TDatabase();

  TDiscreteForm2D *DiscreteForm;
  TDiscreteForm2D *DiscreteFormAuxProbPressSep;
  TDiscreteForm2D *DiscreteFormColetti;
  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormGL00Convolution;
  TDiscreteForm2D *DiscreteFormGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormGL00AuxProblemRHS;
  TDiscreteForm2D *DiscreteFormMatricesA_SUPG_Bulk;
  TDiscreteForm2D *DiscreteFormMatricesA_SUPG_Bulk_Cc;
  TDiscreteForm2D *DiscreteFormMatricesA_SUPG_Bulk_mom;
  TDiscreteForm2D *DiscreteFormRhs_SUPG_Bulk;
  TDiscreteForm2D *DiscreteFormRhs_SUPG_Bulk_Cc;
  TDiscreteForm2D *DiscreteFormRhs_SUPG_Bulk_mom;
  TDiscreteForm2D *DiscreteFormMatricesA_Bulk;
  TDiscreteForm2D *DiscreteFormMatricesA_Bulk_Cc;
  TDiscreteForm2D *DiscreteFormRhs_Bulk;
  TDiscreteForm2D *DiscreteFormRhs_Bulk_Cc;
  TDiscreteForm2D *DiscreteFormMatricesA_Galerkin_Bulk;
  TDiscreteForm2D *DiscreteFormMatricesA_Galerkin_Bulk_Cc;
  TDiscreteForm2D *DiscreteFormRhs_Galerkin_Bulk;
  TDiscreteForm2D *DiscreteFormRhs_Galerkin_Bulk_Cc;
  TDiscreteForm2D *DiscreteFormMatrixAuxProblemU;
  TDiscreteForm2D *DiscreteFormMatrixGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormMatrixMBulk;
  TDiscreteForm2D *DiscreteFormNLColetti;
  TDiscreteForm2D *DiscreteFormNLGalerkin;
  TDiscreteForm2D *DiscreteFormNLGL00AuxProblem;
  TDiscreteForm2D *DiscreteFormNLGL00Convolution;
  TDiscreteForm2D *DiscreteFormNLSDFEM;
  TDiscreteForm2D *DiscreteFormNLSmagorinsky;
  TDiscreteForm2D *DiscreteFormNLUpwind;
  TDiscreteForm2D *DiscreteFormNLVMS_Projection;
  TDiscreteForm2D *DiscreteFormNSRFBRhs;
  TDiscreteForm2D *DiscreteFormPressSep;
  TDiscreteForm2D *DiscreteFormRHS;
  TDiscreteForm2D *DiscreteFormRHSAuxProblemU;
  TDiscreteForm2D *DiscreteFormRHSColetti;
  TDiscreteForm2D *DiscreteFormRHSLESModel;
  TDiscreteForm2D *DiscreteFormRHSSmagorinskyExpl;
  TDiscreteForm2D *DiscreteFormSDFEM;
  TDiscreteForm2D *DiscreteFormSmagorinsky;
  TDiscreteForm2D *DiscreteFormUpwind;
  TDiscreteForm2D *DiscreteFormVMS_Projection;

  TDomain *Domain = new TDomain();

  TFEDatabase2D *FEDatabase = new TFEDatabase2D();

  TFEFunction2D *AuxPArray;
  TFEFunction2D  *c_A, *c_B, *c_C, *c_A_old, *c_B_old, *c_C_old;
  TFEFunction2D *Conv_Vort, *Divergence;
  TFEFunction2D *du1Conv, *du2Conv, *du3Conv;
  TFEFunction2D *GL00AuxProblemSol11, *GL00AuxProblemSol12, *GL00AuxProblemSol22;
  TFEFunction2D *integral_space_c_C_fct;
  TFEFunction2D *old_p, *p, *p_low;
  TFEFunction2D *separated_pressure_fe_funct, *separated_pressure_rhs_fe_funct;
  TFEFunction2D *soldiff_fe1, *soldiff_fe2, *StreamFct;
  TFEFunction2D *u1, *u2, *u1Conv, *u2Conv, *u1_low, *u2_low;
  TFEFunction2D *Vorticity;

  TFEFunction2D *fefct_c_C[9], *fefct_mom[4];
  TFEFunction2D *fefct[7];

  TFEFunction2D **IntegralSpaces_c_C_fct;
  TFEFunction2D **PArray;
  TFEFunction2D **SolArray_c_A, **SolArray_c_B, **SolArray_c_C, **SolArray_other, **SolArray_other_c;
  TFEFunction2D **SolArray_c_A_old, **SolArray_c_B_old, **SolArray_c_C_old;
  TFEFunction2D **SolArray_old, **SolArray_other_old, **SolArray;
  TFEFunction2D **U1Array, **U2Array,  **AuxFEFunctArray;
  TFEVectFunct2D **AuxFEVectFunctArray;

  TFESpace2D *concentration_space, *concentration_space_b;
  TFESpace2D *concentration_space_c;
  TFESpace2D *convolution_space;
  TFESpace2D *integral_space_c_C;
  TFESpace2D *old_u_space, *old_p_space;
  TFESpace2D *pressure_separation_space;
  TFESpace2D *pressure_space, *pressure_space_low;
  TFESpace2D *projection_space;
  TFESpace2D *streamfunction_space;
  TFESpace2D *velocity_space, *velocity_space_low;
  TFESpace2D *vorticity_space;
  TFESpace2D *mom_space;

  TFESpace2D *fesp_c_C[6], *fesp_mom[3];
  TFESpace2D *fesp[4], *ferhs[3];

  TFESpace2D **ConcentrationSpaces;
  TFESpace2D **ConcentrationSpaces_c_A;
  TFESpace2D **ConcentrationSpaces_c_B;
  TFESpace2D **ConcentrationSpaces_c_C;
  TFESpace2D **ConcentrationSpaces_other;
  TFESpace2D **duConvSpaces;
  TFESpace2D **IntegralSpaces_c_C;
  TFESpace2D **PsiSpaces;
  TFESpace2D **PSpaces;
  TFESpace2D **ProjectionSpaces;
  TFESpace2D **uConvSpaces;
  TFESpace2D **USpaces;
  TFESpace2D **VorticitySpaces;
  TFESpace2D **MOMSpaces;

  TFEVectFunct2D *u, *u_low, *old_u, *fe_mom, *fe_mom_old;
  TFEVectFunct2D **UArray, **SolArray_mom, **SolArray_mom_old;

  TItMethod *Auxprec, *Auxitmethod, *itmethod, *itmethod_c, *prec;
  TItMethod *itmethod_c_A, *prec_c_A, *itmethod_c_B, *prec_c_B;
  TItMethod *itmethod_c_C, *prec_c_C;
  TItMethod *itmethod_mom, *prec_mom;

  TMatrix2D *MATRICES[10];

  TMatrix **matrices = (TMatrix **)MATRICES;

  TMatrix2D *matrixB1, *matrixB2;
  TMatrix2D *matrixB1_low, *matrixB2_low;
  TMatrix2D *matrixB1T, *matrixB2T;
  TMatrix2D *matrixB1T_low, *matrixB2T_low;
  TMatrix2D *matrix_G11, *matrix_G22;
  TMatrix2D *matrix_tilde_G11, *matrix_tilde_G22;

  TMatrix2D **MatricesB1, **MatricesB2;
  TMatrix2D **MatricesB1_low, **MatricesB2_low;
  TMatrix2D **MatricesB1T_low, **MatricesB2T_low;
  TMatrix2D **MatricesB1T, **MatricesB2T;
  TMatrix2D **Matrices_G11, **Matrices_G22;
  TMatrix2D **Matrices_tilde_G11, **Matrices_tilde_G22;

  TMGLevel2D *MGLevel_c_A, *MGLevel_c_B, *MGLevel_c_C, *MGLevel_mom;

  TMultiGrid2D *MG_c_A, *MG_c_B, *MG_c_C, *MG_mom;

  TNSE_MGLevel *MGLevel, *MGLevel_low;

  TNSE_MultiGrid *MG;

  TOutput2D *Output;

  TSquareMatrix2D *SQMATRICES[8];
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;

  // matrix structure
  TSquareMatrix2D *mat;
  TSquareMatrix2D *matM, *matM_cons;
  TSquareMatrix2D *sqmatrixA;
  TSquareMatrix2D *sqmatrixA_low;
  TSquareMatrix2D *sqmatrixA11, *sqmatrixA12;
  TSquareMatrix2D *sqmatrixA21, *sqmatrixA22;
  TSquareMatrix2D *sqmatrixA11_low, *sqmatrixA12_low;
  TSquareMatrix2D *sqmatrixA21_low, *sqmatrixA22_low;
  TSquareMatrix2D *sqmatrixK;
  TSquareMatrix2D *sqmatrixK_c_A, *sqmatrixK_c_B, *sqmatrixK_c_C;
  TSquareMatrix2D *sqmatrixS;
  TSquareMatrix2D *sqmatrixS_c_A, *sqmatrixS_c_B, *sqmatrixS_c_C;
  TSquareMatrix2D *sqmatrixL, *sqmatrixM;
  TSquareMatrix2D *sqmatrixM_c_A, *sqmatrixM_c_B, *sqmatrixM_c_C;
  TSquareMatrix2D *sqmatrixM_mom, *sqmatrixK_mom, *sqmatrixS_mom;
  TSquareMatrix2D *sqmatrixM11, *sqmatrixM12;
  TSquareMatrix2D *sqmatrixM21, *sqmatrixM22;
  TSquareMatrix2D *sqmatrixPressSep;

  TSquareMatrix2D **MatricesA, **MatricesA_c;
  TSquareMatrix2D **MatricesA_low;
  TSquareMatrix2D **MatricesA11, **MatricesA12;
  TSquareMatrix2D **MatricesA21, **MatricesA22;
  TSquareMatrix2D **MatricesA11_low, **MatricesA12_low;
  TSquareMatrix2D **MatricesA21_low, **MatricesA22_low;
  TSquareMatrix2D **MatricesA_c_A, **MatricesA_c_B, **MatricesA_c_C;
  TSquareMatrix2D **MatricesK, **MatricesK_c;
  TSquareMatrix2D **MatricesK_c_A, **MatricesK_c_B, **MatricesK_c_C;
  TSquareMatrix2D **MatricesS, **MatricesS_c;
  TSquareMatrix2D **MatricesS_c_A, **MatricesS_c_B, **MatricesS_c_C;
  TSquareMatrix2D **MatricesL, **MatricesM, **MatricesM_c;
  TSquareMatrix2D **MatricesM_c_A, **MatricesM_c_B, **MatricesM_c_C;
  TSquareMatrix2D **MatricesM11, **MatricesM12;
  TSquareMatrix2D **MatricesM21, **MatricesM22;
  TSquareMatrix2D **MatricesA_mom, **MatricesM_mom, **MatricesK_mom;
  TSquareMatrix2D **MatricesS_mom;

  // matrix structure
  TSquareStructure2D *matrix_structure;
  TSquareStructure2D *sqstructureA, *sqstructureA_low, *sqstructureC;
  TSquareStructure2D *sqstructure_c_A, *sqstructure_c_B, *sqstructure_c_C;
  TSquareStructure2D *sqstructureL, *sqstructurePressSep, *sqstructure_mom;

  TStructure2D *structureB, *structureBT;
  TStructure2D *structureB_low, *structureBT_low;
  TStructure2D *structure_tilde_G, *structure_G;

  // char variables

  // strings
  char ReadinDat[] = "readin.dat";
  char NameString[] = "name";
  char DString[] = "d";
  char UString[] = "u";
  char PString[] = "p";
  char c_A_String[] = "c_A";
  char c_B_String[] = "c_B";
  char c_C_String[] = "c_C";
  char PsiString[] = "psi";
  char PsepString[] = "psep";
  char MassMatrix[] = "Mass matrix";
  char Mass[] = "Mass";
  char UConvString[] = "u_conv";
  char UConfString[] = "uconf";
  char AuxProbString[] = "AuxProblem";
  char VorticityString[] = "vorticity";
  char ConvVortString[] = "conv_vort";
  char DivergenceString[] = "divergence";

  // char pointers
  char *PsBaseName, *GrapeBaseName, *GnuBaseName, *ReadGrapeBaseName;
  char *GmvBaseName, *VtkBaseName, *MatlabBaseName;
  char *SaveDataFileName,  *ReadDataFileName;
  char *PRM, *GEO, *MAP;

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
  TDatabase::ParamDB->SDFEM_TYPE = 4;
  OutPut("TDatabase::ParamDB->SDFEM_TYPE set to be 4"<<endl);
  TDatabase::ParamDB->SC_VERBOSE_AMG = 2;
  OutPut("TDatabase::ParamDB->SC_VERBOSE_AMG set initially to be 2"<<endl);

  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();
  //======================================================================
  // copy read parameters into local variables
  //======================================================================
  RE_NR=TDatabase::ParamDB->RE_NR;
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

  //======================================================================
  // definitions for Navier-Stokes equations
  //======================================================================
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

  RhsArray = new double* [LEVELS+1];
  N_Uarray = new int[LEVELS+1];
  N_Parray = new int[LEVELS+1];

  USpaces = new TFESpace2D*[LEVELS+1];
  PSpaces = new TFESpace2D*[LEVELS+1];
  PsiSpaces = new TFESpace2D*[LEVELS+1];
  VorticitySpaces = new TFESpace2D*[LEVELS+1];
  ProjectionSpaces = new TFESpace2D*[LEVELS+1];

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
  }                                               // endswitch

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
    ||(TDatabase::ParamDB->DISCTYPE==GL00_CONVOLUTION)
    || (TDatabase::ParamDB->DISCTYPE==CLASSICAL_LES)
    ||(TDatabase::ParamDB->CONVOLUTE_SOLUTION)
    ||(TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE==4))
  {
    OutPut("TDatabase::ParamDB->DISCTYPE==GL00_AUX_PROBLEM not implemented" << endl);
    OutPut("TDatabase::ParamDB->DISCTYPE==GL00_ONVOLUTION not implemented" << endl);
    OutPut("TDatabase::ParamDB->DISCTYPE==CLASSICAL_LES not implemented" << endl);
    OutPut("TDatabase::ParamDB->CONVOLUTE_SOLUTION not implemented" << endl);
    OutPut("TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE==4 not implemented" << endl);
    exit(4711);
  }

  downwind = new int*[LEVELS+1];

  //======================================================================
  // creating discrete forms
  //======================================================================

  InitializeDiscreteForms(DiscreteFormGalerkin,DiscreteFormUpwind,
    DiscreteFormSmagorinsky,DiscreteFormColetti,
    DiscreteFormGL00Convolution,DiscreteFormGL00AuxProblem,
    DiscreteFormVMS_Projection,
    DiscreteFormNLGalerkin,
    DiscreteFormNLUpwind, DiscreteFormNLSmagorinsky,
    DiscreteFormNLColetti,DiscreteFormNLGL00Convolution,
    DiscreteFormNLGL00AuxProblem,
    DiscreteFormNLVMS_Projection,
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
  BoundaryConditionsAuxProblem[1]  =BoundConditionAuxProblem;
  BoundaryConditionsAuxProblem[2] = BoundConditionAuxProblem;

  BoundValuesAuxProblem[0] = BoundValueAuxProblem;
  BoundValuesAuxProblem[1] = BoundValueAuxProblem;
  BoundValuesAuxProblem[2] = BoundValueAuxProblem;

  Coefficients[0] = LinCoeffs;                    // Navier-Stokes
  // c_A, c_B
  Coefficients[1] = BilinearCoeffs;
  // c_C
  Coefficients[2] = BilinearCoeffs_Cc;
  Coefficients[3] = NoCoeffs;

  for(i=0;i<TDatabase::ParamDB->SC_COARSEST_LEVEL_SADDLE;i++)
    Domain->RegRefineAll();

  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
  omega = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SADDLE;
  divergence = TDatabase::ParamDB->SC_DIV_FACTOR;

  Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
  Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
  Parameters[2] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;

  TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->STARTTIME;

  if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    i=1;
    MG = new TNSE_MultiGrid(i, N_Paramters, Parameters);
  }

  t3 = GetTime();
  total_time = t3 - total_time;
  mg_level = LEVELS+mg_level;

  //======================================================================
  // definitions for convection-reactions equations
  //======================================================================
  // force direct solver for cd equations
  solver_type_reaction = 2;
  ansatz_order_reaction = 1;

  // array for pointers to the solutions on the
  // different levels of the multigrid
  SolArray_c_A = new TFEFunction2D*[LEVELS+1];
  SolArray_c_B = new TFEFunction2D*[LEVELS+1];
  SolArray_c_C = new TFEFunction2D*[LEVELS+1];
  SolArray_c_A_old = new TFEFunction2D*[LEVELS+1];
  SolArray_c_B_old = new TFEFunction2D*[LEVELS+1];
  SolArray_c_C_old = new TFEFunction2D*[LEVELS+1];
  SolArray_other = new TFEFunction2D*[LEVELS+1];
  SolArray = new TFEFunction2D*[LEVELS+1];
  SolArray_old = new TFEFunction2D*[LEVELS+1];
  SolArray_other_old = new TFEFunction2D*[LEVELS+1];
  IntegralSpaces_c_C_fct = new TFEFunction2D*[LEVELS+1];

  // array for pointers to right hand sides on the
  // different levels of the multigrid
  RhsArray_c = new double* [LEVELS+1];
  RhsArray_c_A = new double* [LEVELS+1];
  RhsArray_c_B = new double* [LEVELS+1];
  RhsArray_c_C = new double* [LEVELS+1];
  N_Array_c_A = new int[LEVELS+1];
  N_Array_c_B = new int[LEVELS+1];
  N_Array_c_C = new int[LEVELS+1];
  N_Array_c = new int[LEVELS+1];

  // array which points to the finite element spaces on the
  // different levels of the multigrid
  ConcentrationSpaces = new TFESpace2D*[LEVELS+1];
  ConcentrationSpaces_other = new TFESpace2D*[LEVELS+1];
  ConcentrationSpaces_c_A = new TFESpace2D*[LEVELS+1];
  ConcentrationSpaces_c_B = new TFESpace2D*[LEVELS+1];
  ConcentrationSpaces_c_C = new TFESpace2D*[LEVELS+1];
  IntegralSpaces_c_C = new TFESpace2D*[LEVELS+1];

  // array which points to the system matrices on  the
  // different levels of the multigrid
  MatricesA_c = new TSquareMatrix2D*[LEVELS+1];
  MatricesA_c_A = new TSquareMatrix2D*[LEVELS+1];
  MatricesA_c_B = new TSquareMatrix2D*[LEVELS+1];
  MatricesA_c_C = new TSquareMatrix2D*[LEVELS+1];

  // array which points to the mass matrices on  the
  // different levels of the multigrid
  MatricesM_c = new TSquareMatrix2D*[LEVELS+1];
  MatricesM_c_A = new TSquareMatrix2D*[LEVELS+1];
  MatricesM_c_B = new TSquareMatrix2D*[LEVELS+1];
  MatricesM_c_C = new TSquareMatrix2D*[LEVELS+1];

  if (TDatabase::ParamDB->BULK_REACTION_DISC == SDFEM)
  {
    // array which points to the stabilization  matrices (sdfem) on the
    // different levels of the multigrid
    MatricesK_c = new TSquareMatrix2D*[LEVELS+1];
    MatricesK_c_A = new TSquareMatrix2D*[LEVELS+1];
    MatricesK_c_B = new TSquareMatrix2D*[LEVELS+1];
    MatricesK_c_C = new TSquareMatrix2D*[LEVELS+1];

    if (TDatabase::ParamDB->SOLD_TYPE)
    {
      // array which points to the stabilization  matrices (sold) on the
      // different levels of the multigrid
      MatricesS_c = new TSquareMatrix2D*[LEVELS+1];
      MatricesS_c_A = new TSquareMatrix2D*[LEVELS+1];
      MatricesS_c_B = new TSquareMatrix2D*[LEVELS+1];
      MatricesS_c_C = new TSquareMatrix2D*[LEVELS+1];
    }
  }
  if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
    (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
  {
    // array which points to the stabilization  matrices (sdfem) on the
    // different levels of the multigrid
    MatricesK_c = new TSquareMatrix2D*[LEVELS+1];
    MatricesK_c_A = new TSquareMatrix2D*[LEVELS+1];
    MatricesK_c_B = new TSquareMatrix2D*[LEVELS+1];
    MatricesK_c_C = new TSquareMatrix2D*[LEVELS+1];
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
  DiscreteFormMatrixMBulk = new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_MatrixM_Bulk, Derivatives_MatrixM_Bulk,
    SpacesNumbers_MatrixM_Bulk, N_Matrices_MatrixM_Bulk, N_Rhs_MatrixM_Bulk,
    RowSpace_MatrixM_Bulk, ColumnSpace_MatrixM_Bulk, RhsSpace_MatrixM_Bulk,
    MatrixMAssemble_Bulk, NoCoeffs, NULL);

  // discrete form for assembling stiffness matrix, stabilization matrix and rhs (SDFEM)
  DiscreteFormMatricesA_SUPG_Bulk = new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_MatricesA_SUPG_Bulk, Derivatives_MatricesA_SUPG_Bulk,
    SpacesNumbers_MatricesA_SUPG_Bulk, N_Matrices_MatricesA_SUPG_Bulk, N_Rhs_MatricesA_SUPG_Bulk,
    RowSpace_MatricesA_SUPG_Bulk, ColumnSpace_MatricesA_SUPG_Bulk, RhsSpace_MatricesA_SUPG_Bulk,
    MatricesA_Assemble_SUPG_Bulk, BilinearCoeffs, NULL);

  DiscreteFormMatricesA_SUPG_Bulk_Cc = new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_MatricesA_SUPG_Bulk, Derivatives_MatricesA_SUPG_Bulk,
    SpacesNumbers_MatricesA_SUPG_Bulk, N_Matrices_MatricesA_SUPG_Bulk, N_Rhs_MatricesA_SUPG_Bulk,
    RowSpace_MatricesA_SUPG_Bulk, ColumnSpace_MatricesA_SUPG_Bulk, RhsSpace_MatricesA_SUPG_Bulk,
    MatricesA_Assemble_SUPG_Bulk, BilinearCoeffs_Cc, NULL);

  DiscreteFormMatricesA_SUPG_Bulk_mom = new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_MatricesA_SUPG_Bulk, Derivatives_MatricesA_SUPG_Bulk,
    SpacesNumbers_MatricesA_SUPG_Bulk, N_Matrices_MatricesA_SUPG_Bulk, N_Rhs_MatricesA_SUPG_Bulk,
    RowSpace_MatricesA_SUPG_Bulk, ColumnSpace_MatricesA_SUPG_Bulk, RhsSpace_MatricesA_SUPG_Bulk,
    MatricesA_Assemble_SUPG_Bulk, BilinearCoeffs_mom, NULL);

  DiscreteFormRhs_SUPG_Bulk = new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_Rhs_SUPG_Bulk, Derivatives_Rhs_SUPG_Bulk,
    SpacesNumbers_Rhs_SUPG_Bulk, N_Matrices_Rhs_SUPG_Bulk, N_Rhs_Rhs_SUPG_Bulk,
    RowSpace_Rhs_SUPG_Bulk, ColumnSpace_Rhs_SUPG_Bulk, RhsSpace_Rhs_SUPG_Bulk,
    Rhs_Assemble_SUPG_Bulk, BilinearCoeffs, NULL);

  DiscreteFormRhs_SUPG_Bulk_Cc = new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_Rhs_SUPG_Bulk, Derivatives_Rhs_SUPG_Bulk,
    SpacesNumbers_Rhs_SUPG_Bulk, N_Matrices_Rhs_SUPG_Bulk, N_Rhs_Rhs_SUPG_Bulk,
    RowSpace_Rhs_SUPG_Bulk, ColumnSpace_Rhs_SUPG_Bulk, RhsSpace_Rhs_SUPG_Bulk,
    Rhs_Assemble_SUPG_Bulk, BilinearCoeffs_Cc, NULL);

  DiscreteFormRhs_SUPG_Bulk_mom = new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_Rhs_SUPG_Bulk, Derivatives_Rhs_SUPG_Bulk,
    SpacesNumbers_Rhs_SUPG_Bulk, N_Matrices_Rhs_SUPG_Bulk, N_Rhs_Rhs_SUPG_Bulk,
    RowSpace_Rhs_SUPG_Bulk, ColumnSpace_Rhs_SUPG_Bulk, RhsSpace_Rhs_SUPG_Bulk,
    Rhs_Assemble_SUPG_Bulk, BilinearCoeffs_mom, NULL);

  // discrete form for assembling stiffness matrix, upwinding
  DiscreteFormMatricesA_Bulk = new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_MatricesA_SUPG_Bulk, Derivatives_MatricesA_SUPG_Bulk,
    SpacesNumbers_MatricesA_SUPG_Bulk, N_Matrices_MatricesA_Galerkin_Bulk, N_Rhs_MatricesA_SUPG_Bulk,
    RowSpace_MatricesA_Galerkin_Bulk, ColumnSpace_MatricesA_Galerkin_Bulk, RhsSpace_MatricesA_SUPG_Bulk,
    MatricesA_Assemble_Bulk, BilinearCoeffs, NULL);

  DiscreteFormMatricesA_Bulk_Cc = new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_MatricesA_SUPG_Bulk, Derivatives_MatricesA_SUPG_Bulk,
    SpacesNumbers_MatricesA_SUPG_Bulk, N_Matrices_MatricesA_Galerkin_Bulk, N_Rhs_MatricesA_SUPG_Bulk,
    RowSpace_MatricesA_Galerkin_Bulk, ColumnSpace_MatricesA_Galerkin_Bulk, RhsSpace_MatricesA_SUPG_Bulk,
    MatricesA_Assemble_Bulk, BilinearCoeffs_Cc, NULL);

  DiscreteFormRhs_Bulk =  new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_Rhs_Galerkin_Bulk, Derivatives_Rhs_Galerkin_Bulk,
    SpacesNumbers_Rhs_Galerkin_Bulk, N_Matrices_Rhs_SUPG_Bulk, N_Rhs_Rhs_SUPG_Bulk,
    RowSpace_Rhs_SUPG_Bulk, ColumnSpace_Rhs_SUPG_Bulk, RhsSpace_Rhs_SUPG_Bulk,
    Rhs_Assemble_Bulk, BilinearCoeffs, NULL);

  DiscreteFormRhs_Bulk_Cc =  new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_Rhs_Galerkin_Bulk, Derivatives_Rhs_Galerkin_Bulk,
    SpacesNumbers_Rhs_Galerkin_Bulk, N_Matrices_Rhs_SUPG_Bulk, N_Rhs_Rhs_SUPG_Bulk,
    RowSpace_Rhs_SUPG_Bulk, ColumnSpace_Rhs_SUPG_Bulk, RhsSpace_Rhs_SUPG_Bulk,
    Rhs_Assemble_Bulk, BilinearCoeffs_Cc, NULL);

  // discrete form for assembling stiffness matrix, Galerkin (FEM-FCT)
  DiscreteFormMatricesA_Galerkin_Bulk = new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_MatricesA_SUPG_Bulk, Derivatives_MatricesA_SUPG_Bulk,
    SpacesNumbers_MatricesA_SUPG_Bulk, N_Matrices_MatricesA_Galerkin_Bulk, N_Rhs_MatricesA_SUPG_Bulk,
    RowSpace_MatricesA_Galerkin_Bulk, ColumnSpace_MatricesA_Galerkin_Bulk, RhsSpace_MatricesA_SUPG_Bulk,
    MatricesA_Assemble_Galerkin_Bulk, BilinearCoeffs, NULL);

  DiscreteFormMatricesA_Galerkin_Bulk_Cc = new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_MatricesA_SUPG_Bulk, Derivatives_MatricesA_SUPG_Bulk,
    SpacesNumbers_MatricesA_SUPG_Bulk, N_Matrices_MatricesA_Galerkin_Bulk, N_Rhs_MatricesA_SUPG_Bulk,
    RowSpace_MatricesA_Galerkin_Bulk, ColumnSpace_MatricesA_Galerkin_Bulk, RhsSpace_MatricesA_SUPG_Bulk,
    MatricesA_Assemble_Galerkin_Bulk, BilinearCoeffs_Cc, NULL);

  DiscreteFormRhs_Galerkin_Bulk =  new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_Rhs_Galerkin_Bulk, Derivatives_Rhs_Galerkin_Bulk,
    SpacesNumbers_Rhs_Galerkin_Bulk, N_Matrices_Rhs_SUPG_Bulk, N_Rhs_Rhs_SUPG_Bulk,
    RowSpace_Rhs_SUPG_Bulk, ColumnSpace_Rhs_SUPG_Bulk, RhsSpace_Rhs_SUPG_Bulk,
    Rhs_Assemble_Bulk, BilinearCoeffs, NULL);

  DiscreteFormRhs_Galerkin_Bulk_Cc =  new TDiscreteForm2D
    (MassMatrix, Mass, N_Terms_Rhs_Galerkin_Bulk, Derivatives_Rhs_Galerkin_Bulk,
    SpacesNumbers_Rhs_Galerkin_Bulk, N_Matrices_Rhs_SUPG_Bulk, N_Rhs_Rhs_SUPG_Bulk,
    RowSpace_Rhs_SUPG_Bulk, ColumnSpace_Rhs_SUPG_Bulk, RhsSpace_Rhs_SUPG_Bulk,
    Rhs_Assemble_Bulk, BilinearCoeffs_Cc, NULL);

  // initialize multigrid
  if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
    Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
    i=1;
    MG_c_A = new TMultiGrid2D(i, N_Paramters, Parameters);
    MG_c_B = new TMultiGrid2D(i, N_Paramters, Parameters);
    MG_c_C = new TMultiGrid2D(i, N_Paramters, Parameters);
  }

  BoundaryConditions_Scalar[0] =  BoundCondition_c_A;
  BoundaryConditions_Scalar[1] =  BoundCondition_c_B;
  BoundaryConditions_Scalar[2] =  BoundCondition_c_C;

  BoundValues_Scalar[0] = BoundValue_c_A;
  BoundValues_Scalar[1] = BoundValue_c_B;
  BoundValues_Scalar[2] = BoundValue_c_C;

  //======================================================================
  // definitions for population balance equation
  //======================================================================
  mom = TDatabase::ParamDB->BULK_METHODS_OF_MOMENTS;
  z_min = TDatabase::ParamDB->BULK_D_P_0/TDatabase::ParamDB->BULK_D_P_MAX;
  average_step[0] = 0;
  average_median[0] = 0.0;
  // method of moments
  if (mom)
  {
    N_mom = 3;
    MOMSpaces = new TFESpace2D*[LEVELS+1];
    SolArray_mom = new TFEVectFunct2D*[LEVELS+1];
    SolArray_mom_old = new TFEVectFunct2D*[LEVELS+1];
    RhsArray_mom = new double*[LEVELS+1];
    N_Array_mom = new int[LEVELS+1];
    MatricesA_mom = new TSquareMatrix2D*[LEVELS+1];
    MatricesM_mom = new TSquareMatrix2D*[LEVELS+1];
    if (TDatabase::ParamDB->BULK_MOM_DISC == SDFEM)
    {
      MatricesK_mom = new TSquareMatrix2D*[LEVELS+1];
      if (TDatabase::ParamDB->SOLD_TYPE)
        MatricesS_mom = new TSquareMatrix2D*[LEVELS+1];
    }
    // initialize multigrid
    if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
    {
      Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
      Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
      i=1;
      MG_mom = new TMultiGrid2D(i, N_Paramters, Parameters);
    }
    BoundaryConditions_Scalar[3] =  BoundCondition_mom;
    BoundValues_Scalar[3] = BoundValue_mom;
  }
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
    OutPut("MEMORY: " << setw(10) << GetMemory() << endl);

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
        new double[40];
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
            pressure_space_code,NULL,downwind[i]);
          break;
      }                                           // end switch(NSTYPE)
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

      case VMS_PROJECTION:
        DiscreteForm = DiscreteFormVMS_Projection;
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
      Output = new TOutput2D(5, 4, 1, 2, Domain);
      // das ist zum Einlesen f"ur "altere Daten
      //Output = new TOutput2D(5, 7, 1, 2, Domain);
      Output->AddFEVectFunct(u);
      Output->AddFEFunction(p);
      os.seekp(std::ios::beg);
      Output->AddParameter(TDatabase::TimeDB->CURRENTTIME,os.str().c_str());

      psi = new double[N_V];
      StreamFct = new TFEFunction2D(streamfunction_space, PsiString,  PsiString,  psi, N_V);
      //Output->AddFEFunction(StreamFct);

      vorticity = new double[N_Vort];
      Vorticity = new TFEFunction2D(vorticity_space, VorticityString, VorticityString, vorticity, N_Vort);
      //ComputeVorticity(U1Array[mg_level-1], U2Array[mg_level-1],Vorticity);
      //Output->AddFEFunction(Vorticity);

      div = new double[N_Vort];
      Divergence = new TFEFunction2D(vorticity_space, DivergenceString, DivergenceString, div, N_Vort);
      //Output->AddFEFunction(Divergence);
    }

    // set rhs
    fesp[0] = velocity_space;
    fesp[1] = pressure_space;
    fesp[2] = convolution_space;

    fefct[0] = u1;
    fefct[1] = u2;

    ferhs[0] = velocity_space;
    ferhs[1] = velocity_space;

    // 3 parameters are needed for assembling
    // which are u1_old, u2_old, norm of grad u_old
    switch(CurrentDiscType)
    {
      // turbulent viscosity must be computed
      case SMAGORINSKY:
      case VMS_PROJECTION:
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
    OutPut("rhs " << Ddot(N_Unknowns, rhs, rhs)<<endl);
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
          UpwindForNavierStokes(Coefficients[0], SQMATRICES[0], U1Array[i], U2Array[i]);
          //cout << "UPWINDING DONE : level " << i << endl;
          break;

        case 3:
        case 4:
          // do upwinding with two matrices
          UpwindForNavierStokes(Coefficients[0], SQMATRICES[0], U1Array[i], U2Array[i]);
          UpwindForNavierStokes(Coefficients[0], SQMATRICES[3], U1Array[i], U2Array[i]);
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
    OutPut("MEMORY: " << setw(10) << GetMemory() << endl);

    // if multiple discretization multilevel method is used
    // get space for low order disc on finest geo grid
    if ((mg_type==1)&&(i<mg_level-1))
    {
      // nonconforming velocity space
      concentration_space = new TFESpace2D(coll,NameString,UString,BoundCondition_c_A,-1,NULL);
      concentration_space_b = new TFESpace2D(coll,NameString,UString,BoundCondition_c_B,-1,NULL);
      concentration_space_c = new TFESpace2D(coll,NameString,UString,BoundCondition_c_C,-1,NULL);
    }
    // standard multigrid or finest level
    // get fe space of high order disc on finest geo grid
    else
    {
      concentration_space = new TFESpace2D(coll, NameString,
        UString, BoundCondition_c_A, ansatz_order_reaction, NULL);
      concentration_space_b = new TFESpace2D(coll, NameString,
        UString, BoundCondition_c_B, ansatz_order_reaction, NULL);
      concentration_space_c = new TFESpace2D(coll, NameString,
        UString, BoundCondition_c_C, ansatz_order_reaction, NULL);
    }
    // array of the fe spaces
    ConcentrationSpaces_c_A[i] = concentration_space;
    ConcentrationSpaces_c_B[i] = concentration_space_b;
    ConcentrationSpaces_c_C[i] = concentration_space_c;
    N_Unknowns_c_A = concentration_space->GetN_DegreesOfFreedom();
    N_Unknowns_c_C = concentration_space_c->GetN_DegreesOfFreedom();
    N_Array_c_A[i] = N_Unknowns_c_A;
    N_Array_c_C[i] = N_Unknowns_c_C;
    if (ConcentrationSpaces_c_B[i]->GetN_DegreesOfFreedom()!=N_Unknowns_c_A)
    {
      OutPut("different number of d.o.f. for c_A and c_B, not yet implemented !!!"<<endl);
      exit(4711);
    }
    N_Unknowns_c_B = N_Unknowns_c_A;
    N_Array_c_B[i] = N_Unknowns_c_B;
    //N_Unknowns_c_C = N_Unknowns_c_A;

    // active dof (i.e. dof without Dirichlet dofs)
    N_Active_c_A = concentration_space->GetActiveBound();
    N_Active_c_B = concentration_space_b->GetActiveBound();
    N_Active_c_C = concentration_space_c->GetActiveBound();
    if (ConcentrationSpaces_c_B[i]->GetActiveBound()!=N_Active_c_A)
    {
      OutPut("different number of active d.o.f. for c_A and c_B, not yet implemented !!!"<<endl);
      exit(4711);
    }
    //N_Active_c_C = N_Active_c_A;
    OutPut("dof c_A: "<< setw(10) << N_Unknowns_c_A << " active: " <<
      N_Active_c_A << endl);
    OutPut("dof c_B: "<< setw(10) << N_Unknowns_c_B << " active: " <<
      N_Active_c_B << endl);
    OutPut("dof c_C: "<< setw(10) << N_Unknowns_c_C << " active: " <<
      N_Active_c_C << endl);

    // **********************************************************************************
    // substance A
    // **********************************************************************************
    // build matrices
    // first build matrix structure
    sqstructure_c_A = new TSquareStructure2D(ConcentrationSpaces_c_A[i]);
    sqstructure_c_A->Sort();

    // two matrices used
    // A contains the non time dependent part of the discretization
    sqmatrixA = new TSquareMatrix2D(sqstructure_c_A);
    MatricesA_c_A[i] = sqmatrixA;

    // M is the mass matrix
    // the iterative solver uses M
    sqmatrixM_c_A = new TSquareMatrix2D(sqstructure_c_A);
    MatricesM_c_A[i] = sqmatrixM_c_A;
    if (TDatabase::ParamDB->BULK_REACTION_DISC == SDFEM)
    {
      // stabilisation matrix K
      sqmatrixK_c_A = new TSquareMatrix2D(sqstructure_c_A);
      MatricesK_c_A[i] = sqmatrixK_c_A;
      if (TDatabase::ParamDB->SOLD_TYPE)
      {
        sqmatrixS_c_A = new TSquareMatrix2D(sqstructure_c_A);
        MatricesS_c_A[i] = sqmatrixS_c_A;
      }
    }

    // allocate array for solution
    sol_c_A = new double[N_Unknowns_c_A];
    memset(sol_c_A, 0, N_Unknowns_c_A*SizeOfDouble);
    current_sol = sol_c_A;
    oldsol_c_A = new double[N_Unknowns_c_A];

    // allocate fe function on the current level
    c_A= new TFEFunction2D(concentration_space, c_A_String, c_A_String, sol_c_A, N_Unknowns_c_A);
    SolArray_c_A[i] = c_A;
    if ((TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)&&
      (TDatabase::ParamDB->SOLD_TYPE))
    {
      c_A_old = new TFEFunction2D(concentration_space, c_A_String, c_A_String, oldsol_c_A, N_Unknowns_c_A);
      SolArray_c_A_old[i] = c_A_old;
    }

    // allocate array on which the time scheme works
    B_c = new double [N_Unknowns_c_A];
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
      MGLevel_c_A = new TMGLevel2D(i, sqmatrixM_c_A, current_B, current_sol,  n_aux, NULL);
      MG_c_A->AddLevel(MGLevel_c_A);
    }
    // interpolate initial condition
    c_A->Interpolate(InitialCondition_c_A);

    // **********************************************************************************
    // substance B
    // **********************************************************************************

    // build matrices
    // first build matrix structure
    sqstructure_c_B = new TSquareStructure2D(ConcentrationSpaces_c_B[i]);
    sqstructure_c_B->Sort();

    // two matrices used
    // A contains the non time dependent part of the discretization
    sqmatrixA = new TSquareMatrix2D(sqstructure_c_B);
    MatricesA_c_B[i] = sqmatrixA;

    // M is the mass matrix
    // the iterative solver uses M
    sqmatrixM_c_B = new TSquareMatrix2D(sqstructure_c_B);
    MatricesM_c_B[i] = sqmatrixM_c_B;
    if(TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)
    {
      // stabilisation matrix K
      sqmatrixK_c_B = new TSquareMatrix2D(sqstructure_c_B);
      MatricesK_c_B[i] = sqmatrixK_c_B;
      if (TDatabase::ParamDB->SOLD_TYPE)
      {
        sqmatrixS_c_B = new TSquareMatrix2D(sqstructure_c_B);
        MatricesS_c_B[i] = sqmatrixS_c_B;
      }
    }

    // allocate array for solution
    sol_c_B = new double[N_Unknowns_c_B];
    memset(sol_c_B, 0, N_Unknowns_c_B*SizeOfDouble);
    current_sol = sol_c_B;
    oldsol_c_B = new double[N_Unknowns_c_B];

    // allocate fe function on the current level
    c_B= new TFEFunction2D(concentration_space_b, c_B_String, c_B_String, sol_c_B, N_Unknowns_c_B);
    SolArray_c_B[i] = c_B;
    if ((TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)&&
      (TDatabase::ParamDB->SOLD_TYPE))
    {
      c_B_old = new TFEFunction2D(concentration_space_b, c_B_String, c_B_String, oldsol_c_B, N_Unknowns_c_B);
      SolArray_c_B_old[i] = c_B_old;
    }

    // allocate array on which the time scheme works
    B_c = new double [N_Unknowns_c_B];
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
      MGLevel_c_B = new TMGLevel2D(i, sqmatrixM_c_B, current_B, current_sol,  n_aux, NULL);
      MG_c_B->AddLevel(MGLevel_c_B);
    }
    // interpolate initial condition
    c_B->Interpolate(InitialCondition_c_B);

    // **********************************************************************************
    // substance C
    // **********************************************************************************

    // build matrices
    // first build matrix structure
    sqstructure_c_C = new TSquareStructure2D(ConcentrationSpaces_c_C[i]);
    sqstructure_c_C->Sort();

    // two matrices used
    // A contains the non time dependent part of the discretization
    sqmatrixA = new TSquareMatrix2D(sqstructure_c_C);
    MatricesA_c_C[i] = sqmatrixA;

    // M is the mass matrix
    // the iterative solver uses M
    sqmatrixM_c_C = new TSquareMatrix2D(sqstructure_c_C);
    MatricesM_c_C[i] = sqmatrixM_c_C;
    if(TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)
    {
      // stabilisation matrix K
      sqmatrixK_c_C = new TSquareMatrix2D(sqstructure_c_C);
      MatricesK_c_C[i] = sqmatrixK_c_C;
      if (TDatabase::ParamDB->SOLD_TYPE)
      {
        sqmatrixS_c_C = new TSquareMatrix2D(sqstructure_c_C);
        MatricesS_c_C[i] = sqmatrixS_c_C;
      }
    }

    // allocate array for solution
    sol_c_C = new double[N_Unknowns_c_C];
    memset(sol_c_C, 0, N_Unknowns_c_C*SizeOfDouble);
    current_sol = sol_c_C;
    oldsol_c_C = new double[N_Unknowns_c_C];

    // allocate fe function on the current level
    c_C = new TFEFunction2D(concentration_space_c, c_C_String, c_C_String, sol_c_C, N_Unknowns_c_C);
    SolArray_c_C[i] = c_C;
    if ((TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)&&
      (TDatabase::ParamDB->SOLD_TYPE))
    {
      c_C_old = new TFEFunction2D(concentration_space_c, c_C_String, c_C_String, oldsol_c_C, N_Unknowns_c_C);
      SolArray_c_C_old[i] = c_C_old;
    }

    // allocate array on which the time scheme works
    B_c = new double [N_Unknowns_c_C];
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
      MGLevel_c_C = new TMGLevel2D(i, sqmatrixM_c_C, current_B, current_sol, n_aux, NULL);
      MG_c_C->AddLevel(MGLevel_c_C);
    }
    // interpolate initial condition
    c_C->Interpolate(InitialCondition_c_C);

    // define bilinear fe space which contains the integrals of the decrease of f
    integral_space_c_C = new TFESpace2D(coll, NameString, UString, BoundCondition_c_C, 1, NULL);
    IntegralSpaces_c_C[i] = integral_space_c_C;
    N_Unknowns_Integral_Space = integral_space_c_C->GetN_DegreesOfFreedom();
    integral_val = new double[N_Unknowns_Integral_Space];
    memset(integral_val,0, N_Unknowns_Integral_Space*SizeOfDouble);
    integral_space_c_C_fct = new TFEFunction2D(integral_space_c_C, c_A_String, c_A_String,
      integral_val,N_Unknowns_Integral_Space);
    IntegralSpaces_c_C_fct[i] = integral_space_c_C_fct;

    // **********************************************************************************
    // general definitions
    // **********************************************************************************

    // allocate rhs
    rhs_c = new double[N_Unknowns_c_A];
    memset(rhs_c, 0, N_Unknowns_c_A*SizeOfDouble);
    RhsArray_c_A[i] = rhs_c;
    rhs_c = new double[N_Unknowns_c_B];
    memset(rhs_c, 0, N_Unknowns_c_B*SizeOfDouble);
    RhsArray_c_B[i] = rhs_c;
    rhs_c = new double[N_Unknowns_c_C];
    memset(rhs_c, 0, N_Unknowns_c_C*SizeOfDouble);
    RhsArray_c_C[i] = rhs_c;

    //======================================================================
    // assembling of mass matrices
    //======================================================================

    // set parameters
    N_Rhs = 0;
    N_FESpaces = 1;
    N_SquareMatrices = 1;
    DiscreteForm = DiscreteFormMatrixMBulk;
    aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

    // substance A
    fesp[0] = ConcentrationSpaces_c_A[i];
    SQMATRICES[0] = MatricesM_c_A[i];
    SQMATRICES[0]->Reset();
    
    if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
	(TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
	BoundValues_Scalar[0] = BoundValue_FEM_FCT;

    Assemble2D(N_FESpaces, fesp,
      N_SquareMatrices, SQMATRICES,
      0, NULL,
      0, NULL, NULL,
      DiscreteForm,
      BoundaryConditions_Scalar,
      BoundValues_Scalar,
      aux);

    if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
	(TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
	BoundValues_Scalar[0] = BoundValue_c_A;

    if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
      (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)
      && (i==mg_level-1))
    {
      CheckWrongNeumannNodes_c_A(coll, ConcentrationSpaces_c_A[i], N_neum_to_diri_c_A, neum_to_diri_c_A,
        neum_to_diri_bdry_c_A, neum_to_diri_param_c_A);
    }
    if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
      (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)
      && (i==mg_level-1))
    {
      oldrhs_fem_fct0_c_A = new double[N_Unknowns_c_A];
      memcpy(oldrhs_fem_fct0_c_A,  RhsArray_c_A[i], N_Unknowns_c_A*SizeOfDouble);
      oldrhs_fem_fct1_c_A = new double[N_Unknowns_c_A];
      lump_mass_c_A = new double [N_Unknowns_c_A];
      LumpMassMatrixToVector(MatricesM_c_A[i], lump_mass_c_A);
      matrix_D_Entries_c_A = new double[MatricesA_c_A[i]->GetN_Entries()];
      // matrix K for copy of mass matrix
      sqmatrixK_c_A = new TSquareMatrix2D(sqstructure_c_A);
      MatricesK_c_A[i] = sqmatrixK_c_A;
      tilde_u_c_A = new double [N_Unknowns_c_A];
      // save mass matrix in matricesK
      memcpy(MatricesK_c_A[i]->GetEntries(), MatricesM_c_A[i]->GetEntries(),
        MatricesM_c_A[i]->GetN_Entries() * SizeOfDouble);

    }
    // if (TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING)
    //     LumpMassMatrixToDiag_Bulk(MatricesM_c_A[i]);

    // substance B
    fesp[0] = ConcentrationSpaces_c_B[i];
    SQMATRICES[0] = MatricesM_c_B[i];
    SQMATRICES[0]->Reset();

    if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
	(TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
	BoundValues_Scalar[1] = BoundValue_FEM_FCT;

    Assemble2D(N_FESpaces, fesp,
      N_SquareMatrices, SQMATRICES,
      0, NULL,
      0, NULL, NULL,
      DiscreteForm,
      BoundaryConditions_Scalar+1,
      BoundValues_Scalar+1,
      aux);

    if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
	(TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
	BoundValues_Scalar[1] = BoundValue_c_B;

    if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
      (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)
      && (i==mg_level-1))
    {
      CheckWrongNeumannNodes_c_B(coll, ConcentrationSpaces_c_B[i], N_neum_to_diri_c_B, neum_to_diri_c_B,
        neum_to_diri_bdry_c_B, neum_to_diri_param_c_B);
    }
    if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
      (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)
      && (i==mg_level-1))
    {
      oldrhs_fem_fct0_c_B = new double[N_Unknowns_c_B];
      memcpy(oldrhs_fem_fct0_c_B,  RhsArray_c_B[i], N_Unknowns_c_B*SizeOfDouble);
      oldrhs_fem_fct1_c_B = new double[N_Unknowns_c_B];
      lump_mass_c_B = new double [N_Unknowns_c_B];
      LumpMassMatrixToVector(MatricesM_c_B[i], lump_mass_c_B);
      matrix_D_Entries_c_B = new double[MatricesA_c_B[i]->GetN_Entries()];
      // matrix K for copy of mass matrix
      sqmatrixK_c_B = new TSquareMatrix2D(sqstructure_c_B);
      MatricesK_c_B[i] = sqmatrixK_c_B;
      tilde_u_c_B = new double [N_Unknowns_c_B];
      memcpy(MatricesK_c_B[i]->GetEntries(), MatricesM_c_B[i]->GetEntries(),
        MatricesM_c_B[i]->GetN_Entries() * SizeOfDouble);
    }

    // if (TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING)
    //     LumpMassMatrixToDiag_Bulk(MatricesM_c_B[i]);

    // substance C
    fesp[0] = ConcentrationSpaces_c_C[i];
    SQMATRICES[0] = MatricesM_c_C[i];
    SQMATRICES[0]->Reset();

    if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
	(TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
	BoundValues_Scalar[2] = BoundValue_FEM_FCT;    

    Assemble2D(N_FESpaces, fesp,
      N_SquareMatrices, SQMATRICES,
      0, NULL,
      0, NULL, NULL,
      DiscreteForm,
      BoundaryConditions_Scalar+2,
      BoundValues_Scalar+2,
      aux);

    if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
	(TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
	BoundValues_Scalar[2] = BoundValue_c_C;    

    if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
      (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)
      && (i==mg_level-1))
    {
      CheckWrongNeumannNodes_c_C(coll, ConcentrationSpaces_c_C[i], N_neum_to_diri_c_C, neum_to_diri_c_C,
        neum_to_diri_bdry_c_C, neum_to_diri_param_c_C);
    }
    if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
      (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)
      && (i==mg_level-1))
    {
      oldrhs_fem_fct0_c_C = new double[N_Unknowns_c_C];
      memcpy(oldrhs_fem_fct0_c_C,  RhsArray_c_C[i], N_Unknowns_c_C*SizeOfDouble);
      oldrhs_fem_fct1_c_C = new double[N_Unknowns_c_C];
      lump_mass_c_C = new double [N_Unknowns_c_C];
      LumpMassMatrixToVector(MatricesM_c_C[i], lump_mass_c_C);
      matrix_D_Entries_c_C = new double[MatricesA_c_C[i]->GetN_Entries()];
      // matrix K for copy of mass matrix
      sqmatrixK_c_C = new TSquareMatrix2D(sqstructure_c_C);
      MatricesK_c_C[i] = sqmatrixK_c_C;
      tilde_u_c_C= new double [N_Unknowns_c_C];
      memcpy(MatricesK_c_C[i]->GetEntries(), MatricesM_c_C[i]->GetEntries(),
        MatricesM_c_C[i]->GetN_Entries() * SizeOfDouble);
    }

    delete aux;

    // if (TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING)
    //     LumpMassMatrixToDiag_Bulk(MatricesM_c_C[i]);

    // **********************************************************************************
    // MOM
    // **********************************************************************************
    if (mom)
    {
      // fe space
      mom_space = new TFESpace2D(coll, NameString,
        UString, BoundCondition_mom,
        ansatz_order_reaction, NULL);
      MOMSpaces[i] = mom_space;
      N_Unknowns_mom = MOMSpaces[i]->GetN_DegreesOfFreedom();
      N_Active_mom = MOMSpaces[i]->GetActiveBound();

      // build matrices
      // first build matrix structure
      sqstructure_mom = new TSquareStructure2D(MOMSpaces[i]);
      sqstructure_mom->Sort();

      // two matrices used
      // A contains the non time dependent part of the discretization
      sqmatrixA = new TSquareMatrix2D(sqstructure_mom);
      MatricesA_mom[i] = sqmatrixA;

      // M is the mass matrix
      // the iterative solver uses M
      sqmatrixM_mom = new TSquareMatrix2D(sqstructure_mom);
      MatricesM_mom[i] = sqmatrixM_mom;
      if(TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)
      {
        // stabilisation matrix K
        sqmatrixK_mom = new TSquareMatrix2D(sqstructure_mom);
        MatricesK_mom[i] = sqmatrixK_mom;
        if (TDatabase::ParamDB->SOLD_TYPE)
        {
          sqmatrixS_mom = new TSquareMatrix2D(sqstructure_mom);
          MatricesS_mom[i] = sqmatrixS_mom;
        }
      }

      // allocate array for solution
      sol_mom = new double[N_Unknowns_mom*N_mom];
      memset(sol_mom, 0, N_Unknowns_mom*N_mom*SizeOfDouble);
      current_sol = sol_mom;
      oldsol_mom = new double[N_Unknowns_mom*N_mom];
      memset(oldsol_mom, 0, N_Unknowns_mom*N_mom*SizeOfDouble);

      // allocate fe function on the current level
      fe_mom = new TFEVectFunct2D(MOMSpaces[i], c_C_String, c_C_String, sol_mom, N_Unknowns_mom, N_mom);
      SolArray_mom[i] = fe_mom;
      if ((TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)&&
        (TDatabase::ParamDB->SOLD_TYPE))
      {
        fe_mom_old = new TFEVectFunct2D(MOMSpaces[i], c_C_String, c_C_String, sol_mom, N_Unknowns_mom, N_mom);
        SolArray_mom_old[i] = fe_mom_old;
      }

      // allocate array on which the time scheme works
      B_mom = new double [N_Unknowns_mom];
      current_B = B_mom;
      rhs_single_mom = new double [N_Unknowns_mom];
      memset(rhs_single_mom, 0, N_Unknowns_mom*SizeOfDouble);

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
        MGLevel_mom = new TMGLevel2D(i, sqmatrixM_mom, current_B, current_sol, n_aux, NULL);
        MG_mom->AddLevel(MGLevel_mom);
      }

      // interpolate initial condition
      fe_mom->Interpolate(InitialCondition_mom);

      // allocate rhs
      rhs_mom = new double[N_Unknowns_mom*N_mom];
      memset(rhs_mom, 0, N_Unknowns_mom*N_mom*SizeOfDouble);
      RhsArray_mom[i] = rhs_mom;
      oldrhs_mom = new double[N_Unknowns_mom*N_mom];
      memset(oldrhs_mom, 0, N_Unknowns_mom*N_mom*SizeOfDouble);
      defect_mom = new double[N_Unknowns_mom];
      memset(defect_mom, 0, N_Unknowns_mom*SizeOfDouble);

      //======================================================================
      // assembling of mass matrices
      //======================================================================

      // set parameters
      N_Rhs = 0;
      N_FESpaces = 1;
      N_SquareMatrices = 1;
      DiscreteForm = DiscreteFormMatrixMBulk;
      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

      // same for all moments
      fesp[0] = MOMSpaces[i];
      SQMATRICES[0] = MatricesM_mom[i];
      SQMATRICES[0]->Reset();

      Assemble2D(N_FESpaces, fesp,
        N_SquareMatrices, SQMATRICES,
        0, NULL,
        0, NULL, NULL,
        DiscreteForm,
        BoundaryConditions_Scalar+3,
        BoundValues_Scalar+3,
        aux);

    }                                             // end if (mom)

    // prepare output
    if (i==mg_level-1)
    {
      Output->AddFEFunction(c_A);
      Output->AddFEFunction(c_B);
      Output->AddFEFunction(c_C);
      os.seekp(std::ios::beg);
      Output->AddParameter(real_time,os.str().c_str());
    }

    // read initial solution of finest level from grape file
    if ((i==mg_level-1)&&(TDatabase::ParamDB->READ_GRAPE_FILE))
    {
      // version without vorticity and divergence
      // AuxFEFunctArray = new TFEFunction2D*[2];
      // AuxFEFunctArray[0] = PArray[mg_level-1];
      // AuxFEFunctArray[1] = StreamFct;
      // AuxFEVectFunctArray = new TFEVectFunct2D*[2];
      // AuxFEVectFunctArray[0] = UArray[mg_level-1];
      // AuxFEVectFunctArray[1] = uconf;
      // ReadGrapeFile(ReadGrapeBaseName, 2 , 2, AuxFEFunctArray,AuxFEVectFunctArray);

      // version with vorticity and divergence
      AuxFEFunctArray = new TFEFunction2D*[7];
      AuxFEFunctArray[0] = PArray[mg_level-1];
      AuxFEFunctArray[1] = StreamFct;
      AuxFEFunctArray[2] = Vorticity;
      AuxFEFunctArray[3] = Divergence;
      AuxFEFunctArray[4] = c_A;
      AuxFEFunctArray[5] = c_B;
      AuxFEFunctArray[6] = c_C;
      AuxFEVectFunctArray = new TFEVectFunct2D*[1];
      AuxFEVectFunctArray[0] = UArray[mg_level-1];
      //AuxFEVectFunctArray[1] = uconf;
      ReadGrapeFile(ReadGrapeBaseName, 7 , 1, AuxFEFunctArray,AuxFEVectFunctArray);
      // right hand side of first time step will not be completely correct
      // for GL00--discs and Coletti--disc
      // since to reconstruct the rhs of former time step, the initial
      // velocity of this time step is necessary
      // instead final velocity is used now
      if (TDatabase::TimeDB->RESET_CURRENTTIME > 0)
      {
        TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->RESET_CURRENTTIME_STARTTIME;
        OutPut("start time reset to " << TDatabase::TimeDB->CURRENTTIME << endl);
        memset(sol_c_A, 0, N_Unknowns_c_A*SizeOfDouble);
        memset(sol_c_B, 0, N_Unknowns_c_B*SizeOfDouble);
        memset(sol_c_C, 0, N_Unknowns_c_C*SizeOfDouble);
        OutPut("concentrations set to zero" << endl);
      }
    }
  }                                               // endfor i

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
  if ((TDatabase::ParamDB->WRITE_GRAPE)||(TDatabase::ParamDB->WRITE_GNU)||
    (TDatabase::ParamDB->WRITE_GMV)||(TDatabase::ParamDB->WRITE_VTK))
  {
    StreamFunction(USpaces[mg_level-1], sol, sol+N_Uarray[mg_level-1],
      PsiSpaces[LEVELS-1], psi);
    ComputeVorticityDivergence(USpaces[mg_level-1], U1Array[mg_level-1], U2Array[mg_level-1],
      vorticity_space,vorticity,div);
    comp_vort++;
  }

  defect = new double[N_Unknowns];
  startsol = new double[N_Unknowns];
  frac_step_sol = new double[N_Unknowns];
  oldrhs =  new double[N_Unknowns];

  N_Active = velocity_space->GetActiveBound();

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

  // definitions for convectio-reaction equations and
  // population balance equation

  // allocate arrays for solver
  defect_c = new double[N_Unknowns_c_A];
  startsol_c = new double[N_Unknowns_c_A];
  frac_step_sol_c = new double[N_Unknowns_c_A];
  oldrhs_c_A =  new double[N_Unknowns_c_A];
  oldrhs_c_B =  new double[N_Unknowns_c_B];
  oldrhs_c_C =  new double[N_Unknowns_c_C];
  rhs_c_complete_A  = new double[N_Unknowns_c_A];
  rhs_c_complete_B  = new double[N_Unknowns_c_B];
  memset(oldrhs_c_A, 0, N_Unknowns_c_A*SizeOfDouble);
  memset(oldrhs_c_B, 0, N_Unknowns_c_B*SizeOfDouble);
  memset(oldrhs_c_C, 0, N_Unknowns_c_C*SizeOfDouble);
  memcpy(oldsol_c_A,sol_c_A,N_Unknowns_c_A*SizeOfDouble);
  memcpy(oldsol_c_B,sol_c_B,N_Unknowns_c_B*SizeOfDouble);
  memcpy(oldsol_c_C,sol_c_C,N_Unknowns_c_C*SizeOfDouble);

  // number of active d.o.f.
  N_Active_c_A = ConcentrationSpaces_c_A[mg_level-1]->GetActiveBound();
  N_Unknowns_c_A = ConcentrationSpaces_c_A[mg_level-1]->GetN_DegreesOfFreedom();

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
        prec_c_A = new TMultiGridScaIte(MatVectScalar, DefectScalar, NULL,
          0, N_Unknowns_c_A, MG_c_A, zerostart);
        prec_c_B = new TMultiGridScaIte(MatVectScalar, DefectScalar, NULL,
          0, N_Unknowns_c_B, MG_c_B, zerostart);
        prec_c_C = new TMultiGridScaIte(MatVectScalar, DefectScalar, NULL,
          0, N_Unknowns_c_C, MG_c_C, zerostart);
        if (mom)
        {
          prec_mom = new TMultiGridScaIte(MatVectScalar, DefectScalar, NULL,
            0, N_Unknowns_mom, MG_mom, zerostart);
        }
        break;
      default:
        OutPut("Unknown preconditioner !!!" << endl);
        exit(4711);
    }
    switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
    {
      // fixed point iteration
      case 11:
        itmethod_c_A = new TFixedPointIte(MatVectScalar, DefectScalar, prec_c_A,
          0, N_Unknowns_c_A, 0);
        itmethod_c_B = new TFixedPointIte(MatVectScalar, DefectScalar, prec_c_B,
          0, N_Unknowns_c_B, 0);
        itmethod_c_C = new TFixedPointIte(MatVectScalar, DefectScalar, prec_c_C,
          0, N_Unknowns_c_C, 0);
        if (mom)
        {
          itmethod_mom = new TFixedPointIte(MatVectScalar, DefectScalar, prec_mom,
            0, N_Unknowns_mom, 0);
        }
        if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
        {
          itmethod_sol_c_A = new double[N_Unknowns_c_A];
          itmethod_sol_c_B = itmethod_sol_c_A;
          itmethod_sol_c_C = itmethod_sol_c_A;
          itmethod_rhs_c = new double[N_Unknowns_c_A];
          if (mom)
          {
            itmethod_sol_mom =  new double[N_Unknowns_mom];
            itmethod_rhs_mom =  new double[N_Unknowns_mom];
          }
        }
        else
        {
          itmethod_sol_c_A = sol_c_A;
          itmethod_sol_c_B = sol_c_B;
          itmethod_sol_c_C = sol_c_C;
          itmethod_rhs_c = rhs;
          if (mom)
          {
            itmethod_sol_mom =  sol_mom;
            itmethod_rhs_mom =  rhs_mom;
          }
        }
        break;
      case 16:
        // FGMRES
        itmethod_c_A = new TFgmresIte(MatVectScalar, DefectScalar, prec_c_A,
          0, N_Unknowns_c_A, 1);
        itmethod_c_B = new TFgmresIte(MatVectScalar, DefectScalar, prec_c_B,
          0, N_Unknowns_c_B, 1);
        itmethod_c_C = new TFgmresIte(MatVectScalar, DefectScalar, prec_c_C,
          0, N_Unknowns_c_C, 1);
        if (mom)
        {
          itmethod_mom = new TFgmresIte(MatVectScalar, DefectScalar, prec_mom,
            0, N_Unknowns_mom, 1);
        }
        if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
        {
          itmethod_sol_c_A = new double[N_Unknowns_c_A];
          itmethod_sol_c_B = itmethod_sol_c_A;
          itmethod_sol_c_C = itmethod_sol_c_A;
          itmethod_rhs_c = new double[N_Unknowns_c_A];
          if (mom)
          {
            itmethod_sol_mom =  new double[N_Unknowns_mom];
            itmethod_rhs_mom =  new double[N_Unknowns_mom];
          }
        }
        else
        {
          itmethod_sol_c_A = sol_c_A;
          itmethod_sol_c_B = sol_c_B;
          itmethod_sol_c_C = sol_c_C;
          itmethod_rhs_c = rhs;
          if (mom)
          {
            itmethod_sol_mom =  sol_mom;
            itmethod_rhs_mom =  rhs_mom;
          }
        }
        break;
      default:
        OutPut("Unknown solver !!!" << endl);
        exit(4711);
    }
  }
  else
    low_scalar = mg_level -1;

  if (solver_type_reaction==2)
  {
    // for damping
    itmethod_sol_c_A = new double[N_Unknowns_c_A];
    itmethod_sol_c_B = itmethod_sol_c_A;
    itmethod_sol_c_C = itmethod_sol_c_A;
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
      OutPut("delta is non--constant" << endl);
    }
  }

  limit_c = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR;
  Max_It_c = TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR;

  // solving equation for PSD with population balance equation
  if (!mom)
  {
    //declarations for the population balance equation
    // number of mesh cells in x- and y-direction,  N = (int)sqrt(coll->GetN_Cells()+1e-4)+1;
    //Nx = (int)(sqrt(2.0)/hmin+1e-4);
    Nx = (int)(sqrt(coll->GetN_Cells()+1e-6));
    Ny = Nx;

    //number of layers in z-direction
    Nz = TDatabase::ParamDB->N_CELL_LAYERS;
    OutPut("grid for population balance: " << Nx << " x " << Ny << " x " << Nz << endl);
    // total number of nodes or grid points
    Nodes = (Nx+1)*(Ny+1)*(Nz+1);
    delta_time_step = TDatabase::TimeDB->TIMESTEPLENGTH;

    f_time_space_approx = new double[Nodes];
    memset(f_time_space_approx, 0, Nodes*SizeOfDouble);
    f_time_space_approx_old = new double[Nodes];
    memset(f_time_space_approx_old, 0, Nodes*SizeOfDouble);

    velo1 = new double[(Nx+1)*(Ny+1)];
    memset(velo1, 0, (Nx+1)*(Ny+1)*SizeOfDouble);
    velo1[0] = -4711;
    velo2 = new double[(Nx+1)*(Ny+1)];
    memset(velo2, 0, (Nx+1)*(Ny+1)*SizeOfDouble);
    concent_C_array = new double[(Nx+1)*(Ny+1)];

    x_coord = new double[Nodes];
    memset(x_coord,0,Nodes*SizeOfDouble);
    y_coord = new double[Nodes];
    memset(y_coord,0,Nodes*SizeOfDouble);
    z_coord = new double[Nodes];
    memset(z_coord,0,Nodes*SizeOfDouble);

    z_layers_coord = new double[Nz+1];
    memset(z_layers_coord,0,(Nz+1)*SizeOfDouble);

    // generation of the grid
    // outputs are x_coord, y_coord, z_coord, z_layers_coord
    OutPut("grid generation for population balance equation" << endl);

    grid_generator_3d(x_min, x_max, Nx,
      y_min, y_max, Ny,
      z_min, z_max, Nz,
      x_coord, y_coord, z_coord,
      z_layers_coord);

    if ((TDatabase::ParamDB->BULK_PB_DISC==BULK_BWE_FEM_SUPG) ||
      (TDatabase::ParamDB->BULK_PB_DISC==BULK_BWE_FDM_UPWIND)||
      (TDatabase::ParamDB->BULK_PB_DISC==BULK_FEM_FCT))
    {
      // declaration and initialisation of the column and the row pointer
      col_ptr = new int[27*Nodes];
      memset(col_ptr,0,(27*Nodes)*SizeOfInt);
      row_ptr = new int[Nodes+1];
      memset(row_ptr,0,(Nodes+1)*SizeOfInt);

      // filling of the column and the row pointer
      filling_row_and_col_ptr(&N_Entries, Nodes, Nx, Ny, x_max, x_min, y_max,
        y_min, z_max, z_min, x_coord, y_coord, z_coord,
        row_ptr, col_ptr);

      OutPut("numerate collection" << endl);

      // number of mesh cells -> necessary for the corresponding 2d grid
      N_Cells = coll->GetN_Cells();

      // declaration, initialisation and computing of the corresponding 2d grid
      correspond_2dgrid = new int[N_Cells];
      memset(correspond_2dgrid,0,N_Cells*SizeOfInt);

      generate_correspond_2d_grid(Nx, Ny, x_coord, y_coord, coll, correspond_2dgrid);
      OutPut("matrix structure" << endl);
      matrix_structure = new TSquareStructure2D(Nodes, N_Entries, col_ptr, row_ptr);
      mat = new TSquareMatrix2D(matrix_structure);
      if((TDatabase::ParamDB->BULK_PB_DISC_STAB == GALERKIN) &&
        (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
      {
        matM = new TSquareMatrix2D(matrix_structure);
        matM_cons = new TSquareMatrix2D(matrix_structure);
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
    if (TDatabase::ParamDB->BULK_PB_DISC==BULK_FWE_FDM_UPWIND)
    {
      // number of mesh cells -> necessary for the corresponding 2d grid
      N_Cells = coll->GetN_Cells();

      // declaration, initialisation and computing of the corresponding 2d grid
      correspond_2dgrid = new int[N_Cells];
      memset(correspond_2dgrid,0,N_Cells*SizeOfInt);

      generate_correspond_2d_grid(Nx, Ny, x_coord, y_coord, coll, correspond_2dgrid);
    }

  }                                               // end if (!mom)

  // read initial solution of finest level from grape file
  if (TDatabase::ParamDB->READ_DATA)
  {
    save_sol[0] = UArray[mg_level-1]->GetValues();
    save_sol[1] = SolArray_c_A[mg_level-1]->GetValues();
    save_sol[2] = SolArray_c_B[mg_level-1]->GetValues();
    save_sol[3] = SolArray_c_C[mg_level-1]->GetValues();
    save_sol[4] = f_time_space_approx;
    save_N_Unknowns[0] = 2*UArray[mg_level-1]->GetLength()+PArray[mg_level-1]->GetLength();
    save_N_Unknowns[1] = SolArray_c_A[mg_level-1]->GetLength();
    save_N_Unknowns[2] = SolArray_c_B[mg_level-1]->GetLength();
    save_N_Unknowns[3] = SolArray_c_C[mg_level-1]->GetLength();
    save_N_Unknowns[4] = Nodes;
    ReadData(ReadDataFileName,5,save_sol,save_N_Unknowns);
  }

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
    if (!mom)
    {
      os.seekp(std::ios::beg);
      os << VtkBaseName << "psd." << 0 << ".vtk" << ends;
      write_vtk_psd(Nx, Ny, Nz, x_coord, y_coord, z_coord, f_time_space_approx,os.str().c_str());
    }
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

  // definitions for time stepping scheme
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
        if (very_first_time)
          oldtau=tau;
        // working array for rhs is B, initialize B
        memset(B, 0, N_Unknowns*SizeOfDouble);

        // old rhs multiplied with current subtime step and theta3 on B
        Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, rhs, B);
        Daxpy(N_Active, tau*TDatabase::TimeDB->THETA3, rhs+N_U, B+N_U);

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
        memset(RHSs[0], 0,
          (2*N_Uarray[mg_level-1]+N_Parray[mg_level-1])*SizeOfDouble);

        N_Rhs = 2;
        N_FESpaces = 1;
        N_SquareMatrices = 0;
        N_RectMatrices = 0;
        fefct[0] = U1Array[mg_level-1];
        fefct[1] = U2Array[mg_level-1];

        switch(TDatabase::ParamDB->DISCTYPE)
        {
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
        for(i=0;i<mg_level;i++)
        {
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
            case 3:
              Dscal(MatricesB1[i]->GetN_Entries(),
                tau/oldtau,
                MatricesB1[i]->GetEntries());
              Dscal(MatricesB2[i]->GetN_Entries(),
                tau/oldtau,
                MatricesB2[i]->GetEntries());
              break;

            case 2:
            case 4:
              Dscal(MatricesB1T[i]->GetN_Entries(),
                tau/oldtau,
                MatricesB1T[i]->GetEntries());
              Dscal(MatricesB2T[i]->GetN_Entries(),
                tau/oldtau,
                MatricesB2T[i]->GetEntries());
              Dscal(MatricesB1[i]->GetN_Entries(),
                tau/oldtau,
                MatricesB1[i]->GetEntries());
              Dscal(MatricesB2[i]->GetN_Entries(),
                tau/oldtau,
                MatricesB2[i]->GetEntries());
              break;
          }
        }                                         // endfor

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
          }                                       // endswitch
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
            Daxpy(N_Active, 1, defect, B);
            MatVectActive(MatricesM[mg_level-1], sol+N_U, defect+N_U);
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

        // extrapolate solution to get starting value at next time step
        // current solution sol
        // solution of last time step sol_timestep_m1

        // save sol
        memcpy(oldsol, sol, SizeOfDouble*N_Unknowns);
        //         (take solution of previous discrete time as start solution)
        tau2 = TDatabase::TimeDB->EXTRAPOLATE_WEIGHT*tau/oldtau;
        tau1 = 1 + tau2;
        // at first time step: sol = sol_timestep_m1 -> result is sol
        for (k=0;k<2*N_U;k++)
          sol[k] = tau1*sol[k] - tau2*sol_timestep_m1[k];
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
          }                                       // endswitch
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

            case VMS_PROJECTION:
              DiscreteForm = DiscreteFormNLVMS_Projection;
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

          if(DiscreteForm == DiscreteFormNLUpwind)
          {
            switch(TDatabase::ParamDB->NSTYPE)
            {
              case 1:
              case 2:
                // do upwinding with one matrix
                UpwindForNavierStokes(Coefficients[0], SQMATRICES[0], U1Array[i], U2Array[i]);
                //cout << "UPWINDING DONE : level " << i << endl;
                break;

              case 3:
              case 4:
                // do upwinding with two matrices
                UpwindForNavierStokes(Coefficients[0], SQMATRICES[0], U1Array[i], U2Array[i]);
                UpwindForNavierStokes(Coefficients[0], SQMATRICES[last_sq], U1Array[i], U2Array[i]);
                //cout << "UPWINDING DONE(2) : level " << i << endl;
                //cout << "check correct sqmatrix !!!! " << endl;
                break;
            }                                     // endswitch
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
        }                                         // endfor i

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
          }                                       // endswitch
        }
        // set current factor of steady state matrix
        gamma = tau*TDatabase::TimeDB->THETA1;

        //========================================================================
        // end assembling of system matrix
        //========================================================================

        real_time = TDatabase::TimeDB->CURRENTTIME * TDatabase::ParamDB->BULK_l_infty /
          TDatabase::ParamDB->BULK_u_infty;
        OutPut(endl << "CURRENT TIME: ");
        OutPut(TDatabase::TimeDB->CURRENTTIME);
        OutPut(" (real time: " << real_time << " s)" << endl);

        OutPut(endl << "CURRENT TIME: ");
        OutPut(TDatabase::TimeDB->CURRENTTIME );
        OutPut(" MEMORY: " << setw(10) << GetMemory() << endl);

        //======================================================================
        // nonlinear loop
        //======================================================================
        N_LinIterCurr = 0;
        solver_time_curr = 0;

        for(j=0;j<Max_It;j++)                     // solve nonlinear equation
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
          /*OutPut("solu "<<  Ddot(2*N_U, sol, sol) << " " <<  Ddot(N_P, sol+2*N_U, sol+2*N_U) << endl);
          OutPut("solu "<<  Ddot(N_Unknowns, sol, sol) << endl);
          OutPut("Ma "<<  Ddot(MatricesM[mg_level-1]->GetN_Entries(), MatricesM[mg_level-1]->GetEntries(), MatricesM[mg_level-1]->GetEntries()) << endl);*/

          Defect(sqmatrices,matrices,sol,B,defect);
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
            OutPut("gamma " << gamma << endl);
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
                  Error("AMG does not work for NSTYPE = 4." << endl);
                  return -1;
                  break;
              }
              t2 = GetTime();
              solver_time_curr = t2-t1;
              solver_time += solver_time_curr;
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
              solver_time += solver_time_curr;

              // p1 = 0;
              for(k=0;k<N_Unknowns;k++)
              {
                p2 = sol[k]-oldsol[k];
                sol[k] = oldsol[k] + omega * p2;
                // p1 += p2*p2;
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
                MatAdd(MatricesM21[i], MatricesA21[i], -gamma);
                MatAdd(MatricesM22[i], MatricesA22[i], -gamma);
                break;
            }                                     // endswitch
          }                                       // endfor i
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
                DiscreteForm = DiscreteFormNLVMS_Projection;
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
                  //fefct[2] = u1ConvArray[mg_level-1];
                  //fefct[3] = u2ConvArray[mg_level-1];

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
                  UpwindForNavierStokes(Coefficients[0], SQMATRICES[0], U1Array[i], U2Array[i]);
                  //cout << "UPWINDING DONE : level " << i << endl;
                  break;

                case 3:
                case 4:
                  // do upwinding with two matrices
                  UpwindForNavierStokes(Coefficients[0], SQMATRICES[0], U1Array[i], U2Array[i]);
                  UpwindForNavierStokes(Coefficients[0], SQMATRICES[last_sq], U1Array[i], U2Array[i]);
                  //cout << "UPWINDING DONE(2) : level " << i << endl;
                  //cout << "check correct sqmatrix !!!! " << endl;
                  break;
              }                                   // endswitch
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
          }                                       // endfor i
          // set current factor of steady state matrix
          gamma = tau*TDatabase::TimeDB->THETA1;

        }                                         // endfor Max_It (solution of nonlinear equation)

        //======================================================================
        // end of nonlinear loop for Navier-Stokes equations
        //======================================================================

        IntoL20FEFunction(PArray[mg_level-1]->GetValues(), N_Parray[mg_level-1],
          PSpaces[mg_level-1], velocity_space_code, pressure_space_code);

      }                                           // endfor l (sub steps of fractional step theta)
      //======================================================================
      // computation of reaction
      //======================================================================

      /*real_time = TDatabase::TimeDB->CURRENTTIME * TDatabase::ParamDB->BULK_l_infty /
        TDatabase::ParamDB->BULK_u_infty;
      OutPut(endl << "CURRENT TIME: ");
      OutPut(TDatabase::TimeDB->CURRENTTIME);
      OutPut(" (real time: " << real_time << " s)" << endl);
      */
      nonlin_ite = 0;
      // nonlinear iteration
      // stopping criterion is tested at the end of the cycle
      // therefore gamma_c = 0 when leaving the cycle and the
      // mass matrix is correct for the next time step

      while ( !((nonlin_resid <= limit_c)&& (substance==1)) && (nonlin_ite < Max_It_c))
      {
        switch(substance)
        {
          // substance A was computed, now compute substance B
          case 0:
            N_Unknowns_c = N_Unknowns_c_B;
            N_Active_c = N_Active_c_B;
            sol_c = sol_c_B;
            oldsol_c = oldsol_c_B;
            current_sol_c = sol_c_B;
            itmethod_c = itmethod_c_B;
            itmethod_sol_c = itmethod_sol_c_B;
            RhsArray_c = RhsArray_c_B;
            oldrhs_c = oldrhs_c_B;
            rhs_c_complete = rhs_c_complete_B;
            Coeff_c =  Coefficients[1];
            N_neum_to_diri_c = N_neum_to_diri_c_B;
            neum_to_diri_c = neum_to_diri_c_B;
            neum_to_diri_bdry_c = neum_to_diri_bdry_c_B;
            neum_to_diri_param_c = neum_to_diri_param_c_B;
            oldrhs_fem_fct0_c = oldrhs_fem_fct0_c_B;
            oldrhs_fem_fct1_c = oldrhs_fem_fct1_c_B;
            lump_mass_c = lump_mass_c_B;
            matrix_D_Entries_c = matrix_D_Entries_c_B;
            tilde_u_c = tilde_u_c_B;

            BoundaryConditions_Scalar[0] =  BoundCondition_c_B;
            BoundValues_Scalar[0] = BoundValue_c_B;

            for (i=low_scalar;i<mg_level;i++)
            {
              MatricesA_c[i] = MatricesA_c_B[i];
              if(TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)
              {
                MatricesK_c[i] = MatricesK_c_B[i];
                if (TDatabase::ParamDB->SOLD_TYPE)
                  MatricesS_c[i] = MatricesS_c_B[i];
              }
              MatricesM_c[i] = MatricesM_c_B[i];
              N_Array_c[i] = N_Array_c_B[i];
              ConcentrationSpaces[i] = ConcentrationSpaces_c_B[i];
              ConcentrationSpaces_other[i] = ConcentrationSpaces_c_A[i];
              SolArray_other[i] = SolArray_c_A[i];
              if ((TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)&&
                (TDatabase::ParamDB->SOLD_TYPE))
              {
                SolArray[i] = SolArray_c_B[i];
                SolArray_old[i] = SolArray_c_B_old[i];
                SolArray_other_old[i] = SolArray_c_A_old[i];
              }
              if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
                (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
                MatricesK_c[i] = MatricesK_c_B[i];
            }
            nonlin_ite++;
            substance = 1;

            OutPut("********* Computing substance B **********"<< endl);

            break;
            // substance B was computed, now compute substance A
          case 1:
            N_Unknowns_c = N_Unknowns_c_A;
            N_Active_c = N_Active_c_A;
            sol_c = sol_c_A;
            oldsol_c = oldsol_c_A;
            current_sol_c = sol_c_A;
            itmethod_c = itmethod_c_A;
            itmethod_sol_c = itmethod_sol_c_A;
            RhsArray_c = RhsArray_c_A;
            oldrhs_c = oldrhs_c_A;
            rhs_c_complete = rhs_c_complete_A;
            Coeff_c = Coefficients[1];
            N_neum_to_diri_c = N_neum_to_diri_c_A;
            neum_to_diri_c = neum_to_diri_c_A;
            neum_to_diri_bdry_c = neum_to_diri_bdry_c_A;
            neum_to_diri_param_c = neum_to_diri_param_c_A;
            oldrhs_fem_fct0_c = oldrhs_fem_fct0_c_A;
            oldrhs_fem_fct1_c = oldrhs_fem_fct1_c_A;
            lump_mass_c = lump_mass_c_A;
            matrix_D_Entries_c = matrix_D_Entries_c_A;
            tilde_u_c = tilde_u_c_A;

            BoundaryConditions_Scalar[0] =  BoundCondition_c_A;
            BoundValues_Scalar[0] = BoundValue_c_A;

            for (i=low_scalar;i<mg_level;i++)
            {
              MatricesA_c[i] = MatricesA_c_A[i];
              if(TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)
              {
                MatricesK_c[i] = MatricesK_c_A[i];
                if (TDatabase::ParamDB->SOLD_TYPE)
                  MatricesS_c[i] = MatricesS_c_A[i];
              }
              MatricesM_c[i] = MatricesM_c_A[i];
              N_Array_c[i] = N_Array_c_A[i];
              ConcentrationSpaces[i] = ConcentrationSpaces_c_A[i];
              ConcentrationSpaces_other[i] = ConcentrationSpaces_c_B[i];
              SolArray_other[i] = SolArray_c_B[i];
              if ((TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)&&
                (TDatabase::ParamDB->SOLD_TYPE))
              {
                SolArray[i] = SolArray_c_A[i];
                SolArray_old[i] = SolArray_c_A_old[i];
                SolArray_other_old[i] = SolArray_c_B_old[i];
              }
              if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
                (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
                MatricesK_c[i] = MatricesK_c_A[i];
            }
            nonlin_resid = 0;
            substance = 0;
            OutPut("******** Computing substance A ********"<< endl);
            break;
        }

        // rhs only in first step of the nonlinear iteration
        if (nonlin_ite-substance==0)
        {
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

          ferhs[0] = ConcentrationSpaces[mg_level-1];
          N_FESpaces = 3;
          fesp[0] = ConcentrationSpaces[mg_level-1];
          fesp[1] = ConcentrationSpaces_other[mg_level-1];
          fesp[2] = USpaces[mg_level-1];
          fefct[0] = SolArray_other[mg_level-1];
          fefct[1] = U1Array[mg_level-1];
          fefct[2] = U2Array[mg_level-1];

          aux =  new TAuxParam2D(TimeCDParamsBulkN_FESpaces,
            TimeCDParamsBulkN_Fct,
            TimeCDParamsBulkN_ParamFct,
            TimeCDParamsBulkN_FEValues,
            fesp+1, fefct,
            TimeCDParamsBulkFct,
            TimeCDParamsBulkFEFctIndex,
            TimeCDParamsBulkFEMultiIndex,
            TimeCDParamsBulkN_Params,
            TimeCDParamsBulkBeginParam);

          switch(TDatabase::ParamDB->BULK_REACTION_DISC)
          {
            case UPWIND:
              DiscreteForm = DiscreteFormRhs_Bulk;
              break;
            case SDFEM:
              DiscreteForm = DiscreteFormRhs_SUPG_Bulk;
              break;
            case GALERKIN:
              DiscreteForm = DiscreteFormRhs_Galerkin_Bulk;
              break;
            default:
              OutPut("BULK_REACTION_DISC " << TDatabase::ParamDB->BULK_REACTION_DISC
                <<" not available !!!"<<endl);
              exit(4711);
          }

          if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
            (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
	  {
	      BoundValues_Scalar[1] = BoundValues_Scalar[0];
	      BoundValues_Scalar[0] = BoundValue_FEM_FCT;
	  }

          Assemble2D(N_FESpaces, fesp,
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
          {
	      BoundValues_Scalar[0] = BoundValues_Scalar[1];
	      memcpy(oldrhs_fem_fct1_c, rhs_c, N_Unknowns_c*SizeOfDouble);
          }

          // add rhs from current sub time step to rhs array B_c
          Daxpy(N_Active_c, tau*TDatabase::TimeDB->THETA4, rhs_c, B_c);
          // save rhs for next time step
          memcpy(oldrhs_c,rhs_c, N_Unknowns_c*SizeOfDouble);
        }                                         // end nonlin_ite == 0

        // assembling of matrices on all levels
        for(i=low_scalar;i<mg_level;i++)
        {
          // set parameters
          N_FESpaces = 3;
          fesp[0] = ConcentrationSpaces[i];
          fesp[1] = ConcentrationSpaces_other[i];
          fesp[2] = USpaces[i];
          fefct[0] = SolArray_other[i];
          fefct[1] = U1Array[i];
          fefct[2] = U2Array[i];

          // the case of not using a SOLD method
          if (!((TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)&&
            (TDatabase::ParamDB->SOLD_TYPE)))
          {
            aux =  new TAuxParam2D(TimeCDParamsBulkN_FESpaces,
              TimeCDParamsBulkN_Fct,
              TimeCDParamsBulkN_ParamFct,
              TimeCDParamsBulkN_FEValues,
              fesp+1, fefct,
              TimeCDParamsBulkFct,
              TimeCDParamsBulkFEFctIndex,
              TimeCDParamsBulkFEMultiIndex,
              TimeCDParamsBulkN_Params,
              TimeCDParamsBulkBeginParam);
          }
          else
          {
            fefct[3] = SolArray_old[i];
            fefct[4] = SolArray_other_old[i];
            fefct[5] = SolArray[i];
            aux =  new TAuxParam2D(TimeCDParamsBulk_SOLDN_FESpaces,
              TimeCDParamsBulk_SOLDN_Fct,
              TimeCDParamsBulk_SOLDN_ParamFct,
              TimeCDParamsBulk_SOLDN_FEValues,
              fesp, fefct,
              TimeCDParamsBulk_SOLDFct,
              TimeCDParamsBulk_SOLDFEFctIndex,
              TimeCDParamsBulk_SOLDFEMultiIndex,
              TimeCDParamsBulk_SOLDN_Params,
              TimeCDParamsBulk_SOLDBeginParam);
          }

          N_SquareMatrices = 1;
          SQMATRICES[0] = MatricesA_c[i];
          SQMATRICES[0]->Reset();
          switch(TDatabase::ParamDB->BULK_REACTION_DISC)
          {
            case UPWIND:
              DiscreteForm = DiscreteFormMatricesA_Bulk;
              //SQMATRICES[1] = MatricesA_c[i];
              break;
            case SDFEM:
              DiscreteForm = DiscreteFormMatricesA_SUPG_Bulk;
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
              DiscreteForm = DiscreteFormMatricesA_Galerkin_Bulk;
              break;
          }

          Assemble2D(N_FESpaces, fesp,
            N_SquareMatrices, SQMATRICES,
            0, NULL,
            0, NULL, NULL,
            DiscreteForm,
            BoundaryConditions_Scalar,
            BoundValues_Scalar,
            aux);

          if(TDatabase::ParamDB->BULK_REACTION_DISC  == UPWIND)
          {
            // RHSs[0] will not be used in this routine
            UpwindForConvDiff(Coeff_c, SQMATRICES[0],RHSs[0],
              fesp[0],DiscreteForm,
              U1Array[i], U2Array[i], 1);
          }
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

        if (!((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
          (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)))
        {
          // update rhs by Laplacian and convective term from previous
          // time step
          // scaled by current sub time step length and theta2
          // currently : M := M + gamma A
          // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2) A

          // complete rhs with data from the old time step
          // including mass term
          // no contribution from SOLD method
          if (nonlin_ite-substance==0)
          {
            for(i=low_scalar;i<mg_level;i++)
            {
              MatAdd(MatricesM_c[i], MatricesA_c[i],-gamma_c - tau*TDatabase::TimeDB->THETA2);
            }

            // set current factor of steady state matrix
            gamma_c = -tau*TDatabase::TimeDB->THETA2;

            MatVectActive(MatricesM_c[mg_level-1], oldsol_c, defect_c);
            Daxpy(N_Active_c, 1, defect_c, B_c);
            // contribution of SUPG term from time derivative
            if(TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)
            {
              MatVectActive(MatricesK_c[mg_level-1], oldsol_c, defect_c);
              Daxpy(N_Active_c, 1, defect_c, B_c);
            }
            // set Dirichlet values
            // RHSs[0] still available from assembling
            memcpy(B_c+N_Active_c, RHSs[0]+N_Active_c, (N_Unknowns_c-N_Active_c)*SizeOfDouble);

            // save rhs B_c for the following iterations
            memcpy(rhs_c_complete, B_c, N_Unknowns_c*SizeOfDouble);
          }
          else
          {
            // copy the saved rhs to the working vector
            memcpy(B_c, rhs_c_complete, N_Unknowns_c*SizeOfDouble);
          }
          // copy Dirichlet values from rhs into sol
          memcpy(sol_c+N_Active_c, B_c+N_Active_c, (N_Unknowns_c-N_Active_c)*SizeOfDouble);
          // M = M + (-gamma + tau*TDatabase::TimeDB->THETA1) A
          for(i=low_scalar;i<mg_level;i++)
          {
            MatAdd(MatricesM_c[i], MatricesA_c[i],-gamma_c + tau*TDatabase::TimeDB->THETA1);
          }
          // set current factor of steady state matrix
          gamma_c = tau*TDatabase::TimeDB->THETA1;
          // contribution of SUPG term from time derivative
          if(TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)
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
          //if(j>0)
          //  only_first_time = 0;
	  // enforce Crank-Nicolson
	  theta1 = TDatabase::TimeDB->THETA1;
	  theta2 = TDatabase::TimeDB->THETA2;
	  theta3 = TDatabase::TimeDB->THETA3;
	  theta4 = TDatabase::TimeDB->THETA4;
	  TDatabase::TimeDB->THETA1 = 0.5;
	  TDatabase::TimeDB->THETA2 = 0.5;
	  TDatabase::TimeDB->THETA3 = 0.5;
	  TDatabase::TimeDB->THETA4 = 0.5;
          FEM_FCT_ForConvDiff(MatricesK_c[mg_level-1], MatricesA_c[mg_level-1],
            N_Unknowns_c, N_Active_c,
            lump_mass_c, matrix_D_Entries_c,
            sol_c, oldsol_c,
            B_c, RhsArray_c[mg_level-1], oldrhs_fem_fct0_c, tilde_u_c,
            N_neum_to_diri_c, neum_to_diri_c,
            neum_to_diri_bdry_c, neum_to_diri_param_c,
            only_first_time,
            BoundValues_Scalar[0],NULL);
	  TDatabase::TimeDB->THETA1 = theta1;
	  TDatabase::TimeDB->THETA2 = theta2;
	  TDatabase::TimeDB->THETA3 = theta3;
	  TDatabase::TimeDB->THETA4 = theta4;

          SQMATRICES[0] = MatricesM_c[mg_level-1];
          // only first iteration
          //if (j==0)
          //{
          MatricesM_c[mg_level-1]->Reset();
          // system matrix for FEM-FCT   M_lump + theta1*tau*A
          // A = Galerkin + D

          FEM_FCT_SystemMatrix(MatricesM_c[mg_level-1], MatricesA_c[mg_level-1],
            lump_mass_c, N_Unknowns_c);
          //}

          // set Diriclet nodes
          if (N_neum_to_diri_c)
            SetDirichletNodesFromNeumannNodes(SQMATRICES, B_c, sol_c,
              N_neum_to_diri_c, neum_to_diri_c,
              neum_to_diri_bdry_c, neum_to_diri_param_c,
              BoundValues_Scalar[0]);
        }

        //======================================================================
        // solution of linear system
        //======================================================================

        memset(defect_c, 0, N_Unknowns_c*SizeOfDouble);
        SQMATRICES[0] = MatricesM_c[mg_level-1];
        // compute defect
        DefectScalar(sqmatrices,NULL,sol_c,B_c,defect_c);
        residual =  Ddot(N_Unknowns_c, defect_c, defect_c);
        residual = sqrt(residual);

        nonlin_resid += residual;
        if (substance==0)
          OutPut("nonlinear iteration " << nonlin_ite <<
            " residual " <<  nonlin_resid << endl;);
        OutPut("initial residual ");
        OutPut(setw(14) << residual << endl);

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
            DirectSolver(SQMATRICES[0], B_c, itmethod_sol_c);
            Dsum(N_Unknowns_c,1-TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SCALAR,
              TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SCALAR,
              sol_c, itmethod_sol_c, sol_c);
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
        if(TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)
        {
          for(i=low_scalar;i<mg_level;i++)
          {
            MatAdd(MatricesM_c[i], MatricesK_c[i], -1);
            // contribution from SOLD method only to left hand side
            if (TDatabase::ParamDB->SOLD_TYPE)
            {
              MatAdd(MatricesM_c[i], MatricesS_c[i],-1.0);
            }
          }

        }
        // restore mass matrices by subtracting the A-matrices
        for(i=low_scalar;i<mg_level;i++)
        {
          MatAdd(MatricesM_c[i], MatricesA_c[i], -gamma_c);
        }
        // set current factor of steady state matrix
        gamma_c = 0;

      }                                           // nonlinear loop
      OutPut(TDatabase::TimeDB->CURRENTTIME << " NONLIN ITE: " << nonlin_ite <<
        " (LINITE: " << N_LinIter << "), residual: "
        << nonlin_resid <<endl);
      nonlin_resid = 1e8;

      max_c[0] = 0;
      max_c[1] = 0;
      for (i=0;i<N_Unknowns_c;i++)
      {
        if (sol_c_A[i] < 0)
        {
          if (-sol_c_A[i] > max_c[0])
            max_c[0] = -sol_c_A[i];
          if (!((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
            (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)))
            sol_c_A[i] = 0;
        }
        if (sol_c_A[i] > 1)
        {
          if (sol_c_A[i] -1 >  max_c[1])
            max_c[1] = sol_c_A[i] -1;
          if (!((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
            (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)))
            sol_c_A[i] = 1;
        }
      }
      max_c[2] = 0;
      max_c[3] = 0;
      for (i=0;i<N_Unknowns_c;i++)
      {
        if (sol_c_B[i] < 0)
        {
          if (-sol_c_B[i] > max_c[2])
            max_c[2] = -sol_c_B[i];
          if (!((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
            (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)))
            sol_c_B[i] = 0;
        }
        if (sol_c_B[i] > 1)
        {
          if (sol_c_B[i] -1 >  max_c[3])
            max_c[3] = sol_c_B[i] -1;
          if (!((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
            (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)))
            sol_c_B[i] = 1;
        }
      }

      // compute change in the solution to the previous time step
      Daxpy(N_Unknowns_c,-1,sol_c_A,oldsol_c_A);
      Daxpy(N_Unknowns_c,-1,sol_c_B,oldsol_c_B);
      residual =  Ddot(N_Unknowns_c, oldsol_c_A,oldsol_c_A);
      residual +=  Ddot(N_Unknowns_c, oldsol_c_B,oldsol_c_B);
      //don't compute the residual after substance C
      residual = sqrt(residual);

      if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
        (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
      {
        memcpy(oldrhs_fem_fct0_c_A, oldrhs_fem_fct1_c_A, N_Unknowns_c_A*SizeOfDouble);
        memcpy(oldrhs_fem_fct0_c_B, oldrhs_fem_fct1_c_B, N_Unknowns_c_B*SizeOfDouble);
      }

      // *********************************************************
      // substances A and B were computed, now compute substance C
      // *********************************************************
      // step 1: compute c_C using the current c_A, c_B
      N_Unknowns_c = N_Unknowns_c_C;
      N_Active_c = N_Active_c_C;
      sol_c = sol_c_C;
      oldsol_c = oldsol_c_C;
      current_sol_c = sol_c_C;
      itmethod_c = itmethod_c_C;
      itmethod_sol_c = itmethod_sol_c_C;
      RhsArray_c = RhsArray_c_C;
      oldrhs_c = oldrhs_c_C;
      Coeff_c = Coefficients[2];
      N_neum_to_diri_c = N_neum_to_diri_c_C;
      neum_to_diri_c = neum_to_diri_c_C;
      neum_to_diri_bdry_c = neum_to_diri_bdry_c_C;
      neum_to_diri_param_c = neum_to_diri_param_c_C;
      oldrhs_fem_fct0_c = oldrhs_fem_fct0_c_C;
      oldrhs_fem_fct1_c = oldrhs_fem_fct1_c_C;
      lump_mass_c = lump_mass_c_C;
      matrix_D_Entries_c = matrix_D_Entries_c_C;
      tilde_u_c = tilde_u_c_C;
      BoundaryConditions_Scalar[0] =  BoundCondition_c_C;
      BoundValues_Scalar[0] = BoundValue_c_C;

      for (i=low_scalar;i<mg_level;i++)
      {
        MatricesA_c[i] = MatricesA_c_C[i];
        if(TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)
        {
          MatricesK_c[i] = MatricesK_c_C[i];
          if (TDatabase::ParamDB->SOLD_TYPE)
            MatricesS_c[i] = MatricesS_c_C[i];
        }
        MatricesM_c[i] = MatricesM_c_C[i];
        N_Array_c[i] = N_Array_c_C[i];
        ConcentrationSpaces[i] = ConcentrationSpaces_c_C[i];
        ConcentrationSpaces_other[i] = ConcentrationSpaces_c_B[i];
        //SolArray_other[i] = SolArray_c_B[i];
        if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
          (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
          MatricesK_c[i] = MatricesK_c_C[i];
      }
      OutPut("******** Computing substance C ********"<< endl);

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
      N_FESpaces = 5;
      //for c_C
      fesp_c_C[0] = ConcentrationSpaces[mg_level-1];
      // for C_a
      fesp_c_C[1] = ConcentrationSpaces_c_A[mg_level-1];
      // for C_b
      fesp_c_C[2] = ConcentrationSpaces_c_B[mg_level-1];
      //for the velocity
      fesp_c_C[3] = USpaces[mg_level-1];
      //for c_C
      fesp_c_C[4] = ConcentrationSpaces[mg_level-1];

      fefct_c_C[0] = SolArray_c_A[mg_level-1];
      fefct_c_C[1] = SolArray_c_B[mg_level-1];
      fefct_c_C[2] = U1Array[mg_level-1];
      fefct_c_C[3] = U2Array[mg_level-1];
      fefct_c_C[4] = SolArray_c_C[mg_level-1];

      if (!mom)
      {
        fesp_c_C[5] = IntegralSpaces_c_C[mg_level-1];
        fefct_c_C[5] = IntegralSpaces_c_C_fct[mg_level-1];
      }
      else
      {
        // second moment
        fesp_c_C[5] = MOMSpaces[mg_level-1];
        fefct_c_C[5] = SolArray_mom[mg_level-1]->GetComponent(2);
      }

      aux =  new TAuxParam2D(TimeCDParamsBulk_CcN_FESpaces,
        TimeCDParamsBulk_CcN_Fct,
        TimeCDParamsBulk_CcN_ParamFct,
        TimeCDParamsBulk_CcN_FEValues,
        fesp_c_C+1, fefct_c_C,
        TimeCDParamsBulk_CcFct,
        TimeCDParamsBulk_CcFEFctIndex,
        TimeCDParamsBulk_CcFEMultiIndex,
        TimeCDParamsBulk_CcN_Params,
        TimeCDParamsBulk_CcBeginParam);

      switch(TDatabase::ParamDB->BULK_REACTION_DISC)
      {
        case UPWIND:
          DiscreteForm = DiscreteFormRhs_Bulk_Cc;
          break;
        case SDFEM:
          DiscreteForm = DiscreteFormRhs_SUPG_Bulk_Cc;
          break;
        case GALERKIN:
          DiscreteForm = DiscreteFormRhs_Galerkin_Bulk_Cc;
          break;
      }

      if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
	  (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
      {
	  BoundValues_Scalar[1] = BoundValues_Scalar[0];
	  BoundValues_Scalar[0] = BoundValue_FEM_FCT;
      }
      
      // assembling of rhs
      Assemble2D(N_FESpaces, fesp_c_C,
        0, NULL,
        0, NULL,
        N_Rhs, RHSs, ferhs,
        DiscreteForm,
        BoundaryConditions_Scalar,
        BoundValues_Scalar,
        aux);
      OutPut("rhsnorm " << Ddot(N_Active_c, rhs_c,  rhs_c) << endl);
      delete aux;

      if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
        (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
      {
	  BoundValues_Scalar[0] = BoundValues_Scalar[1];
	  memcpy(oldrhs_fem_fct1_c, rhs_c, N_Unknowns_c*SizeOfDouble);
      }
      // add rhs from current sub time step to rhs array B_c
      Daxpy(N_Active_c, tau*TDatabase::TimeDB->THETA4, rhs_c, B_c);
      // save rhs for next time step
      memcpy(oldrhs_c,rhs_c, N_Unknowns_c*SizeOfDouble);

      // assembling of A
      // for SDFEM: in addition stabilisation matrix K
      for(i=low_scalar;i<mg_level;i++)
      {
        //test space!!!
        ferhs[0] = ConcentrationSpaces[i];
        N_FESpaces = 5;
        //for c_C
        fesp_c_C[0] = ConcentrationSpaces[i];
        // for C_a
        fesp_c_C[1] = ConcentrationSpaces_c_A[i];
        // for C_b
        fesp_c_C[2] = ConcentrationSpaces_c_B[i];
        //for the velocity
        fesp_c_C[3] = USpaces[i];
        //for c_C
        fesp_c_C[4] = ConcentrationSpaces[i];

        fefct_c_C[0] = SolArray_c_A[i];
        fefct_c_C[1] = SolArray_c_B[i];
        fefct_c_C[2] = U1Array[i];
        fefct_c_C[3] = U2Array[i];
        fefct_c_C[4] = SolArray_c_C[i];

        if (!mom)
        {
          fesp_c_C[5] = IntegralSpaces_c_C[mg_level-1];
          fefct_c_C[5] = IntegralSpaces_c_C_fct[mg_level-1];
        }
        else
        {
          // second moment
          fesp_c_C[5] = MOMSpaces[mg_level-1];
          fefct_c_C[5] = SolArray_mom[mg_level-1]->GetComponent(2);
        }

        if (!((TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)&&
          (TDatabase::ParamDB->SOLD_TYPE)))
        {
          aux =  new TAuxParam2D(TimeCDParamsBulk_CcN_FESpaces,
            TimeCDParamsBulk_CcN_Fct,
            TimeCDParamsBulk_CcN_ParamFct,
            TimeCDParamsBulk_CcN_FEValues,
            fesp_c_C+1, fefct_c_C,
            TimeCDParamsBulk_CcFct,
            TimeCDParamsBulk_CcFEFctIndex,
            TimeCDParamsBulk_CcFEMultiIndex,
            TimeCDParamsBulk_CcN_Params,
            TimeCDParamsBulk_CcBeginParam);
        }
        else
        {
          fefct_c_C[6] = SolArray_c_C_old[i];
          fefct_c_C[7] = SolArray_c_A_old[i];
          fefct_c_C[8] = SolArray_c_B_old[i];
          aux =  new TAuxParam2D(TimeCDParamsBulk_SOLD_CcN_FESpaces,
            TimeCDParamsBulk_SOLD_CcN_Fct,
            TimeCDParamsBulk_SOLD_CcN_ParamFct,
            TimeCDParamsBulk_SOLD_CcN_FEValues,
            fesp_c_C+1, fefct_c_C,
            TimeCDParamsBulk_SOLD_CcFct,
            TimeCDParamsBulk_SOLD_CcFEFctIndex,
            TimeCDParamsBulk_SOLD_CcFEMultiIndex,
            TimeCDParamsBulk_SOLD_CcN_Params,
            TimeCDParamsBulk_SOLD_CcBeginParam);
        }

        //======================================================================
        // assembling of matrices
        //======================================================================

        N_SquareMatrices = 1;
        SQMATRICES[0] = MatricesA_c[i];
        SQMATRICES[0]->Reset();
        switch(TDatabase::ParamDB->BULK_REACTION_DISC)
        {
          case UPWIND:
            DiscreteForm = DiscreteFormMatricesA_Bulk_Cc;
            //SQMATRICES[1] =  MatricesA_c[i];
            break;
          case SDFEM:
            DiscreteForm = DiscreteFormMatricesA_SUPG_Bulk_Cc;
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
            DiscreteForm = DiscreteFormMatricesA_Galerkin_Bulk_Cc;
            break;
        }

        Assemble2D(N_FESpaces, fesp_c_C,
          N_SquareMatrices, SQMATRICES,
          0, NULL,
          0, NULL, NULL,
          DiscreteForm,
          BoundaryConditions_Scalar,
          BoundValues_Scalar,
          aux);

        if(TDatabase::ParamDB->BULK_REACTION_DISC  == UPWIND)
        {
          UpwindForConvDiff(Coeff_c, SQMATRICES[0],RHSs[0],
            fesp_c_C[0],DiscreteForm,
            U1Array[i], U2Array[i], 1);
        }
        delete aux;

      }                                           // endfor i

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

      if (!((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
        (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)))
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
        if(TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)
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
        if(TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)
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
        //if(j>0)
        //  only_first_time = 0;

	// enforce Crank-Nicolson
	theta1 = TDatabase::TimeDB->THETA1;
	theta2 = TDatabase::TimeDB->THETA2;
	theta3 = TDatabase::TimeDB->THETA3;
	theta4 = TDatabase::TimeDB->THETA4;
	TDatabase::TimeDB->THETA1 = 0.5;
	TDatabase::TimeDB->THETA2 = 0.5;
	TDatabase::TimeDB->THETA3 = 0.5;
	TDatabase::TimeDB->THETA4 = 0.5;
        FEM_FCT_ForConvDiff(MatricesK_c[mg_level-1], MatricesA_c[mg_level-1],
          N_Unknowns_c, N_Active_c,
          lump_mass_c, matrix_D_Entries_c,
          sol_c, oldsol_c,
          B_c, RhsArray_c[mg_level-1], oldrhs_fem_fct0_c, tilde_u_c,
          N_neum_to_diri_c, neum_to_diri_c,
          neum_to_diri_bdry_c, neum_to_diri_param_c,
          only_first_time,
          BoundValues_Scalar[0],NULL);
	TDatabase::TimeDB->THETA1 = theta1;
	TDatabase::TimeDB->THETA2 = theta2;
	TDatabase::TimeDB->THETA3 = theta3;
	TDatabase::TimeDB->THETA4 = theta4;
	
        SQMATRICES[0] = MatricesM_c[mg_level-1];
        // only first iteration
        //if (j==0)
        //{
        MatricesM_c[mg_level-1]->Reset();
        // system matrix for FEM-FCT   M_lump + theta1*tau*A
        // A = Galerkin + D
        FEM_FCT_SystemMatrix(MatricesM_c[mg_level-1], MatricesA_c[mg_level-1],
          lump_mass_c, N_Unknowns_c);
        //}
       // set Diriclet nodes
        if (N_neum_to_diri_c)
          SetDirichletNodesFromNeumannNodes(SQMATRICES, B_c, sol_c,
            N_neum_to_diri_c, neum_to_diri_c,
            neum_to_diri_bdry_c, neum_to_diri_param_c,
            BoundValues_Scalar[0]);
      }

      //======================================================================
      // solution of linear system
      //======================================================================

      memset(defect_c, 0, N_Unknowns_c*SizeOfDouble);
      SQMATRICES[0] = MatricesM_c[mg_level-1];

      // compute defect
      DefectScalar(sqmatrices,NULL,sol_c,B_c,defect_c);

      residual_Cc =  Ddot(N_Unknowns_c, defect_c, defect_c);
      residual_Cc = sqrt(residual_Cc);

      //nonlin_resid += residual;
      OutPut("initial residual ");
      OutPut(setw(14) << residual_Cc << endl);

      //======================================================================
      // solve linear system
      //======================================================================
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
      }                                           // endswitch SOLVER_TYPE

      //======================================================================
      // end solve linear system
      //======================================================================

      // restore mass matrices by subtracting the K-matrices
      if(TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)
      {
        for(i=low_scalar;i<mg_level;i++)
        {
          MatAdd(MatricesM_c[i], MatricesK_c[i], -1);
          // contribution from SOLD method only to left hand side
          if (TDatabase::ParamDB->SOLD_TYPE)
          {
            MatAdd(MatricesM_c[i], MatricesS_c[i],-1.0);
          }
        }
      }
      // restore mass matrices by subtracting the A-matrices
      for(i=low_scalar;i<mg_level;i++)
      {
        MatAdd(MatricesM_c[i], MatricesA_c[i], -gamma_c);
      }
      // set current factor of steady state matrix
      gamma_c = 0;

      if ((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
        (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
      {
        memcpy(oldrhs_fem_fct0_c_C, oldrhs_fem_fct1_c_C, N_Unknowns_c_C*SizeOfDouble);
      }
      OutPut(TDatabase::TimeDB->CURRENTTIME << " LIN ITE: " << N_LinIter << endl);

      max_c[4] = 0;
      max_c[5] = -1;
      for (i=0;i<N_Unknowns_c;i++)
      {
        if (sol_c_C[i] < 0)
        {
          if (-sol_c_C[i] >  max_c[4])
            max_c[4] = -sol_c_C[i];
          // cut off negative values
          if (!((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
            (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)))
            sol_c_C[i] = 0;
        }
        if (sol_c_C[i] > max_c[5])
          max_c[5] = sol_c_C[i];
        if (!((TDatabase::ParamDB->BULK_REACTION_DISC == GALERKIN) &&
          (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN)))
        {
          if (sol_c_C[i] > TDatabase::ParamDB->BULK_REACTION_C_CUT*TDatabase::ParamDB->BULK_c_C_infty)
            sol_c_C[i] = TDatabase::ParamDB->BULK_REACTION_C_CUT*TDatabase::ParamDB->BULK_c_C_infty;
        }
      }
      OutPut(TDatabase::TimeDB->CURRENTTIME <<
        " -A " << max_c[0] <<  " +A " << max_c[1] <<
        " -B " << max_c[2] <<  " +B " << max_c[3] <<
        " -C " << max_c[4] <<  " maximal_c_C " <<  max_c[5] << endl);
      if (max_c[5] == -1)
      {
        OutPut("indefinite values for c_C computed " << endl);
        exit(4711);
      }

      // compute change in the solution to the previous time step
      Daxpy(N_Unknowns_c_C,-1,sol_c_C,oldsol_c_C);

      // store solutions of this time step
      memcpy(oldsol_c_A,sol_c_A,N_Unknowns_c_A*SizeOfDouble);
      memcpy(oldsol_c_B,sol_c_B,N_Unknowns_c_B*SizeOfDouble);
      memcpy(oldsol_c_C,sol_c_C,N_Unknowns_c_C*SizeOfDouble);

      // stop time stepping scheme
      //OutPut(TDatabase::TimeDB->CURRENTTIME << " change in the solution: "
      //  << residual << endl);

      /*if (residual  <= TDatabase::TimeDB->STEADY_STATE_TOL)
      {
         OutPut("time stepping scheme stopped since stationary state reached"
         << endl);
         TDatabase::TimeDB->CURRENTTIME = end_time+1;
      }*/

      //**************************************************
      // computing population size distribution
      //**************************************************
      if (!mom)
      {
        OutPut("******** Computing f  ********"<< endl);
        OutPut("h_PB " << hmin << " Nx_PB/Ny_PB(+1) " << Nx+1 << " Nz_PB(+1) " << Nz+1 <<
          " dof PB " << Nodes << endl);

        if (TDatabase::ParamDB->BULK_PB_DISC== BULK_FWE_FDM_UPWIND)
        {
          if (m==0)                               //first time step
          {
            //initial condition, since explicit scheme
            memset(f_time_space_approx, 0, Nodes*SizeOfDouble);
          }
          else
          {
            OutPut("===== Begin Forward Euler FD Upwind Method ======="<< endl);
            // ======================> Michael 11.01.2007
            t1 = GetTime();
            Bulk_FWE_FDM_Upwind_3D(coll, u1, u2, c_C, f_time_space_approx,
              Nx, Ny, Nz,
              x_coord, y_coord, z_coord,
              x_min, x_max, y_min, y_max,
              z_min, z_max,
              velo1, velo2, concent_C_array,
              correspond_2dgrid);
            t2 = GetTime();
            OutPut("time for psd " << t2-t1 << endl);
            //OutPut("===== End Forward Euler FD Upwind Method ======="<< endl);

            //OutPut("norm of f " << sqrt(Ddot(Nodes,f_time_space_approx,f_time_space_approx)) << endl);
          }
        }

        if (TDatabase::ParamDB->BULK_PB_DISC== BULK_BWE_FDM_UPWIND)
        {
          OutPut("===== Begin Backward Euler FD Upwind Method ======="<< endl);
          t1 = GetTime();
          Bulk_BWE_FDM_Upwind_3D(coll, u1, u2, c_C, f_time_space_approx,
            correspond_2dgrid, Nx, Ny, Nz,
            x_coord, y_coord, z_coord, mat);

          t2 = GetTime();
          OutPut("time for psd " << t2-t1 << endl);
          //OutPut("===== End Backward Euler FD Upwind Method ======="<< endl);
          //OutPut("norm of f " << sqrt(Ddot(Nodes,f_time_space_approx,f_time_space_approx)) << endl);
        }

        if (TDatabase::ParamDB->BULK_PB_DISC== BULK_BWE_FEM_SUPG)
        {
          OutPut("===== Begin Backward Euler FE-SUPG Method ======="<< endl);
          t1 = GetTime();

          Build_3D_FEM_Matrix_Q1(coll, u1, u2, c_C, f_time_space_approx, f_time_space_approx_old,
            lump_mass_PSD, matrix_D_Entries_PSD,
            correspond_2dgrid, Nx, Ny, Nz,
            x_coord, y_coord, z_coord, mat, matM);

          OutPut("===== Begin Backward Euler FE-SUPG Method ======="<< endl);
          //OutPut("norm of f " << sqrt(Ddot(Nodes,f_time_space_approx,f_time_space_approx)) << endl);
          t2 = GetTime();
          OutPut("time for psd " << t2-t1 << endl);
        }

        if (TDatabase::ParamDB->BULK_PB_DISC== BULK_FEM_FCT)
        {
          OutPut("===== Begin Crank-Nicolson FEM-FCT Method ======="<< endl);
          t1 = GetTime();
          if (first_psd_fem)
          {
            Build_3D_FEM_FCT_MassMatrix_Q1(coll, Nx, Ny, Nz, x_coord, y_coord,
              z_coord, index_test_ansatz,
            //index_test_ansatz_diag,
              matM_cons,
              lump_mass_PSD);
            Compute_Neum_To_Diri_FEM_FCT(Nx, Ny, Nz,
              x_coord, y_coord, z_coord,
              N_neum_to_diri_psd,
              neum_to_diri_psd,
              neum_to_diri_psd_x,
              neum_to_diri_psd_y,
              neum_to_diri_psd_z);
            OutPut("Neum_to_diri " << N_neum_to_diri_psd << endl);
            first_psd_fem = 0;
          }
          Build_3D_FEM_FCT_Matrix_Q1(coll, u1, u2, c_C,
            f_time_space_approx, f_time_space_approx_old,
            lump_mass_PSD, matrix_D_Entries_PSD,
            correspond_2dgrid, Nx, Ny, Nz,
            x_coord, y_coord, z_coord, mat, matM_cons, matM,
            index_test_ansatz,
          //index_test_ansatz_diag,
            N_neum_to_diri_psd,
            neum_to_diri_psd,
            neum_to_diri_psd_x,
            neum_to_diri_psd_y,
            neum_to_diri_psd_z);

          //OutPut("norm of f " << sqrt(Ddot(Nodes,f_time_space_approx,f_time_space_approx)) << endl);
          t2 = GetTime();
          OutPut("time for psd " << t2-t1 << endl);
          //OutPut("===== End Crank-Nicolson FEM-FCT Method ======="<< endl);
        }
#ifdef __BULK_ACAD_TEST__
	// for academic test 
	if(TDatabase::ParamDB->MEASURE_ERRORS)
	{
	    fesp[0] = ConcentrationSpaces_c_C[mg_level-1];
	    
	    aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
	    c_C->GetErrors(Exact_c_C, 3, AllDerivatives, 2, L2H1Errors,
			    BilinearCoeffs_Cc, aux, 1, fesp, errors);
	    delete aux;
	    
	    OutPut("c_C error: time: " << TDatabase::TimeDB->CURRENTTIME);
	    OutPut(" L2: " << errors[0]);
	    OutPut(" H1-semi: " << errors[1] << endl);
	    ErrorPSD(Nx, Ny, Nz,
		     x_coord, y_coord, z_coord,f_time_space_approx);
	}	
#endif
        // compute f at the outflow boundary
        Evalute_f_at_outflow(Nx, Ny, Nz, x_coord, z_layers_coord, f_time_space_approx,
          average_median, average_step);
        // compute integral values which are needed in the next time step
        Integral_For_Particle_Increase_Term(integral_space_c_C, integral_space_c_C_fct,
          Nx, Ny, Nz,
          x_coord, y_coord, z_layers_coord,
          concent_C_array, f_time_space_approx);
      }                                           // end if (!mom)
      else
      {
        // assembling the matrices, they are the same for all moments
        // mass matrix already assembled
        // use the same temporal variables as for the reactions
        N_Unknowns_c = N_Unknowns_mom;
        N_Active_c = N_Active_mom;
        sol_c = sol_mom;
        oldsol_c = oldsol_mom;
        current_sol_c = sol_mom;
        itmethod_c = itmethod_mom;
        itmethod_sol_c = itmethod_sol_mom;
        RhsArray_c = RhsArray_mom;
        oldrhs_c = oldrhs_mom;
        Coeff_c = Coefficients[2];

        BoundaryConditions_Scalar[0] =  BoundCondition_mom;
        BoundValues_Scalar[0] = BoundValue_mom;

        for (i=low_scalar;i<mg_level;i++)
        {
          MatricesA_c[i] = MatricesA_mom[i];
          if(TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)
          {
            MatricesK_c[i] = MatricesK_mom[i];
            if (TDatabase::ParamDB->SOLD_TYPE)
              MatricesS_c[i] = MatricesS_mom[i];
          }
          MatricesM_c[i] = MatricesM_mom[i];
          N_Array_c[i] = N_Array_mom[i];
        }
        gamma_c = 0;
        BoundaryConditions_Scalar[0] =  BoundCondition_mom;
        BoundValues_Scalar[0] = BoundValue_mom;
        TDatabase::ParamDB->INTERNAL_MOMENT = 0;

        // assembling of A
        // for SDFEM: in addition stabilisation matrix K
        for(i=low_scalar;i<mg_level;i++)
        {

          ferhs[0] = MOMSpaces[i];
          N_FESpaces = 3;

          fesp_mom[0] =  USpaces[i];
          fesp_mom[1] =  ConcentrationSpaces_c_C[i];
          fesp_mom[2] =  MOMSpaces[i];

          fefct_mom[0] = U1Array[i];
          fefct_mom[1] = U2Array[i];
          // this if for computing rhs, just dummy here that
          // BilinearCoeff_mom is computed correctly
          fefct_mom[2] = SolArray_c_C[i];
          fefct_mom[3] = SolArray_mom[i];

          if (!((TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)&&
            (TDatabase::ParamDB->SOLD_TYPE)))
          {
            aux =  new TAuxParam2D(TimeCDParamsBulk_momN_FESpaces,
              TimeCDParamsBulk_momN_Fct,
              TimeCDParamsBulk_momN_ParamFct,
              TimeCDParamsBulk_momN_FEValues,
              fesp_mom, fefct_mom,
              TimeCDParamsBulk_momFct,
              TimeCDParamsBulk_momFEFctIndex,
              TimeCDParamsBulk_momFEMultiIndex,
              TimeCDParamsBulk_momN_Params,
              TimeCDParamsBulk_momBeginParam);
          }
          else
          {
            // has to be checked
            //fefct_mom[4] = SolArray_mom_old[i];
            aux =  new TAuxParam2D(TimeCDParamsBulk_SOLD_momN_FESpaces,
              TimeCDParamsBulk_SOLD_momN_Fct,
              TimeCDParamsBulk_SOLD_momN_ParamFct,
              TimeCDParamsBulk_SOLD_momN_FEValues,
              fesp_mom, fefct_mom,
              TimeCDParamsBulk_SOLD_momFct,
              TimeCDParamsBulk_SOLD_momFEFctIndex,
              TimeCDParamsBulk_SOLD_momFEMultiIndex,
              TimeCDParamsBulk_SOLD_momN_Params,
              TimeCDParamsBulk_SOLD_momBeginParam);
          }

          //======================================================================
          // assembling of matrices
          //======================================================================

          N_SquareMatrices = 2;
          SQMATRICES[0] = MatricesA_mom[i];
          SQMATRICES[0]->Reset();
          switch(TDatabase::ParamDB->BULK_REACTION_DISC)
          {
            case UPWIND:
              OutPut("UPWIND not implemented for MOM !!!"<<endl);
              exit(4711);
            case SDFEM:
              DiscreteForm = DiscreteFormMatricesA_SUPG_Bulk_mom;
              SQMATRICES[1] = MatricesK_mom[i];
              SQMATRICES[1]->Reset();
              if (TDatabase::ParamDB->SOLD_TYPE)
              {
                SQMATRICES[2] = MatricesS_mom[i];
                SQMATRICES[2]->Reset();
                N_SquareMatrices = 3;
              }
              break;
          }

          Assemble2D(N_FESpaces, fesp_mom,
            N_SquareMatrices, SQMATRICES,
            0, NULL,
            0, NULL, NULL,
            DiscreteForm,
            BoundaryConditions_Scalar,
            BoundValues_Scalar,
            aux);

          delete aux;

        }                                         // endfor i

        // loop over all moments
        for (j=0;j<N_mom;j++)
        {
          OutPut("computing moment " << j << endl);
          TDatabase::ParamDB->INTERNAL_MOMENT = j;
          // assemble rhs only on finest level
          // use the same temporal variables as for the reactions
          rhs_single_mom = RhsArray_mom[mg_level-1]+j*N_Unknowns_mom;
          // working array for rhs is B, initialize B
          memset(B_mom, 0, N_Unknowns_mom*SizeOfDouble);
          // old rhs multiplied with current subtime step and theta3 on B
          Daxpy(N_Active_mom, tau*TDatabase::TimeDB->THETA3, oldrhs_mom+j*N_Unknowns_mom, B_mom);

          // set parameter for compute term with new rhs
          N_Rhs = 1;
          RHSs[0] = RhsArray_mom[mg_level-1]+j*N_Unknowns_mom;
          memset(RHSs[0],0,N_Array_mom[mg_level-1]*SizeOfDouble);

          ferhs[0] = MOMSpaces[mg_level-1];
          N_FESpaces = 3;

          fesp_mom[0] =  USpaces[mg_level-1];
          fesp_mom[1] =  ConcentrationSpaces_c_C[mg_level-1];
          fesp_mom[2] =  MOMSpaces[mg_level-1];

          // this if for computing matrix, just dummy here that
          // BilinearCoeff_mom is computed correctly
          fefct_mom[0] = U1Array[mg_level-1];
          fefct_mom[1] = U2Array[mg_level-1];
          fefct_mom[2] = SolArray_c_C[mg_level-1];
          if (j>0)
          {
            // previous moment
            fefct_mom[3] = SolArray_mom[mg_level-1]->GetComponent(j-1);
          }
          else
          {
            // previous moment
            fefct_mom[3] = SolArray_mom[mg_level-1]->GetComponent(0);
          }

          aux =  new TAuxParam2D(TimeCDParamsBulk_momN_FESpaces,
            TimeCDParamsBulk_momN_Fct,
            TimeCDParamsBulk_momN_ParamFct,
            TimeCDParamsBulk_momN_FEValues,
            fesp_mom, fefct_mom,
            TimeCDParamsBulk_momFct,
            TimeCDParamsBulk_momFEFctIndex,
            TimeCDParamsBulk_momFEMultiIndex,
            TimeCDParamsBulk_momN_Params,
            TimeCDParamsBulk_momBeginParam);

          switch(TDatabase::ParamDB->BULK_REACTION_DISC)
          {

            case SDFEM:
              DiscreteForm = DiscreteFormRhs_SUPG_Bulk_mom;
              break;
          }

          // assembling of rhs
          Assemble2D(N_FESpaces, fesp_mom,
            0, NULL,
            0, NULL,
            N_Rhs, RHSs, ferhs,
            DiscreteForm,
            BoundaryConditions_Scalar,
            BoundValues_Scalar,
            aux);
          delete aux;
          // add rhs from current sub time step to rhs array B_c
          Daxpy(N_Active_mom, tau*TDatabase::TimeDB->THETA4, rhs_single_mom, B_mom);
          OutPut("rhs " << Ddot(N_Unknowns_mom,rhs_single_mom,rhs_single_mom) << endl);
          // save rhs for next time step
          memcpy(oldrhs_mom+j*N_Unknowns_mom, rhs_single_mom, N_Unknowns_mom*SizeOfDouble);

          // update rhs by Laplacian and convective term from previous
          // time step
          // scaled by current sub time step length and theta2
          // currently : M := M + gamma A
          // M = M + (-gamma - tau*TDatabase::TimeDB->THETA2) A
          for(i=low_scalar;i<mg_level;i++)
          {
            MatAdd(MatricesM_mom[i], MatricesA_mom[i], -gamma_c - tau*TDatabase::TimeDB->THETA2);
          }
          // set current factor of steady state matrix
          gamma_c = -tau*TDatabase::TimeDB->THETA2;

          MatVectActive(MatricesM_mom[mg_level-1], oldsol_mom+j*N_Unknowns_mom, defect_mom);
          Daxpy(N_Active_mom, 1, defect_mom, B_mom);
          // contribution of SUPG term from time derivative
          if(TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)
          {
            MatVectActive(MatricesK_mom[mg_level-1], oldsol_mom+j*N_Unknowns_mom, defect_mom);
            Daxpy(N_Active_mom, 1, defect_mom, B_mom);
          }

          // set Dirichlet values
          // RHSs[0] still available from assembling
          memcpy(B_mom+N_Active_mom, RHSs[0]+N_Active_mom, (N_Unknowns_mom-N_Active_mom)*SizeOfDouble);
          // copy Dirichlet values from rhs into sol
          memcpy(fe_mom->GetComponent(j-1)+N_Active_mom, RHSs[0]+N_Active_mom,
            (N_Unknowns_mom-N_Active_mom)*SizeOfDouble);

          // M = M + (-gamma + tau*TDatabase::TimeDB->THETA1) A
          for(i=low_scalar;i<mg_level;i++)
          {
            MatAdd(MatricesM_mom[i], MatricesA_mom[i], -gamma_c + tau*TDatabase::TimeDB->THETA1);
          }
          // set current factor of steady state matrix
          gamma_c = tau*TDatabase::TimeDB->THETA1;
          // contribution of SUPG term from time derivative
          if(TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)
          {
            for(i=low_scalar;i<mg_level;i++)
            {
              MatAdd(MatricesM_mom[i], MatricesK_mom[i], 1);
              if (TDatabase::ParamDB->SOLD_TYPE)
              {
                MatAdd(MatricesM_mom[i], MatricesS_mom[i],1.0);
              }
            }
          }

          //======================================================================
          // solution of linear system
          //======================================================================

          memset(defect_mom, 0, N_Unknowns_mom*SizeOfDouble);
          SQMATRICES[0] = MatricesM_mom[mg_level-1];

          // compute defect
          DefectScalar(sqmatrices,NULL,sol_mom+j*N_Unknowns_mom,B_mom,defect_mom);

          residual_mom =  Ddot(N_Unknowns_mom, defect_mom, defect_mom);
          residual_mom = sqrt(residual_mom);

          //nonlin_resid += residual;
          OutPut("initial residual ");
          OutPut(setw(14) << residual_mom << endl);

          //======================================================================
          // solve linear system
          //======================================================================

          switch(solver_type_reaction)
          {
            case AMG:
              TDatabase::ParamDB->SC_VERBOSE=0;
              t1 = GetTime();
              Solver(SQMATRICES[0], B_mom, sol_mom+j*N_Unknowns_mom);
              t2 = GetTime();
              solver_time_curr = t2-t1;
              solver_time += solver_time_curr;
              break;

            case DIRECT:
              TDatabase::ParamDB->SC_VERBOSE=0;
              t1 = GetTime();
              DirectSolver(SQMATRICES[0], B_mom, sol_mom+j*N_Unknowns_mom);
              t2 = GetTime();
              solver_time_curr = t2-t1;
              solver_time += solver_time_curr;
              break;

            case GMG:
              t1 = GetTime();
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
              {
                memcpy(itmethod_sol_mom, sol_mom+j*N_Unknowns_mom, N_Unknowns_mom*SizeOfDouble);
                memcpy(itmethod_rhs_mom, B_mom, N_Unknowns_mom*SizeOfDouble);
              }
              N_LinIterCurrIte = itmethod_mom->Iterate(sqmatrices,NULL,itmethod_sol_mom,itmethod_rhs_mom);
              if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
              {
                memcpy(sol_mom+j*N_Unknowns_mom, itmethod_sol_mom, N_Unknowns_mom*SizeOfDouble);
                memcpy(B_mom, itmethod_rhs_mom, N_Unknowns_mom*SizeOfDouble);
              }
              t2 = GetTime();
              solver_time_curr += t2-t1;
              solver_time += solver_time_curr;
              break;
          }                                       // endswitch SOLVER_TYPE

          //======================================================================
          // end solve linear system
          //======================================================================

          // restore mass matrices by subtracting the K-matrices
          if(TDatabase::ParamDB->BULK_REACTION_DISC  == SDFEM)
          {
            for(i=low_scalar;i<mg_level;i++)
            {
              MatAdd(MatricesM_mom[i], MatricesK_mom[i], -1);
              // contribution from SOLD method only to left hand side
              if (TDatabase::ParamDB->SOLD_TYPE)
              {
                MatAdd(MatricesM_mom[i], MatricesS_mom[i],-1.0);
              }
            }
          }
          // restore mass matrices by subtracting the A-matrices
          for(i=low_scalar;i<mg_level;i++)
          {
            MatAdd(MatricesM_mom[i], MatricesA_mom[i], -gamma_c);
          }
          // set current factor of steady state matrix
          gamma_c = 0;

          OutPut(TDatabase::TimeDB->CURRENTTIME << " LIN ITE: " << N_LinIterCurrIte << endl);

        }                                         // end for j (number of moments)
      }
    }                                             // endfor two time discs of adaptive time step control

    if (time_discs==2)
    {
      // compute difference of solutions
      for (i=0;i<N_Unknowns_c;i++)
        sol_c[i]-=frac_step_sol[i];

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
      memcpy(sol_c,frac_step_sol,N_Unknowns_c*SizeOfDouble);
    }                                             // adaptive time step control

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
      delete aux;
    }                                             // endif MEASURE_ERRORS

    if ((TDatabase::ParamDB->WRITE_GRAPE)||(TDatabase::ParamDB->WRITE_GNU)
      ||(TDatabase::ParamDB->WRITE_GMV)||(TDatabase::ParamDB->WRITE_VTK))
    {
      if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      {
        StreamFunction(USpaces[mg_level-1], sol_c, sol_c+N_Uarray[mg_level-1],
          PsiSpaces[LEVELS-1], psi);

        //Prolongate(USpaces[mg_level-1], PsiSpaces[LEVELS-1], sol_c,
        //  uconf->GetValues(), app);

        // Prolongate(USpaces[mg_level-1], PsiSpaces[LEVELS-1],
        //  sol_c+N_Uarray[mg_level-1], uconf->GetValues()+N_V, app);
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

    if(TDatabase::ParamDB->WRITE_VTK)
    {
      if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      {
        os.seekp(std::ios::beg);
        os << VtkBaseName << m << ".vtk" << ends;
        Output->WriteVtk(os.str().c_str());
        os.seekp(std::ios::beg);
        os << VtkBaseName << "psd." << m << ".vtk" << ends;
        write_vtk_psd(Nx, Ny, Nz, x_coord, y_coord, z_coord, f_time_space_approx,os.str().c_str());
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

    if (TDatabase::ParamDB->SAVE_DATA)
    {
      if(m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
      {
        save_sol[0] = UArray[mg_level-1]->GetValues();
        save_sol[1] = SolArray_c_A[mg_level-1]->GetValues();
        save_sol[2] = SolArray_c_B[mg_level-1]->GetValues();
        save_sol[3] = SolArray_c_C[mg_level-1]->GetValues();
        save_sol[4] = f_time_space_approx;
        save_N_Unknowns[0] = 2*UArray[mg_level-1]->GetLength()+PArray[mg_level-1]->GetLength();
        save_N_Unknowns[1] = SolArray_c_A[mg_level-1]->GetLength();
        save_N_Unknowns[2] = SolArray_c_B[mg_level-1]->GetLength();
        save_N_Unknowns[3] = SolArray_c_C[mg_level-1]->GetLength();
        save_N_Unknowns[4] = Nodes;
        SaveData(SaveDataFileName,5,save_sol,save_N_Unknowns);
      }
    }

    // for BMBF project
    /*os.seekp(std::ios::beg);
    os << VtkBaseName << "psd" << ends;
    write_psd(Nx, Ny, Nz, x_coord, y_coord, z_layers_coord, f_time_space_approx,os.str().c_str());
    */
    comp_vort =0;
  }                                               // while

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
    os.seekp(std::ios::beg);
    os << VtkBaseName << "end.psd." << m << ".vtk" << ends;
    write_vtk_psd(Nx, Ny, Nz, x_coord, y_coord, z_coord, f_time_space_approx,os.str().c_str());
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
