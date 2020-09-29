// =======================================================================
//
// Purpose:     main program
//
// Author:      Volker John 2000/08/25
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
#include <VMS.h>
#include <Windtunnel_3d4d.h>

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

//#include "../Examples/TNSE_3D/DrivenCavity3D_Bulk.h"
#include "../Examples/TNSE_3D/windtunnel_fine.h"

// =======================================================================
// start of main programm
// =======================================================================

int main(int argc, char* argv[])
{
  //======================================================================
  // begin of the declaration of the variables
  //======================================================================

  // integer variables
  int i, j, k, l, m, n, ii, ll;
   int CurrentDiscType, comp_vort;
  int LEVELS;
  int Max_It, methods, mg_level, mg_type, mid_sq, last_sq;
  int N, N_, N_Active;
  int n_aux, N_Cells, N_Entries, N_FESpaces, N_L, N_LinIter, N_LinIterCurr, N_LinIterCurrIte;
  int Nodes, nonlinite, nonlin_ite, N_P, N_RectMatrices, N_Rhs, N_SquareMatrices, N_SubSteps;
  int N_U, N_Unknowns;
  int Nx, Ny, Nz, Na;
  int pressure_space_code;
  int ret;
  int time_discs, velocity_space_code, zerostart;

  // initialised integer variables
  int first_psd_fem = 1;
  int mixing_layer_galerkin = 0;
  int N_neum_to_diri_psd = 0;
  int N_Paramters = 1, very_first_time = 0,only_first_time=1;
  

  // integer pointers
  int *col_ptr, *correspond_3dgrid, *N_Parray, *N_Uarray, *RowPtr;
  int *N_Array, *row_ptr;
  int *neum_to_diri, N_neum_to_diri;
  int *neum_to_diri_psd;
  int *index_test_ansatz;
 

  int **downwind;
  int ***mean_value_outflow_indices;

  // double variables
  double end_time, gamma, h;
  double hmin, hmax, impuls_residual;
  double limit;
  double oldresidual, oldtau, omega, p2;
  double real_time, res, residual;
  double solver_time, solver_time_curr;
  double t, t1, t11, t2, t22, t3, t4, tau, tau1, tau2, theta1, theta2, theta3, theta4, total_time;
  double x, y;

  // initialised double variables
  double l_infty_l_2 = 0, l_infty_l_2_time = -4711.0, l_2_h_1u = 0, l_2_l_2u = 0;
  double  olderror_l_2_l_2u = 0, olderror_l_2_h_1u = 0;
 

  // interval limits for the coordinates
  double x_min, x_max, y_min, y_max, z_min, z_max, a_min, a_max;
  // coordinate for the 3d cut -> visulization in paraview

  // declaration and initialisation of the coordinate vectors
  double *x_coord, *y_coord, *z_coord, *a_coord;
  // declaration and initialisation of the vector with the differences between the z layers
  double *a_layers_coord;

  // double pointers
  double *B;
  double *defect,*div, *dmean_velocity;
  double *frac_step_sol, *sol_psd_help, *sol_psd, *h1p;
  double *integral_val, *itmethod_rhs, *itmethod_sol;
  double *l_inf, *l2p, *LESModelRhs;
  double *lump_mass_PSD;
  double *matrix_D_Entries_PSD;
  double *mean_velocity, *mean_velocity_u2, *mean_velocity_u3;
  double *neum_to_diri_x, *neum_to_diri_y, *neum_to_diri_z;
  double *neum_to_diri_psd_x, *neum_to_diri_psd_y, *neum_to_diri_psd_z, *neum_to_diri_psd_a;
  double *newton;
  double *oldrhs;
  double *oldsol;
  double *projection_u1x, *psi;
  double *R_xx, *R_yy, *R_zz, *R_xy, *R_xz, *R_yz, *rhsGL00AuxProblem, *rhs_vms_expl;
  double *rhs;
  double *rms_velocity, *rms_velocity2, *rms_velocity3, *rms_velocity1_type1, *rms_velocity2_type1, *rms_velocity3_type1;
  double *sol, *solGL00AuxProblem, *sol_timestep_m1, *sol_vort_tmp, *startsol;
  double *u1x, *u_conv, *u_uConv;
  double *velo1, *velo2, *velo3, *vms_projection, *vorticity;
  double *x_dof, *y_dof, *z_dof;
  double *size_small_scales,*label_space;
  double cut_coord = 0.25;
  double **RhsArray;
  double ***mean_value_outflow;
  double ***diff_velo_air_drops;
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

  CoeffFct3D *Coefficients[4];

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
  TDiscreteForm3D *DiscreteFormMatricesA_Bulk;
  //TDiscreteForm3D *DiscreteFormMatricesA_Bulk_Cc;
  TDiscreteForm3D *DiscreteFormMatricesA_Galerkin_Bulk;
  //TDiscreteForm3D *DiscreteFormMatricesA_Galerkin_Bulk_Cc;
  TDiscreteForm3D *DiscreteFormMatricesA_SUPG_Bulk;
  // TDiscreteForm3D *DiscreteFormMatricesA_SUPG_Bulk_Cc;
  TDiscreteForm3D *DiscreteFormMatrixAuxProblemU;
  TDiscreteForm3D *DiscreteFormMatrixGL00AuxProblem;
  TDiscreteForm3D *DiscreteFormMatrixMBulk;
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
  TDiscreteForm3D *DiscreteFormRhs_Bulk;
  //TDiscreteForm3D *DiscreteFormRhs_Bulk_Cc;
  TDiscreteForm3D *DiscreteFormRhs_Galerkin_Bulk;
  //TDiscreteForm3D *DiscreteFormRhs_Galerkin_Bulk_Cc;
  TDiscreteForm3D *DiscreteFormRHSAuxProblemU;
  TDiscreteForm3D *DiscreteFormRHSClassicalLES;
  TDiscreteForm3D *DiscreteFormRHSLESModel;
  TDiscreteForm3D *DiscreteFormRHSNewton;
  TDiscreteForm3D *DiscreteFormRHSNewtonNL;
  TDiscreteForm3D *DiscreteFormRhs_SUPG_Bulk;
  // TDiscreteForm3D *DiscreteFormRhs_SUPG_Bulk_Cc;
  TDiscreteForm3D *DiscreteFormSmagorinsky;
  TDiscreteForm3D *DiscreteFormUpwind;
  TDiscreteForm3D *DiscreteFormUpwindNC;
  TDiscreteForm3D *DiscreteFormVMS_Projection;
  TDiscreteForm3D *DiscreteFormVMS_SUPG;
  TDiscreteForm3D *DiscreteFormNLVMS_SUPG;
  TDiscreteForm3D *DiscreteFormRHSSUPG;

  TDomain *Domain = new TDomain();

  TFEDatabase3D *FEDatabase = new TFEDatabase3D();

  TFEFunction3D *Approx, *AuxPArray, *Divergence;
  TFEFunction3D *du11Conv, *du12Conv, *du13Conv, *du22Conv, *du23Conv, *du33Conv;
  TFEFunction3D *GL00AuxProblemSol11, *GL00AuxProblemSol12;
  TFEFunction3D *GL00AuxProblemSol13, *GL00AuxProblemSol22;
  TFEFunction3D *GL00AuxProblemSol23, *GL00AuxProblemSol33;
  TFEFunction3D *old_p, *old_u, *p, *pConv;
  TFEFunction3D *separated_pressure_fe_funct, *separated_pressure_rhs_fe_funct;
  TFEFunction3D *soldiff_fe1,*soldiff_fe2, *StreamFct;
  TFEFunction3D *u1, *u2, *u3, *u1Conv, *u2Conv, *u3Conv, *u4Conv, *u5Conv, *u6Conv;
  TFEFunction3D *vort1, *vort2, *vort3, *Vort_x, *Vort_y, *Vort_z;
  TFEFunction3D *vms_proj_11, *vms_proj_12, *vms_proj_13, *vms_proj_22;
  TFEFunction3D *vms_proj_23, *vms_proj_33, *size_small_scales_fefct, *label_space_fefct;

  TFEFunction3D *fefct[12];

  TFEFunction3D **AuxFEFunctArray;
  TFEFunction3D **du11ConvArray, **du12ConvArray, **du13ConvArray;
  TFEFunction3D **du22ConvArray, **du23ConvArray, **du33ConvArray;
  TFEFunction3D **GL00AuxProblemSol11Array, **GL00AuxProblemSol12Array;
  TFEFunction3D **GL00AuxProblemSol13Array, **GL00AuxProblemSol22Array;
  TFEFunction3D **GL00AuxProblemSol23Array, **GL00AuxProblemSol33Array;
  TFEFunction3D **PArray, **pConvArray;
  //TFEFunction3D **SolArray_other;
  // TFEFunction3D **SolArray_old, **SolArray_other_old, **SolArray;
  TFEFunction3D **u1ConvArray, **u2ConvArray, **u3ConvArray;
  TFEFunction3D **u4ConvArray, **u5ConvArray, **u6ConvArray;
  TFEFunction3D **U1Array, **U2Array, **U3Array;

  TFEVectFunct3D *u, **UArray, *uconf, *duConv, **duConvArray;
  TFEVectFunct3D *uConv, **uConvArray, **AuxFEVectFunctArray;
  TFEVectFunct3D *Vorticity, *vms_projection_fe;
  TFEVectFunct3D *GL00AuxProblemSol, **GL00AuxProblemSolArray;

  TFESpace3D *convolution_space;
  TFESpace3D *old_p_space, *old_u_space;
  TFESpace3D *pressure_separation_space, *pressure_space, *projection_space;
  TFESpace3D *streamfunction_space;
  TFESpace3D *velocity_space, *vorticity_space, *size_small_scales_fesp, *label_space_fesp;

  TFESpace3D *ferhs[6], *fesp[4];

  TFESpace3D **ConcentrationSpaces, **ConcentrationSpaces_other;
  TFESpace3D **duConvSpaces;
  TFESpace3D **pConvSpaces, **ProjectionSpaces, **PsiSpaces, **PSpaces;
  TFESpace3D **uConvSpaces, **USpaces;
  FE3D *fes, *fes1;

  TItMethod *itmethod;
  TItMethod *prec;

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

  TMGLevel3D  *MGLevelGL00AuxProblem;

  TMultiGrid3D *MGGL00AuxProblem;

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

  TSquareMatrix3D **MatricesGL00AuxProblem;
  TSquareMatrix3D **MatricesM11, **MatricesM12, **MatricesM13;
  TSquareMatrix3D **MatricesM21, **MatricesM22, **MatricesM23;
  TSquareMatrix3D **MatricesM31, **MatricesM32, **MatricesM33;
  TSquareMatrix3D **MatricesA, **MatricesK, **MatricesL, **MatricesM;
  TSquareMatrix3D **MatricesA11, **MatricesA12, **MatricesA13;
  TSquareMatrix3D **MatricesA21, **MatricesA22, **MatricesA23;
  TSquareMatrix3D **MatricesA31, **MatricesA32, **MatricesA33;

  TSquareStructure3D *matrix_structure, *sqstructureA, *sqstructureC;
  TSquareStructure3D  *sqstructureL, *sqstructurePressSep;

  TStructure3D *structureB, *structureBT;
  TStructure3D *structure_G, *structure_tilde_G;

  // strings
  char AuxProblemString[] = "AuxProblem";
  // char c_A_String[] = "c_A";
  //char c_B_String[] = "c_B";
  //char c_C_String[] = "c_C";
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

  InitializeDiscreteForms(DiscreteFormGalerkin,
    DiscreteFormUpwind,
    DiscreteFormUpwindNC,
    DiscreteFormSmagorinsky,
    DiscreteFormClassicalLES,
    DiscreteFormGL00Convolution,
    DiscreteFormGL00AuxProblem,
    DiscreteFormVMS_Projection,
    DiscreteFormVMS_SUPG,
    DiscreteFormNLGalerkin,
    DiscreteFormNLUpwind,
    DiscreteFormNLUpwindNC,
    DiscreteFormNLSmagorinsky,
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
    DiscreteFormC,
    DiscreteFormJ,
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

  double xv, yv, zv;
  int N_Vertices;
  coll=Domain->GetCollection(It_Finest, 0);
  for (i=0;i< coll->GetN_Cells();i++)
    {
      cell = coll->GetCell(i);
      // number of vertices
      N_Vertices = cell->GetN_Vertices();
      for ( j=0 ; j<N_Vertices ; j++ )
	{
	  cell->GetVertex(j)->SetClipBoard(0);
	}
    }
  for (i=0;i< coll->GetN_Cells();i++)
    {
      cell = coll->GetCell(i);
      // number of vertices
      N_Vertices = cell->GetN_Vertices();
      for ( j=0 ; j<N_Vertices ; j++ )
	{
	   cell->GetVertex(j)->GetCoords(xv, yv, zv);
	   if ((zv >=0)&&( cell->GetVertex(j)->GetClipBoard()!=4711))
	     {
	       cell->GetVertex(j)->SetClipBoard(4711);
	       cell->GetVertex(j)->SetCoords(xv, yv, zv-0.18);
	     }
	}
    }

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
  // c_A, c_B
  //Coefficients[1] = BilinearCoeffs;
  // c_C
  //Coefficients[2] = BilinearCoeffs_Cc;
  //Coefficients[3] = NoCoeffs;

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
  // definitions for population balance equation
  //======================================================================
  a_min = TDatabase::ParamDB->WINDTUNNEL_R_MIN/TDatabase::ParamDB->WINDTUNNEL_R_INFTY;
  a_max = 1.0;
  cut_coord = 0;

#ifdef  __WINDTUNNEL__
  srand48(time(NULL));
  ReadExperimentalBoundaryConditions();
 #endif
  
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

    if ((i==mg_level-1)&&(TDatabase::TimeDB->EXTRAPOLATE_VELOCITY))
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

      switch(TDatabase::ParamDB->NSTYPE)          //in meinem Fall 4
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
            B, sol, n_aux, alpha, velocity_space_code, pressure_space_code,coll,downwind[i]);
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

    // set discrete forms
    // read initial solution of finest level from grape file to
    // generate the correct matrix for the NSE
    if ((i==mg_level-1)&&(TDatabase::ParamDB->READ_GRAPE_FILE))
    {
      // only velocity and pressure is read
      AuxFEFunctArray = new TFEFunction3D*[1];
      AuxFEFunctArray[0] = PArray[mg_level-1];
      AuxFEVectFunctArray = new TFEVectFunct3D*[1];
      AuxFEVectFunctArray[0] = UArray[mg_level-1];
      ReadGrapeFile3D(ReadGrapeBaseName, 1 , 1, AuxFEFunctArray,AuxFEVectFunctArray);
    }

    // read initial solution of finest level from grape file
    if ((i==mg_level-1)&&(TDatabase::ParamDB->READ_DATA))
    {
      save_sol[0] = sol;
      save_N_Unknowns[0] = N_Unknowns;
      ReadData(ReadDataFileName,1,save_sol,save_N_Unknowns);
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
            //OutPut("fine"<<endl);
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
           // OutPut("coarse"<<endl);
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

    //========================================================================
    //general definitions
    //
    //========================================================================

    // prepare output
    if (i==mg_level-1)
    {
      if ((TDatabase::ParamDB->WRITE_GRAPE)||(TDatabase::ParamDB->WRITE_GMV)||(TDatabase::ParamDB->WRITE_VTK))
      {
        // output velocity, pressure, c_A, c_B, c_C
        Output = new TOutput3D(5, 4, 1, 2, Domain);
        Output->AddFEVectFunct(u);
        Output->AddFEFunction(p);
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
      AuxFEVectFunctArray = new TFEVectFunct3D*[1];
      AuxFEVectFunctArray[0] = UArray[mg_level-1];
      ReadGrapeFile3D(ReadGrapeBaseName, 1 , 1, AuxFEFunctArray,AuxFEVectFunctArray);
      if (TDatabase::TimeDB->RESET_CURRENTTIME > 0)
      {
        TDatabase::TimeDB->CURRENTTIME = TDatabase::TimeDB->RESET_CURRENTTIME_STARTTIME;
        OutPut("start time reset to " << TDatabase::TimeDB->CURRENTTIME << endl);
      }
      // read initial solution of finest level from grape file
      if ((i==mg_level-1)&&(TDatabase::ParamDB->READ_DATA))
      {
        save_sol[0] = sol;
        save_N_Unknowns[0] = N_Unknowns;
        ReadData(ReadDataFileName,1,save_sol,save_N_Unknowns);
      }
    }
  }

  t4 =  GetTime();
  total_time += t4 - t3;
  t3 = t4;

  //======================================================================
  // end of space cycle, finest grid reached
  // everything happens on the same grid
  //======================================================================
  // definitions for Navier-Stokes equations
  // copy sol for extrapolation after time step
  if (TDatabase::TimeDB->EXTRAPOLATE_VELOCITY)
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
      //os << VtkBaseName << "psd." << 0 << ".vtk" << ends;
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

  // definitions for turbulence models
  coll->GetHminHmax(&hmin,&hmax);
  OutPut("h_min : " << hmin << " h_max : " << hmax << endl);

  /* #ifdef __WINDTUNNEL__
   CheckNeumannNodesForVelocity(coll, USpaces[mg_level-1], N_neum_to_diri, neum_to_diri);
   #endif*/

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

  OutPut("grid generation for population balance equation " << endl);

  grid_generator_4d(coll, Nx, Ny, Nz,
    x_min, x_max, y_min, y_max, z_min, z_max,
    a_min, a_max, Na,
    x_coord, y_coord, z_coord, a_coord,
    a_layers_coord);

  OutPut("grid for population balance: " << Nx << " x " << Ny << " x "
    << Nz << " x " << Na << endl);
  OutPut("size of 3D domain: ["<<x_min<<","<<x_max<<"] x ["
    <<y_min<<","<<y_max<<"] x [" <<z_min<<","<<z_max<<"]" <<endl);

  // total number of nodes or grid points
  Nodes = (Nx+1)*(Ny+1)*(Nz+1)*(Na+1);

  sol_psd_help = new double[Nodes];
  memset(sol_psd_help, 0, Nodes*SizeOfDouble);
  sol_psd = new double[Nodes];
  memset(sol_psd, 0, Nodes*SizeOfDouble);

  velo1 = new double[(Nx+1)*(Ny+1)*(Nz+1)];
  memset(velo1, 0, (Nx+1)*(Ny+1)*(Nz+1)*SizeOfDouble);
  velo1[0] = -4711;
  velo2 = new double[(Nx+1)*(Ny+1)*(Nz+1)];
  memset(velo2, 0, (Nx+1)*(Ny+1)*(Nz+1)*SizeOfDouble);
  velo3 = new double[(Nx+1)*(Ny+1)*(Nz+1)];
  memset(velo3, 0, (Nx+1)*(Ny+1)*(Nz+1)*SizeOfDouble);
 
 //mean values at the end of the wind tunnel
  alloc_cubix_int(&mean_value_outflow_indices, Ny+1, Nz+1, Na+1);
  alloc_cubix(&mean_value_outflow, Ny+1, Nz+1, Na+1);
  alloc_cubix(&diff_velo_air_drops, (TDatabase::ParamDB->WINDTUNNEL_LAYER_NUMBER_X)+1, Ny+2, Nz+2);
 
#ifdef  __WINDTUNNEL__
  ReadExperimentalBoundaryConditionsDrops();
  ReadExperimentalVelocityDrops(diff_velo_air_drops);
  //accumulate_bound_condition_drops_neu(Ny, Nz, Na, a_layers_coord);
#endif
/*for(i=0;i<=3;i++)
    for(j=0;j<=Ny+1;j++)
        for(n=0;n<=Nz+1;n++)
            {
             cout << diff_velo_air_drops[i][j][n]<< " ";
             }
*/
  if ((TDatabase::ParamDB->BULK_PB_DISC==BULK_BWE_FEM_SUPG) ||
    (TDatabase::ParamDB->BULK_PB_DISC==WINDTUNNEL_BWE_FDM_UPWIND) ||
    (TDatabase::ParamDB->BULK_PB_DISC==WINDTUNNEL_FEM_FCT))
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

    if((TDatabase::ParamDB->BULK_PB_DISC_STAB == GALERKIN) &&
      (TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE == FEM_FCT_LIN))
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

  if (TDatabase::ParamDB->BULK_PB_DISC==WINDTUNNEL_FWE_FDM_UPWIND)
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
              break;
          }                                       // endswitch
        }
        // set current factor of steady state matrix
        gamma = tau*theta1;

        //========================================================================
        // end assembling of system matrix
        //========================================================================

        real_time = TDatabase::TimeDB->CURRENTTIME * TDatabase::ParamDB->WINDTUNNEL_U_INFTY/TDatabase::ParamDB->WINDTUNNEL_L_INFTY ;
        OutPut(endl << "CURRENT TIME: ");
        OutPut(TDatabase::TimeDB->CURRENTTIME);
        OutPut(" (real time: " << real_time << " s)" << endl);
        OutPut(" MEMORY: " << setw(10) << GetMemory()/(1048576.0));
        OutPut(" MB" << endl);

        //======================================================================
        // nonlinear loop
        //======================================================================

        N_LinIterCurr = 0;
        solver_time_curr = 0;

        // solve nonlinear equation
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

        //======================================================================
        // end of nonlinear loop
        //======================================================================
        if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
        {
          IntoL20FEFunction3D(PArray[mg_level-1]->GetValues(),
            N_Parray[mg_level-1], PSpaces[mg_level-1]);
        }

        // **************************************************
        // computing population size distribution
        // **************************************************

        if ( TDatabase::TimeDB->CURRENTTIME > TDatabase::TimeDB->T1)
        {
          OutPut("******** Computing f  ********"<< endl);
          if ( m==1 )
            OutPut("h_PB " << hmin << " Nx_PB " << Nx+1 << " Ny_PB " << Ny+1
              << " Nz_PB " << Nz+1 << " Na_PB " << Na+1 << " dof PB " << Nodes << endl);

          if (TDatabase::ParamDB->BULK_PB_DISC== WINDTUNNEL_FWE_FDM_UPWIND)
          {
            //first time step
            if ( m==1 )
            {
              //initial condition, since explicit scheme
              memset(sol_psd, 0, Nodes*SizeOfDouble);
            }
            else
            {
              OutPut("===== Begin Forward Euler FD Upwind Method ======="<< endl);

              Windtunnel_FWE_FDM_Upwind_4D(coll, u1, u2, u3, sol_psd,
                Nx, Ny, Nz, Na,
                x_coord, y_coord, z_coord, a_coord,
                x_min, x_max, y_min, y_max, z_min, z_max, a_min, a_max,
                velo1, velo2, velo3, correspond_3dgrid, diff_velo_air_drops);

              OutPut("===== End Forward Euler FD Upwind Method ======="<< endl);
            }
          }

          if (TDatabase::ParamDB->BULK_PB_DISC== WINDTUNNEL_BWE_FDM_UPWIND)
          {
            OutPut("===== Begin Backward Euler FD Upwind Method ======="<< endl);

            Windtunnel_BWE_FDM_Upwind_4D(coll, u1, u2, u3, sol_psd,
              correspond_3dgrid, Nx, Ny, Nz, Na,
              x_min, x_max, y_min, y_max, z_min, z_max, a_min, a_max,
              x_coord, y_coord, z_coord, a_coord, diff_velo_air_drops, mat);

            OutPut("===== End Backward Euler FD Upwind Method ======="<< endl);
          }

          if (TDatabase::ParamDB->BULK_PB_DISC== WINDTUNNEL_FEM_FCT)
          {
            OutPut("===== Begin Crank-Nicolson FEM-FCT Method ======="<< endl);

            if (first_psd_fem)
            {
              Build_4D_FEM_FCT_MassMatrix_Q1_Windtunnel(coll, Nx, Ny, Nz, Na, x_coord, y_coord,
                z_coord, a_coord, index_test_ansatz, matM_cons,
                lump_mass_PSD);

	      Compute_Neum_To_Diri_FEM_FCT_Windtunnel(Nx, Ny, Nz, Na,
                x_coord, y_coord, z_coord, a_coord,
                N_neum_to_diri_psd,
                neum_to_diri_psd,
                neum_to_diri_psd_x,
                neum_to_diri_psd_y,
                neum_to_diri_psd_z,
                neum_to_diri_psd_a);
		OutPut("Neum_to_diri_psd " << N_neum_to_diri_psd << endl);
              first_psd_fem = 0;
            }

            FEM_FCT_Matrix_Q1_4D_Windtunnel(coll, u1, u2, u3, diff_velo_air_drops, sol_psd_help, sol_psd,
              lump_mass_PSD, matrix_D_Entries_PSD,
              correspond_3dgrid, Nx, Ny, Nz, Na,
              x_coord, y_coord, z_coord, a_coord, mat, matM_cons, matM,
              index_test_ansatz, N_neum_to_diri_psd, neum_to_diri_psd,
              neum_to_diri_psd_x, neum_to_diri_psd_y,
              neum_to_diri_psd_z, neum_to_diri_psd_a);

            OutPut("===== End Crank-Nicolson FEM-FCT Method ======="<< endl);
          }
         midpointflow(x_coord, y_coord,z_coord, a_coord,  Nx, Ny,
		       Nz,  Na, sol_psd );
        
         compute_mean_value_outflow(mean_value_outflow_indices, mean_value_outflow, sol_psd, x_coord, 
                       Nx,  Ny, Nz,  Na, &only_first_time); 
         }
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
    }                                             // adaptive time step control

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
          os.seekp(std::ios::beg);
          os << VtkBaseName << "psd.0." << m << ".vtk" << ends;
          write_vtk_file_yzlayer(Nx, Ny, Nz, Na, 0, x_coord, y_coord, z_coord, a_coord, sol_psd,os.str().c_str());
          os.seekp(std::ios::beg);
          os << VtkBaseName << "psd." << m << ".vtk" << ends;
          write_vtk_file_yzlayer(Nx, Ny, Nz, Na, cut_coord, x_coord, y_coord, z_coord, a_coord, sol_psd,os.str().c_str());
          os.seekp(std::ios::beg);
          os << VtkBaseName << m << ".txt" << ends;
          write_data_file_meanvalue(mean_value_outflow,  a_layers_coord, Nx, Ny, Nz, Na, os.str().c_str()); 
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
          //save_sol[1] = sol_c_A;
          //save_sol[2] = sol_c_B;
          //save_sol[3] = sol_c_C;
          save_sol[4] = sol_psd;
          save_N_Unknowns[0] = N_Unknowns;
          //save_N_Unknowns[1] = N_Unknowns_c_A;
          // save_N_Unknowns[2] = N_Unknowns_c_B;
          //save_N_Unknowns[3] = N_Unknowns_c_C;
          save_N_Unknowns[4] = Nodes;
          SaveData(SaveDataFileName,5,save_sol,save_N_Unknowns);
        }
      }
    }
    comp_vort =0;
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
      os.seekp(std::ios::beg);
      os << VtkBaseName << "psd." << m << ".vtk" << ends;
     write_vtk_file_yzlayer(Nx, Ny, Nz, Na, cut_coord, x_coord, y_coord, z_coord, a_coord, sol_psd,os.str().c_str());
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
