// =======================================================================
//
// Purpose:     main program
//
// Author:      Volker John   August 2000
//
// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <DirectSolver.h>
#include <DiscreteForm3D.h>
#include <FEDatabase3D.h>
#include <FESpace3D.h>
#include <FEVectFunct3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <AuxParam3D.h>
#include <Solver.h>
#include <Assemble3D.h>
#include <Output3D.h>
#include <LinAlg.h>
#include <NSE3D_ParamRout.h>
#include <BoundFace.h>
//#include <TNSE3D_ParamRout.h>
#include <RFB.h>
#include <VMS.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
// #include <malloc.h>
#include <math.h>

double bound = 0;

// ======================================================================
// utilities for main program
// ======================================================================
#include <MainUtilities.h>
#include <Upwind3D.h>
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
#include <MultiGrid3D.h>
#include <MGLevel3D.h>
#include <MultiGridScaIte.h>

#define AMG 0
#define GMG 1
#define DIRECT 2

// =======================================================================
// include current example
// =======================================================================
//#include "Examples/NSE_3D/DrivenCavity3D.h"
#include "../Examples/NSE_3D/AnsatzLinConst.h"
//#include "Examples/NSE_3D/AnsatzQuadLin.h"
//#include "Examples/NSE_3D/BSExample_alpha.h"
//#include "Examples/NSE_3D/BSExample.h"
//#include "Examples/NSE3D/StatBSExample.h"
//#include "Examples/NSE_3D/Bench3DQuaderStat.h"
//#include "Examples/NSE_3D/Bench3DQuaderStatNeum.h"
//#include "Examples/NSE_3D/Bench3DCylinderStatNeum.h"
//#include "Examples/NSE_3D/Bench3DQuaderStatNeum.h"
//#include "Examples/NSE_3D/ChannelQuadLong.h"
//#include "Examples/NSE_3D/LES_Bench3DQuaderStatNeum.h"
//#include "../Examples/TNSE_3D/LES_Bench3DQuaderNeum.h"
//#include "../Examples/TNSE_3D/ChannelStepSlip3D.h"
//#include "../Examples/TNSE_3D/ChannelStepSlipDiri3D.h"
//#include "Examples/NSE_3D/BrennstoffzelleParabol.h"
//#include "Examples/NSE_3D/BrennstoffzelleParabolGravi.h"
//#include "Examples/NSE_3D/BrennStutzenNSE.h"
//#include "Examples/NSE_3D/NoSlipExample.h"
//#include "../Examples/TNSE_3D/ChannelSlipDiri3D.h"
//#include "../Examples/TNSE_3D/ChannelStepSlipTopFreeSlipLatDiri.h"
//#include "Examples/NSE_3D/NoFlow.h"
//#include "Examples/NSE_3D/CircularChannel.h"
//#include "Examples/TNSE_3D/DrivenCavity3D_SFB_Periodic.h"
//#include "Examples/TNSE_3D/ChannelStepDiri3D.h"
//#include "../Examples/TNSE_3D/WallMountedCube.h"
//#include "../Examples/TNSE_3D/UnitCubeSlip3D.h"
//#include "../Examples/TNSE_3D/windchannel.00_bilin_inter.h"
//#include "../Examples/TNSE_3D/Calotte.h"
// #include "../Examples/TNSE_3D/Urea.h"

int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase3D *FEDatabase = new TFEDatabase3D();
  TCollection *coll;
  TBaseCell *cell;
  TFESpace3D *velocity_space, *pressure_space, *streamfunction_space;
  TFESpace3D *velocity_space_low, *pressure_space_low;
  TFESpace3D *old_u_space, *old_p_space, *projection_space;
  TFESpace3D *pressure_separation_space;
  TFESpace3D **USpaces, **PSpaces, **PsiSpaces, **ProjectionSpaces;
  TOutput3D *Output;

  double *rhs, *sol, *oldsol, tol, tolerance, *psi, *defect, *fesol, *soldiff;
  double *rhs_low, *sol_low, *old_sol, *itmethod_sol, *itmethod_rhs;
  double *rhs_high, *nosep_p;
  int i,j,k,l, N_, Len, low,ii;
  int N_Rows, N_Columns, N_U, N_P, N_Unknowns, N_V;
  int N_Active, N_NonActive, N_U_low,  N_P_low, N_Unknowns_low;
  double *l2u1, *l2u2, *l2u3, *h1u1, *h1u2, *h1u3;
  double *l2p, *h1p, *sd, *l_inf, vect[3], exactvect[3];
  int which;
  double DiffL2, DiffH1, hmin, hmax;
  char *PRM, *GEO, *MAP;
  int LEVELS, BASELEVEL;
  int ret, pde;
  double negPower;
  double x,y,max,min,sum;
  double RE_NR;
  double tau1, tau2;
  double errors[4], p1, p2, errors_mg[4],velo_l2, *save_sol[1];
  double t1, t2, res, res2, oldres, solver_time,residual, assemble_time;
  double impuls_residual,limit, total_time1, total_time2;
  int N_LinIter, save_N_Unknowns[1];

  std::ostringstream os;
  char *PsBaseName, *GrapeBaseName, *GmvBaseName, *ReadGrapeBaseName;
  char *VTKBaseName;
  char *SaveDataFileName,  *ReadDataFileName;

  double *val;
  TFEFunction3D *u1, *u2, *u3, *p, *fefct[3], *StreamFct;
  TFEFunction3D *u1_low, *u2_low, *u3_low, *p_low, *old_p, *soldiff_fe1;
  TFEFunction3D *soldiff_fe2, *soldiff_fe3;
  TFEFunction3D **U1Array, **U2Array, **U3Array, **AuxFEFunctArray;
  TFEFunction3D **PArray,  *AuxPArray, *separated_pressure_fe_funct;
  TFEFunction3D *separated_pressure_rhs_fe_funct;
  TFEVectFunct3D *u, **UArray, *u_low, *old_u, **AuxFEVectFunctArray;
  TFESpace3D *fesp[3], *ferhs[3];

  TAuxParam3D *aux;

  TSquareStructure3D *sqstructureA, *sqstructurePressSep;
  TStructure3D *structureB, *structureBT;
  TSquareMatrix3D *sqmatrixA, *SQMATRICES[9];
  TSquareMatrix3D *sqmatrixA11, *sqmatrixA12, *sqmatrixA13;
  TSquareMatrix3D *sqmatrixA21, *sqmatrixA22, *sqmatrixA23;
  TSquareMatrix3D *sqmatrixA31, *sqmatrixA32, *sqmatrixA33;
  TSquareMatrix3D **MatricesA;
  TSquareMatrix3D **MatricesA11, **MatricesA12, **MatricesA13;
  TSquareMatrix3D **MatricesA21, **MatricesA22, **MatricesA23;
  TSquareMatrix3D **MatricesA31, **MatricesA32, **MatricesA33;
  TMatrix3D *matrixB1, *matrixB2, *matrixB3, *MATRICES[9];
  TMatrix3D *matrixB1T, *matrixB2T, *matrixB3T;
  TMatrix3D **MatricesB1, **MatricesB2, **MatricesB3;
  TMatrix3D **MatricesB1T, **MatricesB2T, **MatricesB3T;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TMatrix **matrices = (TMatrix **)MATRICES;
  //TMatrix **matrices = MATRICES;

  TSquareStructure3D *sqstructureA_low;
  TStructure3D *structureB_low, *structureBT_low;
  TSquareMatrix3D *sqmatrixA_low;
  TSquareMatrix3D *sqmatrixA11_low, *sqmatrixA12_low, *sqmatrixA13_low;
  TSquareMatrix3D *sqmatrixA21_low, *sqmatrixA22_low, *sqmatrixA23_low;
  TSquareMatrix3D *sqmatrixA31_low, *sqmatrixA32_low, *sqmatrixA33_low;
  TSquareMatrix3D **MatricesA_low;
  TSquareMatrix3D **MatricesA11_low, **MatricesA12_low, **MatricesA13_low;
  TSquareMatrix3D **MatricesA21_low, **MatricesA22_low, **MatricesA23_low;
  TSquareMatrix3D **MatricesA31_low, **MatricesA32_low, **MatricesA33_low;
  TMatrix3D *matrixB1_low, *matrixB2_low, *matrixB3_low;
  TMatrix3D *matrixB1T_low, *matrixB2T_low, *matrixB3T_low;
  TMatrix3D **MatricesB1_low, **MatricesB2_low, **MatricesB3_low;
  TMatrix3D **MatricesB1T_low, **MatricesB2T_low, **MatricesB3T_low;
  MatVecProc *MatVect;
  DefectProc *Defect;
  TSquareMatrix3D *sqmatrixPressSep;

  TSquareStructure3D *sqstructureL;
  TStructure3D *structure_tilde_G, *structure_G;
  TSquareMatrix3D *sqmatrixL, **MatricesL;
  TMatrix3D *matrix_tilde_G11, *matrix_tilde_G22, *matrix_tilde_G33;
  TMatrix3D *matrix_G11, *matrix_G22, *matrix_G33;
  TMatrix3D **Matrices_tilde_G11, **Matrices_tilde_G22, **Matrices_tilde_G33;
  TMatrix3D **Matrices_G11, **Matrices_G22, **Matrices_G33;
  int N_L;

  int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;
  double **RhsArray;

  TNSE_MGLevel *MGLevel, *MGLevel_low;
  TNSE_MultiGrid *MG;

  double *RHSs[4];
  int *N_Uarray, *N_Parray;

  TDiscreteForm3D *DiscreteFormGalerkin;
  TDiscreteForm3D *DiscreteFormSDFEM;
  TDiscreteForm3D *DiscreteFormUpwind;
  TDiscreteForm3D *DiscreteFormSmagorinsky;
  TDiscreteForm3D *DiscreteFormVMSProjection;

  TDiscreteForm3D *DiscreteFormNLGalerkin;
  TDiscreteForm3D *DiscreteFormNLSDFEM;
  TDiscreteForm3D *DiscreteFormNLUpwind;
  TDiscreteForm3D *DiscreteFormNLSmagorinsky;
  TDiscreteForm3D *DiscreteFormNLSDFEM_DivDiv;
  TDiscreteForm3D *DiscreteFormNLVMSProjection;

  TDiscreteForm3D *DiscreteFormPressSep;
  TDiscreteForm3D *DiscreteFormAuxProbPressSep;

  TDiscreteForm3D *DiscreteFormNSRFBRhs;

  TDiscreteForm3D *DiscreteForm, *DiscreteForm0;

  BoundCondFunct3D *BoundaryConditions[3], *BoundaryConditionsPressureSeparation[1];
  BoundValueFunct3D *BoundValues[3], *BoundaryValuesPressureSeparation[1];
  CoeffFct3D *Coefficients[1];
  double average;

  TItMethod *itmethod, *prec, *Auxprec, *Auxitmethod;
  int Max_It, FirstSolve, zerostart;
  double omega, alpha[2], alpha_fine[2],cd,cl,dp;
  int N_Paramters=2,n_aux, **downwind;
  double Parameters[4],delta0,delta1, reatt_pt;
  double convergence_speed, residuals[10], firsterror,lasterror;
  double firsterrorl2,lasterrorl2,p3,p4;
  int last_digit_ite,slow_conv,last_sq,mg_level,mg_type, pre_calculation=1;
  int calculations=1,ll;
  int velocity_space_code, pressure_space_code;
  int pressure_separation = 0, N_P_sep;
  double *separated_pressure_array, *separated_pressure_aux_array;
  double *pressure_aux_array;
  double *rhsPressSep;

  TMultiGrid3D *AuxMG;
  TMGLevel3D *AuxMGLevel;
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
  double Cd, Cl, dP1[4], dP2[4];
  #endif
  #ifdef  __CHANNEL_OBSTACLE__
  double Cd, Cl, dP1[4], dP2[4];
  #endif

#ifdef __STUTZEN__
  TBaseCell **Cells;
  int N_CellsOld, N_CellsNew;
  int N_Vertices;
  double xm, ym, zm;
  double xp, yp, zp;
  double r1, r2;
#endif

#ifdef __WALL_MOUNTED_CUBE__
  TBaseCell **Cells;
  int N_CellsOld, N_CellsNew, N_Vertices;
  double xm, ym, zm, xp, yp, zp; 
#endif

  os << " ";

  //======================================================================
  // read parameter file
  //======================================================================

  if (argc >= 2)
    ret = Domain->ReadParam(argv[1]);
  else
    ret = Domain->ReadParam(Readin);

  #ifdef  __BENCH3_CYLINDER__
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 1;
  #endif

  if (ret == -1)
  {
    exit(-1);
  }

  OpenFiles();
  OutFile.setf(std::ios::scientific);

#ifdef __STUTZEN__
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = 1234;
#endif

  RE_NR = TDatabase::ParamDB->RE_NR;
  delta0 = TDatabase::ParamDB->DELTA0;
  delta1 = TDatabase::ParamDB->DELTA1;
  convergence_speed = TDatabase::ParamDB->SC_DIV_FACTOR;

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
  GmvBaseName = TDatabase::ParamDB->GMVBASENAME;
  VTKBaseName = TDatabase::ParamDB->VTKBASENAME;
  ReadGrapeBaseName = TDatabase::ParamDB->READGRAPEBASENAME;
  SaveDataFileName = TDatabase::ParamDB->SAVE_DATA_FILENAME;
  ReadDataFileName = TDatabase::ParamDB->READ_DATA_FILENAME;

  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SADDLE;
  if (mg_type)
    mg_level = 0;
  else
    mg_level = -1;
  LEVELS = TDatabase::ParamDB->LEVELS;
  BASELEVEL = TDatabase::ParamDB->UNIFORM_STEPS;
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

  RhsArray = new double* [LEVELS+1];
  N_Uarray = new int[LEVELS+1];
  N_Parray = new int[LEVELS+1];

  USpaces = new TFESpace3D*[LEVELS+1];
  PSpaces = new TFESpace3D*[LEVELS+1];
  PsiSpaces = new TFESpace3D*[LEVELS+1];

  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      MatricesA = new TSquareMatrix3D*[LEVELS+1];

      MatricesB1 = new TMatrix3D*[LEVELS+1];
      MatricesB2 = new TMatrix3D*[LEVELS+1];
      MatricesB3 = new TMatrix3D*[LEVELS+1];
      MatVect = MatVect_NSE1;
      Defect = Defect_NSE1;
      break;

    case 2:
      MatricesA = new TSquareMatrix3D*[LEVELS+1];

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

      MatricesB1 = new TMatrix3D*[LEVELS+1];
      MatricesB2 = new TMatrix3D*[LEVELS+1];
      MatricesB3 = new TMatrix3D*[LEVELS+1];
      MatricesB1T = new TMatrix3D*[LEVELS+1];
      MatricesB2T = new TMatrix3D*[LEVELS+1];
      MatricesB3T = new TMatrix3D*[LEVELS+1];
      MatVect = MatVect_NSE4;
      Defect = Defect_NSE4;
      break;
  }                              // endswitch

  // matrices for VMS_PROJECTION
  if (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION)
  {
    MatricesL = new TSquareMatrix3D*[LEVELS+1];
    Matrices_tilde_G11 = new TMatrix3D*[LEVELS+1];
    Matrices_tilde_G22 = new TMatrix3D*[LEVELS+1];
    Matrices_tilde_G33 = new TMatrix3D*[LEVELS+1];
    Matrices_G11 = new TMatrix3D*[LEVELS+1];
    Matrices_G22 = new TMatrix3D*[LEVELS+1];
    Matrices_G33 = new TMatrix3D*[LEVELS+1];
    ProjectionSpaces = new TFESpace3D*[LEVELS+1];
  }

  downwind = new int*[LEVELS+1];

  //======================================================================
  // read boundary parameterization and initialize coarse grid
  //======================================================================
  Domain->Init(PRM, GEO);

  //======================================================================
  // initialize all discrete forms
  //======================================================================

  InitializeDiscreteForms(
    DiscreteFormGalerkin, DiscreteFormSDFEM,
    DiscreteFormUpwind, DiscreteFormSmagorinsky,
    DiscreteFormVMSProjection,
    DiscreteFormNLGalerkin, DiscreteFormNLSDFEM,
    DiscreteFormNLUpwind, DiscreteFormNLSmagorinsky, 
    DiscreteFormNLVMSProjection,
    DiscreteFormNLSDFEM_DivDiv,
    DiscreteFormPressSep,
    DiscreteFormAuxProbPressSep,
    DiscreteFormNSRFBRhs,
    LinCoeffs, TDatabase::ParamDB->NSTYPE);

  BoundaryConditions[0] = BoundCondition;
  BoundaryConditions[1] = BoundCondition;
  BoundaryConditions[2] = BoundCondition;

  BoundValues[0] = U1BoundValue;
  BoundValues[1] = U2BoundValue;
  BoundValues[2] = U3BoundValue;

  BoundaryConditionsPressureSeparation[0] = BoundaryConditionPressSep3D;
  BoundaryValuesPressureSeparation[0] = BoundaryValuePressSep3D;

  Coefficients[0] = LinCoeffs;
  // refine up to user defined coarsest level

  for(i=0;i<TDatabase::ParamDB->SC_COARSEST_LEVEL_SADDLE;i++)
  {
    Domain->RegRefineAll();
    #ifdef __DC3D_SFB_PERIODIC__
    coll=Domain->GetCollection(It_EQ, i);
    SetPeriodicFaceJoints(coll);
    #endif
  }
  // initialize solver parameters

  limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
  Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
  omega = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SADDLE;
  velocity_space_code =   TDatabase::ParamDB->VELOCITY_SPACE;

  Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
  Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;

  if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    i=1;
    MG = new TNSE_MultiGrid(i, N_Paramters, Parameters);
  }
  Parameters[2] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
  Parameters[3] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
  if ((TDatabase::ParamDB->PRESSURE_SEPARATION==3)||
    (TDatabase::ParamDB->PRESSURE_SEPARATION==4)||
    (TDatabase::ParamDB->PRESSURE_SEPARATION==5))
  {
    AuxMG = new TMultiGrid3D(i, N_Paramters, Parameters+2);
  }

  FirstSolve = TDatabase::ParamDB->SC_FIRST_SOLUTION_LEVEL_SADDLE;

  #ifdef  __WINDCHANNEL__
    ReadExperimentalBoundaryConditions();
  #endif

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
    assemble_time = 0.0;
    N_LinIter = 0;
    OutPut("memory before: " << setw(10) << GetMemory() << endl);
    for (j=0;j<10;j++)
      residuals[j]=1e10;
    slow_conv = 0;

    // refine grid if level is greater than 0
    if (i)
      Domain->RegRefineAll();

    coll=Domain->GetCollection(It_Finest, 0);
#ifdef __STUTZEN__
      N_CellsOld = coll->GetN_Cells();
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

        xm -= 17;
        ym -= 120;
        r1 = xm*xm + ym*ym;

        xm += 17; xm -= 48;
        ym += 120; ym -= 20;
        r2 = xm*xm + ym*ym;

        if( (zm > 15.32) || (r1<16) || (r2<16) )
        {
          N_CellsNew++;
        }
      } // endfor j
      OutPut("N_CellsNew: " << N_CellsNew << endl);

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

        xm -= 17;
        ym -= 120;
        r1 = xm*xm + ym*ym;

        xm += 17; xm -= 48;
        ym += 120; ym -= 20;
        r2 = xm*xm + ym*ym;

        if( (zm > 15.32) || (r1<16) || (r2<16) )
        {
          Cells[N_CellsNew] = cell;
          N_CellsNew++;
        }
      } // endfor j

      delete coll;
      coll = new TCollection(N_CellsNew, Cells);
#endif
#ifdef __WALL_MOUNTED_CUBE__
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
  } // endfor j
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
  } // endfor j
  
  delete coll;
  coll = new TCollection(N_CellsNew, Cells);
#endif
    OutPut("cells " << coll->GetN_Cells()<< endl);

    Output = new TOutput3D(1, 1, 1, 1,Domain);
    cout << endl << endl;

    #ifdef __DC3D_SFB_PERIODIC__
    SetPeriodicFaceJoints(coll);
    #endif

    if(TDatabase::ParamDB->WRITE_PS)
    {
      // write grid into an Postscript file
      os.seekp(std::ios::beg);
      os << PsBaseName << i << ".ps" << ends;
      Domain->PS(os.str().c_str(),It_Finest,0);
    }

    // get spaces for low order disc on finest geo grid
    if (mg_type==1)
    {
      velocity_space_low = new TFESpace3D(coll,Name,UString,BoundCondition,
        Non_USpace,1);
      pressure_space_low = new TFESpace3D(coll,Name,PString,BoundCondition,
        DiscP_PSpace,0);
    }
    // get spaces of high order disc on finest geo grid
    if ((i>=FirstSolve)||(mg_type==0))
      GetVelocityAndPressureSpace3D(coll,BoundCondition,
        velocity_space,
        pressure_space, &pressure_space_code,
        TDatabase::ParamDB->VELOCITY_SPACE,
        TDatabase::ParamDB->PRESSURE_SPACE);

    TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
    // fe space for stream function
    streamfunction_space = new TFESpace3D(coll,Name,PsiString,BoundCondition,
      1);
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
    {
      if (i<FirstSolve)
        pressure_space_code=0;
      OutPut("pressure_space_code " << pressure_space_code << endl);
      // allocate finite element space for separated pressure
      switch (pressure_space_code)
      {
        case 0:
          pressure_separation_space = new TFESpace3D(coll,Name,PsiString,BoundaryConditionPressSep3D,1);
          break;
        case 1:
        case -11:
          pressure_separation_space = new TFESpace3D(coll,Name,PsiString,BoundaryConditionPressSep3D,
            2);
          break;
        default:
          OutPut("case for pressure_space_code not implemented" << endl);
          break;
      }
      N_P_sep = pressure_separation_space->GetN_DegreesOfFreedom();
      separated_pressure_array = new double[N_P_sep];
      memset(separated_pressure_array,0, N_P_sep*SizeOfDouble);
      separated_pressure_aux_array = new double[N_P_sep];
      memset(separated_pressure_aux_array,0, N_P_sep*SizeOfDouble);
      separated_pressure_fe_funct = new TFEFunction3D(pressure_separation_space,
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
        sqstructurePressSep = new TSquareStructure3D(pressure_separation_space);
        sqmatrixPressSep = new TSquareMatrix3D(sqstructurePressSep);
        // rhs for auxiliary problem
        rhsPressSep = new double[N_P_sep];
        memset(rhsPressSep, 0, N_P_sep*SizeOfDouble);
        separated_pressure_rhs_fe_funct = new TFEFunction3D(pressure_separation_space,
          PsepString, PsepString,rhsPressSep,
          N_P_sep);
        pressure_aux_array = new double[N_P];
        memset(pressure_aux_array,0, N_P*SizeOfDouble);
        n_aux = 4;
        AuxMGLevel = new TMGLevel3D(i, sqmatrixPressSep,
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
          aux = new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
          RHSs[0] = rhsPressSep;

          // assemble
          Assemble3D(N_FESpaces, fesp,
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

    // build matrices for high order disc
    if ((i>=FirstSolve)||(mg_type==0))
    {
      // matrix structures
      structureB = new TStructure3D(pressure_space, velocity_space);
      structureBT = new TStructure3D(velocity_space, pressure_space);
      sqstructureA = new TSquareStructure3D(velocity_space);
      sqstructureA->Sort();

      // allocate matrices
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          matrixB1 = new TMatrix3D(structureB);
          matrixB2 = new TMatrix3D(structureB);
          matrixB3 = new TMatrix3D(structureB);

          MatricesB1[mg_level] = matrixB1;
          MatricesB2[mg_level] = matrixB2;
          MatricesB3[mg_level] = matrixB3;

          sqmatrixA = new TSquareMatrix3D(sqstructureA);

          MatricesA[mg_level] = sqmatrixA;
          break;

        case 2:
          matrixB1 = new TMatrix3D(structureB);
          matrixB2 = new TMatrix3D(structureB);
          matrixB3 = new TMatrix3D(structureB);
          matrixB1T = new TMatrix3D(structureBT);
          matrixB2T = new TMatrix3D(structureBT);
          matrixB3T = new TMatrix3D(structureBT);

          MatricesB1[mg_level] = matrixB1;
          MatricesB2[mg_level] = matrixB2;
          MatricesB3[mg_level] = matrixB3;
          MatricesB1T[mg_level] = matrixB1T;
          MatricesB2T[mg_level] = matrixB2T;
          MatricesB3T[mg_level] = matrixB3T;

          sqmatrixA = new TSquareMatrix3D(sqstructureA);

          MatricesA[mg_level] = sqmatrixA;
          break;

        case 3:
          matrixB1 = new TMatrix3D(structureB);
          matrixB2 = new TMatrix3D(structureB);
          matrixB3 = new TMatrix3D(structureB);

          MatricesB1[mg_level] = matrixB1;
          MatricesB2[mg_level] = matrixB2;
          MatricesB3[mg_level] = matrixB3;

          sqmatrixA11 = new TSquareMatrix3D(sqstructureA);
          sqmatrixA12 = new TSquareMatrix3D(sqstructureA);
          sqmatrixA13 = new TSquareMatrix3D(sqstructureA);
          sqmatrixA21 = new TSquareMatrix3D(sqstructureA);
          sqmatrixA22 = new TSquareMatrix3D(sqstructureA);
          sqmatrixA23 = new TSquareMatrix3D(sqstructureA);
          sqmatrixA31 = new TSquareMatrix3D(sqstructureA);
          sqmatrixA32 = new TSquareMatrix3D(sqstructureA);
          sqmatrixA33 = new TSquareMatrix3D(sqstructureA);

          MatricesA11[mg_level] = sqmatrixA11;
          MatricesA12[mg_level] = sqmatrixA12;
          MatricesA13[mg_level] = sqmatrixA13;
          MatricesA21[mg_level] = sqmatrixA21;
          MatricesA22[mg_level] = sqmatrixA22;
          MatricesA23[mg_level] = sqmatrixA23;
          MatricesA31[mg_level] = sqmatrixA31;
          MatricesA32[mg_level] = sqmatrixA32;
          MatricesA33[mg_level] = sqmatrixA33;
          break;

        case 4:
          matrixB1 = new TMatrix3D(structureB);
          matrixB2 = new TMatrix3D(structureB);
          matrixB3 = new TMatrix3D(structureB);
          matrixB1T = new TMatrix3D(structureBT);
          matrixB2T = new TMatrix3D(structureBT);
          matrixB3T = new TMatrix3D(structureBT);

          MatricesB1[mg_level] = matrixB1;
          MatricesB2[mg_level] = matrixB2;
          MatricesB3[mg_level] = matrixB3;
          MatricesB1T[mg_level] = matrixB1T;
          MatricesB2T[mg_level] = matrixB2T;
          MatricesB3T[mg_level] = matrixB3T;

          sqmatrixA11 = new TSquareMatrix3D(sqstructureA);
          sqmatrixA12 = new TSquareMatrix3D(sqstructureA);
          sqmatrixA13 = new TSquareMatrix3D(sqstructureA);
          sqmatrixA21 = new TSquareMatrix3D(sqstructureA);
          sqmatrixA22 = new TSquareMatrix3D(sqstructureA);
          sqmatrixA23 = new TSquareMatrix3D(sqstructureA);
          sqmatrixA31 = new TSquareMatrix3D(sqstructureA);
          sqmatrixA32 = new TSquareMatrix3D(sqstructureA);
          sqmatrixA33 = new TSquareMatrix3D(sqstructureA);

          MatricesA11[mg_level] = sqmatrixA11;
          MatricesA12[mg_level] = sqmatrixA12;
          MatricesA13[mg_level] = sqmatrixA13;
          MatricesA21[mg_level] = sqmatrixA21;
          MatricesA22[mg_level] = sqmatrixA22;
          MatricesA23[mg_level] = sqmatrixA23;
          MatricesA31[mg_level] = sqmatrixA31;
          MatricesA32[mg_level] = sqmatrixA32;
          MatricesA33[mg_level] = sqmatrixA33;
          break;
      }
    }

    // build matrices for low order disc
    if (mg_type==1)
    {
      // matrix structures
      structureB_low = new TStructure3D(pressure_space_low, velocity_space_low);
      structureBT_low = new TStructure3D(velocity_space_low, pressure_space_low);
      sqstructureA_low = new TSquareStructure3D(velocity_space_low);
      sqstructureA_low->Sort();

      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          matrixB1_low = new TMatrix3D(structureB_low);
          matrixB2_low = new TMatrix3D(structureB_low);
          matrixB3_low = new TMatrix3D(structureB_low);

          MatricesB1[i] = matrixB1_low;
          MatricesB2[i] = matrixB2_low;
          MatricesB3[i] = matrixB3_low;

          sqmatrixA_low = new TSquareMatrix3D(sqstructureA_low);

          MatricesA[i] = sqmatrixA_low;
          break;

        case 2:
          matrixB1_low = new TMatrix3D(structureB_low);
          matrixB2_low = new TMatrix3D(structureB_low);
          matrixB3_low = new TMatrix3D(structureB_low);
          matrixB1T_low = new TMatrix3D(structureBT_low);
          matrixB2T_low = new TMatrix3D(structureBT_low);
          matrixB3T_low = new TMatrix3D(structureBT_low);

          MatricesB1[i] = matrixB1_low;
          MatricesB2[i] = matrixB2_low;
          MatricesB3[i] = matrixB3_low;
          MatricesB1T[i] = matrixB1T_low;
          MatricesB2T[i] = matrixB2T_low;
          MatricesB3T[i] = matrixB3T_low;

          sqmatrixA_low = new TSquareMatrix3D(sqstructureA_low);

          MatricesA[i] = sqmatrixA_low;
          break;

        case 3:
          matrixB1_low = new TMatrix3D(structureB_low);
          matrixB2_low = new TMatrix3D(structureB_low);
          matrixB3_low = new TMatrix3D(structureB_low);

          MatricesB1[i] = matrixB1_low;
          MatricesB2[i] = matrixB2_low;
          MatricesB3[i] = matrixB3_low;

          sqmatrixA11_low = new TSquareMatrix3D(sqstructureA_low);
          sqmatrixA12_low = new TSquareMatrix3D(sqstructureA_low);
          sqmatrixA13_low = new TSquareMatrix3D(sqstructureA_low);
          sqmatrixA21_low = new TSquareMatrix3D(sqstructureA_low);
          sqmatrixA22_low = new TSquareMatrix3D(sqstructureA_low);
          sqmatrixA23_low = new TSquareMatrix3D(sqstructureA_low);
          sqmatrixA31_low = new TSquareMatrix3D(sqstructureA_low);
          sqmatrixA32_low = new TSquareMatrix3D(sqstructureA_low);
          sqmatrixA33_low = new TSquareMatrix3D(sqstructureA_low);

          MatricesA11[i] = sqmatrixA11_low;
          MatricesA12[i] = sqmatrixA12_low;
          MatricesA13[i] = sqmatrixA13_low;
          MatricesA21[i] = sqmatrixA21_low;
          MatricesA22[i] = sqmatrixA22_low;
          MatricesA23[i] = sqmatrixA23_low;
          MatricesA31[i] = sqmatrixA31_low;
          MatricesA32[i] = sqmatrixA32_low;
          MatricesA33[i] = sqmatrixA33_low;
          break;

        case 4:
          matrixB1_low = new TMatrix3D(structureB_low);
          matrixB2_low = new TMatrix3D(structureB_low);
          matrixB3_low = new TMatrix3D(structureB_low);
          matrixB1T_low = new TMatrix3D(structureBT_low);
          matrixB2T_low = new TMatrix3D(structureBT_low);
          matrixB3T_low = new TMatrix3D(structureBT_low);

          MatricesB1[i] = matrixB1_low;
          MatricesB2[i] = matrixB2_low;
          MatricesB3[i] = matrixB3_low;
          MatricesB1T[i] = matrixB1T_low;
          MatricesB2T[i] = matrixB2T_low;
          MatricesB3T[i] = matrixB3T_low;

          sqmatrixA11_low = new TSquareMatrix3D(sqstructureA_low);
          sqmatrixA12_low = new TSquareMatrix3D(sqstructureA_low);
          sqmatrixA13_low = new TSquareMatrix3D(sqstructureA_low);
          sqmatrixA21_low = new TSquareMatrix3D(sqstructureA_low);
          sqmatrixA22_low = new TSquareMatrix3D(sqstructureA_low);
          sqmatrixA23_low = new TSquareMatrix3D(sqstructureA_low);
          sqmatrixA31_low = new TSquareMatrix3D(sqstructureA_low);
          sqmatrixA32_low = new TSquareMatrix3D(sqstructureA_low);
          sqmatrixA33_low = new TSquareMatrix3D(sqstructureA_low);

          MatricesA11[i] = sqmatrixA11_low;
          MatricesA12[i] = sqmatrixA12_low;
          MatricesA13[i] = sqmatrixA13_low;
          MatricesA21[i] = sqmatrixA21_low;
          MatricesA22[i] = sqmatrixA22_low;
          MatricesA23[i] = sqmatrixA23_low;
          MatricesA31[i] = sqmatrixA31_low;
          MatricesA32[i] = sqmatrixA32_low;
          MatricesA33[i] = sqmatrixA33_low;
          break;
      }
    }                            // end if (mg_type==1)

    if ((i>=FirstSolve)||(mg_type==0))
    {
      N_Unknowns = 3*N_U + N_P;
      OutPut("dof velocity : "<< setw(10) << 3*N_U << endl);
      OutPut("dof pressure : "<< setw(10) << N_P << endl);
      OutPut("dof all      : "<< setw(10) << N_Unknowns  << endl);
    }

    // matrices for VMS_PROJECTION
    if ((i>=FirstSolve)||(mg_type==0))
    {
    if (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION)
    {
	switch(TDatabase::ParamDB->VMS_LARGE_VELOCITY_SPACE)
	{
	    case -1:
		projection_space = new TFESpace3D(coll, Name, UString, BoundCondition,
					      DiscP_PSpace,102);
		break;
	    case 0:
		projection_space = new TFESpace3D(coll, Name, UString, BoundCondition,
					      DiscP_PSpace,0);
		break;
	    case 1:
		projection_space = new TFESpace3D(coll, Name, UString, BoundCondition,
						  DiscP_PSpace,1);
		break;
	    default:
		OutPut("VMS_LARGE_VELOCITY_SPACE: " <<
		       TDatabase::ParamDB->VMS_LARGE_VELOCITY_SPACE <<
		       " not available !!!" << endl);
		exit(4711);
	}

      ProjectionSpaces[mg_level] = projection_space;
      sqstructureL = new TSquareStructure3D(projection_space);
      sqstructureL->Sort();
      structure_tilde_G = new TStructure3D(velocity_space, projection_space);
      structure_G = new TStructure3D(projection_space, velocity_space);
      sqmatrixL = new TSquareMatrix3D(sqstructureL);
      MatricesL[mg_level] = sqmatrixL;
      matrix_tilde_G11 = new TMatrix3D(structure_tilde_G);
      Matrices_tilde_G11[mg_level] = matrix_tilde_G11;
      matrix_tilde_G22 = new TMatrix3D(structure_tilde_G);
      Matrices_tilde_G22[mg_level] = matrix_tilde_G22;
      matrix_tilde_G33 = new TMatrix3D(structure_tilde_G);
      Matrices_tilde_G33[mg_level] = matrix_tilde_G33;
      matrix_G11 = new TMatrix3D(structure_G);
      Matrices_G11[mg_level] = matrix_G11;
      matrix_G22 = new TMatrix3D(structure_G);
      Matrices_G22[mg_level] = matrix_G22;
      matrix_G33 = new TMatrix3D(structure_G);
      Matrices_G33[mg_level] = matrix_G33;
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
      N_Unknowns_low = 3*N_U_low + N_P_low;
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
          n_aux= 4;
        else
          n_aux= 2;

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
                matrixB3_low,
                structureBT_low,
                rhs_low,
                sol_low,
                n_aux, alpha, -1, 0, NULL,downwind[i]);
              if (i==0)
                MG->AddLevel(MGLevel_low);
              else
                MG->ReplaceLevel(i,MGLevel_low);
            }
            // high order disc
            if ((i>=FirstSolve)||(mg_type==0))
            {
              MGLevel = new TNSE_MGLevel1(mg_level, sqmatrixA,
                matrixB1, matrixB2, matrixB3,
                structureBT,
                rhs_high, sol,
                n_aux, alpha_fine,
                velocity_space_code,
                pressure_space_code,
                NULL,downwind[i]);
              MG->AddLevel(MGLevel);
            }
            break;

          case 2:
            if (mg_type==1)
            {
              MGLevel_low = new TNSE_MGLevel2(i, sqmatrixA_low,
                matrixB1_low, matrixB2_low,
                matrixB3_low,
                matrixB1T_low, matrixB2T_low,
                matrixB3T_low,
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
                matrixB1, matrixB2, matrixB3,
                matrixB1T, matrixB2T,  matrixB3T,
                rhs_high, sol,
                n_aux, alpha_fine,
                velocity_space_code,
                pressure_space_code, NULL,downwind[i]);
              MG->AddLevel(MGLevel);
            }
            break;

          case 3:
            if (mg_type==1)
            {
              MGLevel_low = new TNSE_MGLevel3(i,
                sqmatrixA11_low, sqmatrixA12_low,
                sqmatrixA13_low,
                sqmatrixA21_low, sqmatrixA22_low,
                sqmatrixA23_low,
                sqmatrixA31_low, sqmatrixA32_low,
                sqmatrixA33_low,
                matrixB1_low, matrixB2_low,
                matrixB3_low,
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
                sqmatrixA13,
                sqmatrixA21, sqmatrixA22,
                sqmatrixA23,
                sqmatrixA31, sqmatrixA32,
                sqmatrixA33,
                matrixB1, matrixB2,
                matrixB3,
                structureBT,
                rhs_high,  sol,
                n_aux, alpha_fine,
                velocity_space_code,
                pressure_space_code, NULL,downwind[i]);
              MG->AddLevel(MGLevel);
            }
            break;

          case 4:
            // low order disc
            if (mg_type==1)
            {
              MGLevel_low = new TNSE_MGLevel4(i,
                sqmatrixA11_low, sqmatrixA12_low,
                sqmatrixA13_low,
                sqmatrixA21_low, sqmatrixA22_low,
                sqmatrixA23_low,
                sqmatrixA31_low, sqmatrixA32_low,
                sqmatrixA33_low,
                matrixB1_low, matrixB2_low,
                matrixB3_low,
                matrixB1T_low, matrixB2T_low,
                matrixB3T_low,
                rhs_low, sol_low,
                n_aux, alpha, -1, 0, coll,downwind[i]);
              if (i==0)
                MG->AddLevel(MGLevel_low);
              else
                MG->ReplaceLevel(i,MGLevel_low);
            }
            // high order disc
            if ((i>=FirstSolve)||(mg_type==0))
            {
              MGLevel = new TNSE_MGLevel4(mg_level,
                sqmatrixA11, sqmatrixA12,
                sqmatrixA13,
                sqmatrixA21, sqmatrixA22,
                sqmatrixA23,
                sqmatrixA31, sqmatrixA32,
                sqmatrixA33,
                matrixB1, matrixB2, matrixB3,
                matrixB1T, matrixB2T, matrixB3T,
                rhs_high, sol,
                n_aux, alpha_fine,
                velocity_space_code,
                pressure_space_code, coll,downwind[i]);
              MG->AddLevel(MGLevel);
            }
            break;
        }                        // end switch(NSTYPE)
        break;
    }

    // build new fe functions
    // high order fe space
    if ((i>=FirstSolve)||(mg_type==0))
    {
      u = new TFEVectFunct3D(velocity_space, UString, UString, sol, N_U, 3);
      u1 = u->GetComponent(0);
      u2 = u->GetComponent(1);
      u3 = u->GetComponent(2);
      p = new TFEFunction3D(pressure_space, PString, PString, sol+3*N_U, N_P);

      U1Array[mg_level] = u1;
      U2Array[mg_level] = u2;
      U3Array[mg_level] = u3;
      PArray[mg_level] = p;
      UArray[mg_level] = u;
    }

    // low order fe space
    if (mg_type==1)
    {
      u_low = new TFEVectFunct3D(velocity_space_low, UString, UString, sol_low, N_U_low, 3);
      u1_low = u_low->GetComponent(0);
      u2_low = u_low->GetComponent(1);
      u3_low = u_low->GetComponent(2);
      p_low = new TFEFunction3D(pressure_space_low, PString, PString, sol_low+3*N_U_low, N_P_low);
      U1Array[i] = u1_low;
      U2Array[i] = u2_low;
      U3Array[i] = u3_low;
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
      Prolongate(old_u_space, USpaces[mg_level],
        old_u->GetComponent(2)->GetValues(), U3Array[mg_level]->GetValues(),
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
      for (k=0;k<3*N_U;k++)
        sol[k] = 0.0;
      for (k=0;k<N_P;k++)
        sol[3*N_U+k] = 0.0;
      fesol = new double[N_Unknowns];
      memset(fesol, 0, N_Unknowns*SizeOfDouble);
      soldiff = new double[N_Unknowns];
      memset(soldiff, 0, N_Unknowns*SizeOfDouble);
      soldiff_fe1 = new TFEFunction3D(velocity_space, DString, DString, soldiff,N_U);
      soldiff_fe2 = new TFEFunction3D(velocity_space, DString, DString, soldiff+N_U,N_U);
      soldiff_fe3 = new TFEFunction3D(velocity_space, DString, DString, soldiff+2*N_U,N_U);
      pre_calculation = 1;
      calculations = 2;
      TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SADDLE = 5e-13;
    }

    // memset(sol, 0, N_Unknowns*SizeOfDouble);
    //for (k=0;k<3*N_U;k++)
    // sol[k] = 1.0;

    // restrict solution to all grids
    if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
      MG->RestrictToAllGrids();

    // if no solution on this grid, continue
    if(FirstSolve>i)
      continue;

    if(TDatabase::ParamDB->READ_GRAPE_FILE)
    {
      AuxFEFunctArray = new TFEFunction3D*[1];
      AuxFEFunctArray[0] = PArray[mg_level];
      AuxFEVectFunctArray = new TFEVectFunct3D*[1];
      AuxFEVectFunctArray[0] = UArray[mg_level];
      ReadGrapeFile3D(ReadGrapeBaseName, 1 , 1 , AuxFEFunctArray,AuxFEVectFunctArray);
      TDatabase::ParamDB->READ_GRAPE_FILE = 0;
      OutPut("read u " << Ddot(3*N_U,sol,sol)<< endl);
      memcpy(oldsol,sol, N_Unknowns*SizeOfDouble);
      if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
        MG->RestrictToAllGrids();
    }

    if (TDatabase::ParamDB->READ_DATA)
    {
	save_sol[0] = sol;
	save_N_Unknowns[0] = N_Unknowns;
	ReadData(ReadDataFileName,1,save_sol,save_N_Unknowns);
	memcpy(oldsol,sol,N_Unknowns*SizeOfDouble);
    }
//     u1->Interpolate(InitialU1);
    memcpy(oldsol,sol,N_Unknowns*SizeOfDouble);
 
    // build the discretizations
    t1 = GetTime();
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
      RHSs[3] = rhs + 3*N_U;
      memset(rhs, 0, (3*N_U+N_P)*SizeOfDouble);

      fesp[0] = USpaces[k];
      fesp[1] = PSpaces[k];

      fefct[0] = U1Array[k];
      fefct[1] = U2Array[k];
      fefct[2] = U3Array[k];
      ferhs[0] = USpaces[k];
      ferhs[1] = USpaces[k];
      ferhs[2] = USpaces[k];

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
          case SDFEM_DIVDIV:
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
          MATRICES[2] = MatricesB3[k];

          SQMATRICES[0]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();

          N_SquareMatrices = 1;
          N_RectMatrices = 3;

          N_Rhs = 3;
          N_FESpaces = 2;
	  if (DiscreteForm == DiscreteFormVMSProjection)
	  {
	      SQMATRICES[1] = MatricesL[k];
	      MATRICES[3] = Matrices_tilde_G11[k];
	      MATRICES[4] = Matrices_tilde_G22[k];
	      MATRICES[5] = Matrices_tilde_G33[k];
	      MATRICES[6] = Matrices_G11[k];
	      MATRICES[7] = Matrices_G22[k];
	      MATRICES[8] = Matrices_G33[k];
	      
	      SQMATRICES[1]->Reset();
	      MATRICES[3]->Reset();
	      MATRICES[4]->Reset();
	      MATRICES[5]->Reset();
	      MATRICES[6]->Reset();
	      MATRICES[7]->Reset();
	      MATRICES[8]->Reset();
	      
	      N_SquareMatrices = 2;
	      N_RectMatrices = 9;
	      
	      N_FESpaces = 3;
	      fesp[2] = ProjectionSpaces[k];
	  }
          break;

        case 2:
          SQMATRICES[0] = MatricesA[k];
          MATRICES[0] = MatricesB1[k];
          MATRICES[1] = MatricesB2[k];
          MATRICES[2] = MatricesB3[k];
          MATRICES[3] = MatricesB1T[k];
          MATRICES[4] = MatricesB2T[k];
          MATRICES[5] = MatricesB3T[k];

          SQMATRICES[0]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();
          MATRICES[4]->Reset();
          MATRICES[5]->Reset();

          N_SquareMatrices = 1;
          N_RectMatrices = 6;

          N_Rhs = 3;
          N_FESpaces = 2;
          break;

        case 3:
          SQMATRICES[0] = MatricesA11[k];
          SQMATRICES[1] = MatricesA12[k];
          SQMATRICES[2] = MatricesA13[k];
          SQMATRICES[3] = MatricesA21[k];
          SQMATRICES[4] = MatricesA22[k];
          SQMATRICES[5] = MatricesA23[k];
          SQMATRICES[6] = MatricesA31[k];
          SQMATRICES[7] = MatricesA32[k];
          SQMATRICES[8] = MatricesA33[k];
          MATRICES[0] = MatricesB1[k];
          MATRICES[1] = MatricesB2[k];
          MATRICES[2] = MatricesB3[k];

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          SQMATRICES[4]->Reset();
          SQMATRICES[5]->Reset();
          SQMATRICES[6]->Reset();
          SQMATRICES[7]->Reset();
          SQMATRICES[8]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();

          N_SquareMatrices = 9;
          N_RectMatrices = 3;

          N_Rhs = 3;
          N_FESpaces = 2;
          break;

        case 4:
          SQMATRICES[0] = MatricesA11[k];
          SQMATRICES[1] = MatricesA12[k];
          SQMATRICES[2] = MatricesA13[k];
          SQMATRICES[3] = MatricesA21[k];
          SQMATRICES[4] = MatricesA22[k];
          SQMATRICES[5] = MatricesA23[k];
          SQMATRICES[6] = MatricesA31[k];
          SQMATRICES[7] = MatricesA32[k];
          SQMATRICES[8] = MatricesA33[k];
          MATRICES[0] = MatricesB1[k];
          MATRICES[1] = MatricesB2[k];
          MATRICES[2] = MatricesB3[k];
          MATRICES[3] = MatricesB1T[k];
          MATRICES[4] = MatricesB2T[k];
          MATRICES[5] = MatricesB3T[k];

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          SQMATRICES[4]->Reset();
          SQMATRICES[5]->Reset();
          SQMATRICES[6]->Reset();
          SQMATRICES[7]->Reset();
          SQMATRICES[8]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();
          MATRICES[4]->Reset();
          MATRICES[5]->Reset();

          N_SquareMatrices = 9;
          N_RectMatrices = 6;

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
        aux =  new TAuxParam3D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo,
          NSN_ParamFctVelo_GradVelo,
          NSN_FEValuesVelo_GradVelo,
          fesp, fefct,
          NSFctVelo_GradVelo,
          NSFEFctIndexVelo_GradVelo, NSFEMultiIndexVelo_GradVelo,
          NSN_ParamsVelo_GradVelo, NSBeginParamVelo_GradVelo);
      }
      else
      {
        aux =  new TAuxParam3D(NSN_FESpacesVelo, NSN_FctVelo, NSN_ParamFctVelo,
          NSN_FEValuesVelo,
          fesp, fefct,
          NSFctVelo,
          NSFEFctIndexVelo, NSFEMultiIndexVelo,
          NSN_ParamsVelo, NSBeginParamVelo);
      }
      else                       // Newton method
      {
        aux =  new TAuxParam3D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo,
          NSN_ParamFctVelo_GradVelo,
          NSN_FEValuesVelo_GradVelo,
          fesp, fefct,
          NSFctVelo_GradVelo,
          NSFEFctIndexVelo_GradVelo, NSFEMultiIndexVelo_GradVelo,
          NSN_ParamsVelo_GradVelo, NSBeginParamVelo_GradVelo);
      }

      // assemble
      Assemble3D(N_FESpaces, fesp,
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
            UpwindForNavierStokes3D(SQMATRICES[0], U1Array[k], U2Array[k], U3Array[k]);
            cout << "UPWINDING DONE : level " << k << endl;
            break;

          case 3:
          case 4:
            // do upwinding with two matrices
            UpwindForNavierStokes3D(SQMATRICES[0], U1Array[k], U2Array[k], U3Array[k]);
            UpwindForNavierStokes3D(SQMATRICES[4], U1Array[k], U2Array[k], U3Array[k]);
            UpwindForNavierStokes3D(SQMATRICES[8], U1Array[k], U2Array[k], U3Array[k]);
            cout << "UPWINDING DONE : level " << k << endl;
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

        // prepare everything for the assembling of slip with friction bc
        // on all levels
        N_FESpaces = 1;
        N_SquareMatrices = 9;
        N_RectMatrices = 3;
        N_Rhs = 3;
        DiscreteForm0 = NULL;

        SQMATRICES[0] = MatricesA11[k];
        SQMATRICES[1] = MatricesA22[k];
        SQMATRICES[2] = MatricesA33[k];
        SQMATRICES[3] = MatricesA12[k];
        SQMATRICES[4] = MatricesA13[k];
        SQMATRICES[5] = MatricesA21[k];
        SQMATRICES[6] = MatricesA23[k];
        SQMATRICES[7] = MatricesA31[k];
        SQMATRICES[8] = MatricesA32[k];

        MATRICES[0] = MatricesB1T[k];
        MATRICES[1] = MatricesB2T[k];
        MATRICES[2] = MatricesB3T[k];

        fesp[0] = USpaces[k];
        ferhs[0] = USpaces[k];
        ferhs[1] = USpaces[k];
        ferhs[2] = USpaces[k];

        RHSs[0] = RhsArray[k];
        RHSs[1] = RhsArray[k]+N_Uarray[k];
        RHSs[2] = RhsArray[k]+2*N_Uarray[k];

        Assemble3DSlipBC(N_FESpaces, fesp,
          N_SquareMatrices, SQMATRICES,
          N_RectMatrices, MATRICES,
          N_Rhs, RHSs, ferhs,
          DiscreteForm0,
          BoundaryConditions,
          BoundValues,
          aux);
	
	switch (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY)
	{
	    case 40000:
		SQMATRICES[0] = MatricesA11[k];
		SQMATRICES[1] = MatricesA12[k];
		SQMATRICES[2] = MatricesA13[k];
		SQMATRICES[3] = MatricesA21[k];
		SQMATRICES[4] = MatricesA22[k];
		SQMATRICES[5] = MatricesA23[k];
		SQMATRICES[6] = MatricesA31[k];
		SQMATRICES[7] = MatricesA32[k];
		SQMATRICES[8] = MatricesA33[k];
		MATRICES[0] = MatricesB1T[k];
		MATRICES[1] = MatricesB2T[k];
		MATRICES[2] = MatricesB3T[k];
		ModifyMatrixSlipBC(SQMATRICES,MATRICES,N_U, RhsArray[k]);
		break;
	}
        // reset MATRICES for solver
        SQMATRICES[0] = MatricesA11[k];
        SQMATRICES[1] = MatricesA12[k];
        SQMATRICES[2] = MatricesA13[k];
        SQMATRICES[3] = MatricesA21[k];
        SQMATRICES[4] = MatricesA22[k];
        SQMATRICES[5] = MatricesA23[k];
        SQMATRICES[6] = MatricesA31[k];
        SQMATRICES[7] = MatricesA32[k];
        SQMATRICES[8] = MatricesA33[k];
        MATRICES[0] = MatricesB1[k];
        MATRICES[1] = MatricesB2[k];
        MATRICES[2] = MatricesB3[k];
        MATRICES[3] = MatricesB1T[k];
        MATRICES[4] = MatricesB2T[k];
        MATRICES[5] = MatricesB3T[k];
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
          aux = new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
        }
        else
        {
          N_FESpaces = 2;
          fesp[1] = USpaces[mg_level];
          fefct[0] = U1Array[mg_level];
          fefct[1] = U2Array[mg_level];
          fefct[2] = U3Array[mg_level];
          aux =  new TAuxParam3D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo,
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
        Assemble3D(N_FESpaces, fesp,
          N_SquareMatrices, SQMATRICES,
          0, NULL,
          N_Rhs, RHSs, ferhs,
          DiscreteFormAuxProbPressSep,
          BoundaryConditionsPressureSeparation,
          BoundaryValuesPressureSeparation,
          aux);

        // solve linear system
        //  Solver(sqmatrixPressSep, RHSs[0],separated_pressure_array ,1);
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
        //   this has to be said to the assembling routine by an aux object
        fesp[0] = USpaces[mg_level];
        fesp[1] = pressure_separation_space;

        fefct[0] = separated_pressure_fe_funct;

        aux =  new TAuxParam3D(NSN_FESpacesPressSep, NSN_FctPressSep, NSN_ParamFctPressSep,
          NSN_FEValuesPressSep,
          fesp+1, fefct,
          NSFctPressSep,
          NSFEFctIndexPressSep, NSFEMultiIndexPressSep,
          NSN_ParamsPressSep, NSBeginParamPressSep);

        // assemble the right hand side
        N_FESpaces = 2;
        N_SquareMatrices = 0;
        N_RectMatrices = 0;
        N_Rhs = 3;
        RHSs[0] = rhs;
        RHSs[1] = rhs + N_U;
        RHSs[2] = rhs + 2*N_U;
        memset(rhs, 0, (3*N_U+N_P)*SizeOfDouble);
        ferhs[0] = USpaces[mg_level];
        ferhs[1] = USpaces[mg_level];
        ferhs[2] = USpaces[mg_level];
        DiscreteForm = DiscreteFormPressSep;

        Assemble3D(N_FESpaces, fesp,
          N_SquareMatrices, SQMATRICES,
          N_RectMatrices, MATRICES,
          N_Rhs, RHSs, ferhs,
          DiscreteForm,
          BoundaryConditions,
          BoundValues,
          aux);

        // initialize solution array for separated pressure
        memset(sol+3*N_U,0, N_P*SizeOfDouble);
      }

      if (DiscreteForm == DiscreteFormVMSProjection)
      {
	  SQMATRICES[0] = MatricesA[k];
	  SQMATRICES[1] = MatricesL[k];
	  MATRICES[3] = Matrices_tilde_G11[k];
	  MATRICES[4] = Matrices_tilde_G22[k];
	  MATRICES[5] = Matrices_tilde_G33[k];
	  MATRICES[6] = Matrices_G11[k];
	  MATRICES[7] = Matrices_G22[k];
	  MATRICES[8] = Matrices_G33[k];

	  LumpMassMatrixToDiag(MatricesL[k]);
	  VMS_ProjectionUpdateMatrices(N_Uarray[k],USpaces[k]->GetActiveBound(),
	  			      ProjectionSpaces[k]->GetN_DegreesOfFreedom(),
	  			      SQMATRICES,MATRICES);
      }
    }                            // endfor, assembling done
    
   // modify rhs for RFB stabilization
    if (TDatabase::ParamDB->DISCTYPE==NSE_RFB)
    {
	ApproximateRFBSolutionQuadNSE3D(coll, U1Array[mg_level], U2Array[mg_level],
					U3Array[mg_level], Coefficients[0], rhs);
	OutPut("RFB DONE"<<endl);
    }

    t2 = GetTime();
    assemble_time += (t2-t1);
    OutPut("time for assembling " << t2-t1 << " total on this level " << assemble_time << endl);
    coll->GetHminHmax(&hmin,&hmax);
    OutPut("h_min : " << hmin << " h_max : " << hmax << endl);

    // set rhs for Dirichlet nodes
    memcpy(sol+N_Active, rhs_high+N_Active, N_NonActive*SizeOfDouble);
    memcpy(sol+N_U+N_Active, rhs_high+N_U+N_Active, N_NonActive*SizeOfDouble);
    memcpy(sol+2*N_U+N_Active, rhs_high+2*N_U+N_Active, N_NonActive*SizeOfDouble);

    // compute defect
    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
      IntoL20FEFunction3D(sol+3*N_U, N_P,PSpaces[mg_level]);
    defect = new double[N_Unknowns];
    memset(defect,0,N_Unknowns*SizeOfDouble);
    Defect(sqmatrices,matrices,sol,rhs_high,defect);

    // for (k=0;k<N_Unknowns; k++)
    //  if ( abs(defect[k]) > 10)
    //    OutPut(k << " def " << defect[k] << " rhs " << rhs[k] << " sol " << sol[k] << endl);

    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
      IntoL20Vector3D(defect+3*N_U, N_P,pressure_space_code);
    residual =  Ddot(N_Unknowns,defect,defect);
    impuls_residual = Ddot(3*N_U,defect,defect);
    if (pressure_separation == 1)
      OutPut("p sep ");
    OutPut("nonlinear iteration step   0");
    OutPut(setw(14) << impuls_residual);
    OutPut(setw(14) << residual-impuls_residual);
    OutPut(setw(14) << sqrt(residual) << endl);

    /*if ((sqrt(residual)<=limit))
    {
      j = 0;
      OutPut("ITE : " << setw(3) << j);
      OutPut(" (" << setw(3) << N_LinIter << " LINITE)");
      OutPut("  TIME FOR SOLVER : " << solver_time << "s");
      OutPut("  RES : " <<  sqrt(residual));
      OutPut(endl);
      continue;
      }*/
    // solve system

    switch(TDatabase::ParamDB->SOLVER_TYPE)
    {
      case AMG:
        TDatabase::ParamDB->SC_VERBOSE=1;
        TDatabase::ParamDB->CC_VERBOSE=1;
        switch(TDatabase::ParamDB->NSTYPE)
        {
          case 1:
            Solver(sqmatrixA, matrixB1, matrixB2, matrixB3,
              rhs_high, sol);
            break;

          default:
            OutPut("AMG not implemented "<< endl);
            exit(4711);
            break;
        }
        /*
                case 2:
                  Solver(sqmatrixA, matrixB1T, matrixB2T,
                         matrixB1, matrixB2, rhs, sol);
                break;

                case 3:
                  Solver(sqmatrixA11, sqmatrixA12, sqmatrixA21,
                         sqmatrixA22, matrixB1, matrixB2, rhs, sol);
                break;

                case 4:
                  Solver(sqmatrixA11, sqmatrixA12, sqmatrixA21,
                         sqmatrixA22, matrixB1T, matrixB2T,
                         matrixB1, matrixB2, rhs, sol);
                break;
                }*/
        break;

      case DIRECT:
	  t1 = GetTime();
        switch(TDatabase::ParamDB->NSTYPE)
        {
          case 2:
	      DirectSolver(sqmatrixA, matrixB1T, matrixB2T, matrixB3T,
		     matrixB1, matrixB2, matrixB3, rhs_high, sol);
	      break;
          case 4:
	      
	      /*DirectSolver(sqmatrixA11, sqmatrixA12, sqmatrixA13, 
			   sqmatrixA21, sqmatrixA22, sqmatrixA23, 
			   sqmatrixA31, sqmatrixA32, sqmatrixA33, 
			   matrixB1T, matrixB2T, matrixB3T,
			   matrixB1, matrixB2, matrixB3, rhs_high, sol);*/
	      OutPut("Da stimmt noch etwas nicht !!!"<<endl);
	      exit(1);
	      break;

          default:
            OutPut("Direct solver not implemented !!!"<< endl);
            exit(4711);
            break;
        }
        t2 = GetTime();
        solver_time += (t2-t1);
	break;
      case GMG:
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
	OutPut("norm sol " << Ddot(3*N_U+N_P,sol,sol) << endl);
        for (ll=0;ll<calculations;ll++)
        {
          // for two level method of Stokes equation
          if(TDatabase::ParamDB->P9 == 123456789)
          {
            fesp[0] = USpaces[mg_level];
            fefct[0] = U1Array[mg_level];
            fefct[1] = U2Array[mg_level];
            fefct[2] = U3Array[mg_level];
            aux =  new TAuxParam3D(NSN_FESpacesVelo, NSN_FctVelo,
              NSN_ParamFctVelo,
              NSN_FEValuesVelo,
              fesp, fefct,
              NSFctVelo,
              NSFEFctIndexVelo, NSFEMultiIndexVelo,
              NSN_ParamsVelo, NSBeginParamVelo);

            // errors in first velocity component
            U1Array[mg_level]->GetErrors(ExactU1, 4, NSAllDerivatives, 2,
              L2H1Errors,
              NULL, aux, 1, USpaces+mg_level, errors_mg);
            p1 = errors_mg[0] * errors_mg[0];

            // errors in second velocity component
            U2Array[mg_level]->GetErrors(ExactU2, 4, NSAllDerivatives, 2,
              L2H1Errors,
              NULL, aux, 1, USpaces+mg_level, errors_mg);
            p1 += errors_mg[0] * errors_mg[0];
            // errors in second velocity component
            U3Array[mg_level]->GetErrors(ExactU3, 4, NSAllDerivatives, 2,
              L2H1Errors,
              NULL, aux, 1, USpaces+mg_level, errors_mg);
            p1 += errors_mg[0] * errors_mg[0];

            // errors in first velocity component
            soldiff_fe1->GetErrors(ExactNull, 4, NSAllDerivatives, 2,
              L2H1Errors,
              NULL, aux, 1, USpaces+mg_level, errors_mg);
            p1 = errors_mg[0] * errors_mg[0];

            // errors in second velocity component
            soldiff_fe2->GetErrors(ExactNull, 4, NSAllDerivatives, 2,
              L2H1Errors,
              NULL, aux, 1, USpaces+mg_level, errors_mg);
            p1 += errors_mg[0] * errors_mg[0];

            soldiff_fe3->GetErrors(ExactNull, 4, NSAllDerivatives, 2,
              L2H1Errors,
              NULL, aux, 1, USpaces+mg_level, errors_mg);
            p1 += errors_mg[0] * errors_mg[0];
            firsterror =  p1 = sqrt(p1);
            delete aux;

            for(l=0;l<N_Unknowns;l++)
              soldiff[l] = sol[l]-fesol[l];

            p3 = sqrt(Ddot(3*N_U,soldiff,soldiff));
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
            fefct[2] = U3Array[mg_level];
            aux =  new TAuxParam3D(NSN_FESpacesVelo, NSN_FctVelo,
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
            // errors in second velocity component
            U3Array[mg_level]->GetErrors(ExactU3, 3, NSAllDerivatives, 2,
              L2H1Errors,
              NULL, aux, 1, USpaces+mg_level, errors_mg);
            // errors in first velocity component
            soldiff_fe1->GetErrors(ExactNull, 3, NSAllDerivatives, 2,
              L2H1Errors,
              NULL, aux, 1, USpaces+mg_level, errors_mg);
            p1 = errors_mg[0] * errors_mg[0];

            // errors in second velocity component
            soldiff_fe2->GetErrors(ExactNull, 3, NSAllDerivatives, 2,
              L2H1Errors,
              NULL, aux, 1, USpaces+mg_level, errors_mg);

            p1 += errors_mg[0] * errors_mg[0];
            soldiff_fe3->GetErrors(ExactNull, 3, NSAllDerivatives, 2,
              L2H1Errors,
              NULL, aux, 1, USpaces+mg_level, errors_mg);

            p1 += errors_mg[0] * errors_mg[0];
            p1 = sqrt(p1);
            delete aux;

            for(l=0;l<N_Unknowns;l++)
              soldiff[l] = sol[l]-fesol[l];
            p3 = sqrt(Ddot(3*N_U,soldiff,soldiff));

            p2 = p1;
            p4 = p3;
            lasterror = p1;
            lasterrorl2 = p3;
          }                      // endif MEASURE_ERRORS

          for(l=0;l<N_Unknowns;l++)
          {
            p2 = sol[l]-oldsol[l];
            if ((TDatabase::ParamDB->P9 == 1.2389) && (i==FirstSolve))
              p2 = 0.0;
            sol[l] = oldsol[l] + omega * p2;
          }

          if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
            IntoL20FEFunction3D(sol+3*N_U, N_P,PSpaces[mg_level]);

          if(TDatabase::ParamDB->P9 == 123456789)
          {
            if (!pre_calculation)
            {
              OutPut("average error reduction rate (L2/l2) " << pow(lasterror/firsterror,1.0/N_LinIter));
              OutPut(" " << pow(lasterrorl2/firsterrorl2,1.0/N_LinIter) << endl);
            }
            if (pre_calculation)
            {
              for (k=0;k<3*N_U+N_P;k++)
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

    OutPut("MEMORY: " << setw(10) << GetMemory()/(1048576.0));
    OutPut(" MB" << endl);

    // ************************************************************* //
    // the nonlinear iteration
    // ************************************************************* //

    for(j=1;j<=Max_It;j++)
    {
      if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
        MG->RestrictToAllGrids();

      memcpy(oldsol, sol, SizeOfDouble*N_Unknowns);

      t1 = GetTime();
      for(k=low;k<=mg_level;k++)
      {
        fesp[0] = USpaces[k];
        fesp[1] = PSpaces[k];

        fefct[0] = U1Array[k];
        fefct[1] = U2Array[k];
        fefct[2] = U3Array[k];

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

            case SDFEM_DIVDIV:
              DiscreteForm = DiscreteFormNLSDFEM_DivDiv;
              break;
          }                      // endswitch
        if (DiscreteForm == NULL)
	{
	    OutPut("Discrete form not implemented " << endl);
	    exit(4711);
	}
        // this can only happen if pressure separation is applied
        if (TDatabase::ParamDB->STOKES_PROBLEM)
          DiscreteForm = DiscreteFormNLUpwind;

        switch (TDatabase::ParamDB->NSTYPE)
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
              N_Rhs = 3;
              rhs = RhsArray[k];
              RHSs[0] = rhs;
              RHSs[1] = rhs + N_Uarray[k];
              RHSs[2] = rhs + 2*N_Uarray[k];
              memset(rhs, 0, (3*N_Uarray[k]+N_Parray[k])*SizeOfDouble);
              N_FESpaces = 2;
              ferhs[0] = USpaces[k];
              ferhs[1] = USpaces[k];
              ferhs[2] = USpaces[k];
            }
            if (DiscreteForm == DiscreteFormNLVMSProjection)
            {
              MATRICES[0] =  Matrices_tilde_G11[k];
              MATRICES[1] =  Matrices_tilde_G22[k];
              MATRICES[2] =  Matrices_tilde_G33[k];

              MATRICES[0]->Reset();
              MATRICES[1]->Reset();
              MATRICES[2]->Reset();
	      N_RectMatrices = 3;

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
              N_RectMatrices = 3;
              MATRICES[0] = MatricesB1T[k];
              MATRICES[1] = MatricesB2T[k];
              MATRICES[2] = MatricesB3T[k];

              MATRICES[0]->Reset();
              MATRICES[1]->Reset();
              MATRICES[2]->Reset();

              N_Rhs = 3;
              rhs = RhsArray[k];
              RHSs[0] = rhs;
              RHSs[1] = rhs + N_Uarray[k];
              RHSs[2] = rhs + 2*N_Uarray[k];
              memset(rhs, 0, (3*N_Uarray[k]+N_Parray[k])*SizeOfDouble);
              N_FESpaces = 2;
              ferhs[0] = USpaces[k];
              ferhs[1] = USpaces[k];
              ferhs[2] = USpaces[k];
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
		if (DiscreteForm == DiscreteFormNLSmagorinsky)
		    {
                SQMATRICES[0] = MatricesA11[k];
                SQMATRICES[1] = MatricesA12[k];
                SQMATRICES[2] = MatricesA13[k];
                SQMATRICES[3] = MatricesA21[k];
                SQMATRICES[4] = MatricesA22[k];
                SQMATRICES[5] = MatricesA23[k];
                SQMATRICES[6] = MatricesA31[k];
                SQMATRICES[7] = MatricesA32[k];
                SQMATRICES[8] = MatricesA32[k];
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
                N_RectMatrices = 0;
                last_sq = 8;
              }
              else
              {
                SQMATRICES[0] = MatricesA11[k];
                SQMATRICES[1] = MatricesA22[k];
                SQMATRICES[2] = MatricesA33[k];
                SQMATRICES[0]->Reset();
                SQMATRICES[1]->Reset();
                SQMATRICES[2]->Reset();

                N_SquareMatrices = 3;
                N_RectMatrices = 0;
                last_sq = 2;
              }
              N_Rhs = 0;
              N_FESpaces = 1;
            }
            else                 // Newton method
            {
              SQMATRICES[0] = MatricesA11[k];
              SQMATRICES[1] = MatricesA12[k];
              SQMATRICES[2] = MatricesA13[k];
              SQMATRICES[3] = MatricesA21[k];
              SQMATRICES[4] = MatricesA22[k];
              SQMATRICES[5] = MatricesA23[k];
              SQMATRICES[6] = MatricesA31[k];
              SQMATRICES[7] = MatricesA32[k];
              SQMATRICES[8] = MatricesA33[k];
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
              N_RectMatrices = 0;
              last_sq = 8;

              N_Rhs = 3;
              N_FESpaces = 1;
              rhs = RhsArray[k];
              RHSs[0] = rhs;
              RHSs[1] = rhs + N_Uarray[k];
              RHSs[2] = rhs + 2*N_Uarray[k];
              memset(rhs, 0, (3*N_Uarray[k]+N_Parray[k])*SizeOfDouble);
              ferhs[0] = USpaces[k];
              ferhs[1] = USpaces[k];
              ferhs[2] = USpaces[k];
            }
            break;

          case 4:
            if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
            {
		if (DiscreteForm==DiscreteFormNLSDFEM)
		{
		    N_SquareMatrices = 3;
		    SQMATRICES[0] = MatricesA11[k];
		    SQMATRICES[1] = MatricesA22[k];
		    SQMATRICES[2] = MatricesA33[k];
		    SQMATRICES[0]->Reset();
		    SQMATRICES[1]->Reset();
		    SQMATRICES[2]->Reset();
		    
		    N_RectMatrices = 3;
		    MATRICES[0] = MatricesB1T[k];
		    MATRICES[1] = MatricesB2T[k];
		    MATRICES[2] = MatricesB3T[k];
		    MATRICES[0]->Reset();
		    MATRICES[1]->Reset();
		    MATRICES[2]->Reset();
		    
		    N_Rhs = 3;
		    rhs = RhsArray[k];
		    RHSs[0] = rhs;
		    RHSs[1] = rhs + N_Uarray[k];
		    RHSs[2] = rhs + 2*N_Uarray[k];
		    memset(rhs, 0, (3*N_Uarray[k]+N_Parray[k])*SizeOfDouble);
		    N_FESpaces = 2;
		    
		    ferhs[0] = USpaces[k];
		    ferhs[1] = USpaces[k];
		    ferhs[2] = USpaces[k];
		    last_sq = 2;
		}
		else
		{
		    if (DiscreteForm == DiscreteFormNLSmagorinsky)
		    {
			SQMATRICES[0] = MatricesA11[k];
			SQMATRICES[1] = MatricesA12[k];
			SQMATRICES[2] = MatricesA13[k];
			SQMATRICES[3] = MatricesA21[k];
			SQMATRICES[4] = MatricesA22[k];
			SQMATRICES[5] = MatricesA23[k];
			SQMATRICES[6] = MatricesA31[k];
			SQMATRICES[7] = MatricesA32[k];
			SQMATRICES[8] = MatricesA32[k];
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
			N_RectMatrices = 0;
			last_sq = 8;
			N_Rhs = 0;
			N_FESpaces = 1;
		    }
		    else
		    {
			if (DiscreteForm == DiscreteFormNLSDFEM_DivDiv)
			{
			    SQMATRICES[0] = MatricesA11[k];
			    SQMATRICES[1] = MatricesA12[k];
			    SQMATRICES[2] = MatricesA13[k];
			    SQMATRICES[3] = MatricesA21[k];
			    SQMATRICES[4] = MatricesA22[k];
			    SQMATRICES[5] = MatricesA23[k];
			    SQMATRICES[6] = MatricesA31[k];
			    SQMATRICES[7] = MatricesA32[k];
			    SQMATRICES[8] = MatricesA32[k];
			    SQMATRICES[0]->Reset();
			    SQMATRICES[1]->Reset();
			    SQMATRICES[2]->Reset();
			    SQMATRICES[3]->Reset();
			    SQMATRICES[4]->Reset();
			    SQMATRICES[5]->Reset();
			    SQMATRICES[6]->Reset();
			    SQMATRICES[7]->Reset();
			    SQMATRICES[8]->Reset();
			    
			    MATRICES[0] = MatricesB1T[k];
			    MATRICES[1] = MatricesB2T[k];
			    MATRICES[2] = MatricesB3T[k];
			    MATRICES[0]->Reset();
			    MATRICES[1]->Reset();
			    MATRICES[2]->Reset();
			    
			    N_SquareMatrices = 9;
			    N_RectMatrices = 3;
			    N_Rhs = 3;
			    rhs = RhsArray[k];
			    RHSs[0] = rhs;
			    RHSs[1] = rhs + N_Uarray[k];
			    RHSs[2] = rhs + 2*N_Uarray[k];
			    memset(rhs, 0, (3*N_Uarray[k]+N_Parray[k])*SizeOfDouble);
			    N_FESpaces = 2;
			    
			    ferhs[0] = USpaces[k];
			    ferhs[1] = USpaces[k];
			    ferhs[2] = USpaces[k];
			    last_sq = 8;
			}
			else //default
			{
			    N_SquareMatrices = 3;
			    SQMATRICES[0] = MatricesA11[k];
			    SQMATRICES[1] = MatricesA22[k];
			    SQMATRICES[2] = MatricesA33[k];
			    SQMATRICES[0]->Reset();
			    SQMATRICES[1]->Reset();
			    SQMATRICES[2]->Reset();
			    
			    N_RectMatrices = 0;
			    
			    N_Rhs = 0;
			    N_FESpaces = 1;
			    last_sq = 2;
			}
		    }
		}
	    }
	    else                 // Newton method
            {
              SQMATRICES[0] = MatricesA11[k];
              SQMATRICES[1] = MatricesA12[k];
              SQMATRICES[2] = MatricesA13[k];
              SQMATRICES[3] = MatricesA21[k];
              SQMATRICES[4] = MatricesA22[k];
              SQMATRICES[5] = MatricesA23[k];
              SQMATRICES[6] = MatricesA31[k];
              SQMATRICES[7] = MatricesA32[k];
              SQMATRICES[8] = MatricesA33[k];
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
              N_RectMatrices = 0;
              last_sq = 8;

              N_Rhs = 3;
              N_FESpaces = 1;
              rhs = RhsArray[k];
              RHSs[0] = rhs;
              RHSs[1] = rhs + N_Uarray[k];
              RHSs[2] = rhs + 2*N_Uarray[k];
              memset(rhs, 0, (3*N_Uarray[k]+N_Parray[k])*SizeOfDouble);
              ferhs[0] = USpaces[k];
              ferhs[1] = USpaces[k];
              ferhs[2] = USpaces[k];
            }
            break;
        }                        // endswitch

        if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
          if ((DiscreteForm == DiscreteFormNLSmagorinsky)||
	      (DiscreteForm == DiscreteFormNLVMSProjection))
        {
          aux =  new TAuxParam3D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo,
            NSN_ParamFctVelo_GradVelo,
            NSN_FEValuesVelo_GradVelo,
            fesp, fefct,
            NSFctVelo_GradVelo,
            NSFEFctIndexVelo_GradVelo, NSFEMultiIndexVelo_GradVelo,
            NSN_ParamsVelo_GradVelo, NSBeginParamVelo_GradVelo);
        }
        else
        {
          aux =  new TAuxParam3D(NSN_FESpacesVelo, NSN_FctVelo, NSN_ParamFctVelo,
            NSN_FEValuesVelo,
            fesp, fefct,
            NSFctVelo,
            NSFEFctIndexVelo, NSFEMultiIndexVelo,
            NSN_ParamsVelo, NSBeginParamVelo);
        }
        else                     // Newton method
        {
          aux =  new TAuxParam3D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo,
            NSN_ParamFctVelo_GradVelo,
            NSN_FEValuesVelo_GradVelo,
            fesp, fefct,
            NSFctVelo_GradVelo,
            NSFEFctIndexVelo_GradVelo, NSFEMultiIndexVelo_GradVelo,
            NSN_ParamsVelo_GradVelo, NSBeginParamVelo_GradVelo);
        }

        // assembling
        Assemble3D(N_FESpaces, fesp,
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
              UpwindForNavierStokes3D(SQMATRICES[0], U1Array[k], U2Array[k], U3Array[k]);
              cout << "UPWINDING DONE : level " << k << endl;
              break;

            case 3:
            case 4:
              // do upwinding with three matrices
              if (N_SquareMatrices == 9)
              {
                UpwindForNavierStokes3D(SQMATRICES[0], U1Array[k], U2Array[k], U3Array[k]);
                UpwindForNavierStokes3D(SQMATRICES[4], U1Array[k], U2Array[k], U3Array[k]);
                UpwindForNavierStokes3D(SQMATRICES[8], U1Array[k], U2Array[k], U3Array[k]);
              }
              else
              {
                UpwindForNavierStokes3D(SQMATRICES[0], U1Array[k], U2Array[k], U3Array[k]);
                UpwindForNavierStokes3D(SQMATRICES[1], U1Array[k], U2Array[k], U3Array[k]);
                UpwindForNavierStokes3D(SQMATRICES[2], U1Array[k], U2Array[k], U3Array[k]);
              }
              cout << "UPWINDING DONE(2) : level " << k << endl;
              break;
          }                      // endswitch
        }                        // endif

        if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
        {
          // prepare everything for the assembling of slip with friction bc
          // on all levels
          N_FESpaces = 1;
          N_SquareMatrices = 9;
          N_RectMatrices = 0;
          N_Rhs = 3;
          DiscreteForm0 = NULL;

          SQMATRICES[0] = MatricesA11[k];
          SQMATRICES[1] = MatricesA22[k];
          SQMATRICES[2] = MatricesA33[k];
          SQMATRICES[3] = MatricesA12[k];
          SQMATRICES[4] = MatricesA13[k];
          SQMATRICES[5] = MatricesA21[k];
          SQMATRICES[6] = MatricesA23[k];
          SQMATRICES[7] = MatricesA31[k];
          SQMATRICES[8] = MatricesA32[k];

          fesp[0] = USpaces[k];
          ferhs[0] = USpaces[k];
          ferhs[1] = USpaces[k];
          ferhs[2] = USpaces[k];

          RHSs[0] = RhsArray[k];
          RHSs[1] = RhsArray[k]+N_Uarray[k];
          RHSs[2] = RhsArray[k]+2*N_Uarray[k];

          Assemble3DSlipBC(N_FESpaces, fesp,
            N_SquareMatrices, SQMATRICES,
            N_RectMatrices, MATRICES,
            N_Rhs, RHSs, ferhs,
            DiscreteForm0,
            BoundaryConditions,
            BoundValues,
            aux);

	switch (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY)
	{
	    case 40000:

        SQMATRICES[0] = MatricesA11[k];
        SQMATRICES[1] = MatricesA12[k];
        SQMATRICES[2] = MatricesA13[k];
        SQMATRICES[3] = MatricesA21[k];
        SQMATRICES[4] = MatricesA22[k];
        SQMATRICES[5] = MatricesA23[k];
        SQMATRICES[6] = MatricesA31[k];
        SQMATRICES[7] = MatricesA32[k];
        SQMATRICES[8] = MatricesA33[k];

	MATRICES[0] = MatricesB1T[k];
	MATRICES[1] = MatricesB2T[k];
	MATRICES[2] = MatricesB3T[k];
	ModifyMatrixSlipBC(SQMATRICES,MATRICES,N_Uarray[k], RhsArray[k]);
	break;
	}
        }
        delete aux;
	if (DiscreteForm == DiscreteFormNLVMSProjection)
	{
	    SQMATRICES[0] = MatricesA[k];
	    SQMATRICES[1] = MatricesL[k];
	    MATRICES[3] = Matrices_tilde_G11[k];
	    MATRICES[4] = Matrices_tilde_G22[k];
	    MATRICES[5] = Matrices_tilde_G33[k];
	    MATRICES[6] = Matrices_G11[k];
	    MATRICES[7] = Matrices_G22[k];
	    MATRICES[8] = Matrices_G33[k];

	    LumpMassMatrixToDiag(MatricesL[k]);
	    VMS_ProjectionUpdateMatrices(N_Uarray[k],USpaces[k]->GetActiveBound(),
					ProjectionSpaces[k]->GetN_DegreesOfFreedom(),
					SQMATRICES,MATRICES);
	}
	
      }                          // endfor k, assembling done
      // RFB stabilization
      if (TDatabase::ParamDB->DISCTYPE==NSE_RFB)
      {
	  // assemble rhs
	  fesp[0] = USpaces[mg_level];
          aux = new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      
	  // assemble the right hand side
	  N_FESpaces = 1;
	  N_SquareMatrices = 0;
	  N_RectMatrices = 0;
	  N_Rhs = 3;
	  RHSs[0] = rhs;
	  RHSs[1] = rhs + N_U;
	  RHSs[2] = rhs + 2* N_U;
	  memset(rhs, 0, (3*N_U+N_P)*SizeOfDouble);
	  ferhs[0] = USpaces[mg_level];
	  ferhs[1] = USpaces[mg_level];
	  ferhs[2] = USpaces[mg_level];
	  DiscreteForm = DiscreteFormNSRFBRhs;

	  Assemble3D(N_FESpaces, fesp,
		     N_SquareMatrices, SQMATRICES,
		     N_RectMatrices, MATRICES,
		     N_Rhs, RHSs, ferhs,
		     DiscreteForm,
		     BoundaryConditions,
		     BoundValues,
		     aux);
	  delete aux;
	  // modify rhs with RFB stabilization
	  ApproximateRFBSolutionQuadNSE3D(coll, U1Array[mg_level], U2Array[mg_level],
					  U3Array[mg_level], Coefficients[0], rhs);
	  OutPut("RFB DONE"<<endl);
      }
      t2 = GetTime();
      assemble_time += (t2-t1);
      OutPut("time for assembling " << t2-t1 << " total on this level " << assemble_time << endl);

      if ((TDatabase::ParamDB->PRESSURE_SEPARATION==2)&&(j==1))
      {
        OutPut("apply pressure separation"<<endl);
        // save original pressure
        if ((j==1)||(1))
          memcpy(nosep_p,sol+3*N_U, N_P*SizeOfDouble);
        else
        {
          Daxpy(N_P, 1.0, sol+3*N_U, nosep_p);
          memcpy(sol+3*N_U, nosep_p, N_P*SizeOfDouble);
        }
        // assemble separated pressure entry on right hand side
        // first : interpolate discrete normal pressure to the
        //         pressure separation space
        Prolongate(pressure_space, pressure_separation_space,
          sol+3*N_U, separated_pressure_array,
          separated_pressure_aux_array);

        // second : assemble
        // the gradient of the separated pressure is needed for assembling
        // this has to be said to the assembling routine by an aux object
        fesp[0] = USpaces[mg_level];
        fesp[1] = pressure_separation_space;

        fefct[0] = separated_pressure_fe_funct;

        aux =  new TAuxParam3D(NSN_FESpacesPressSep, NSN_FctPressSep, NSN_ParamFctPressSep,
          NSN_FEValuesPressSep,
          fesp+1, fefct,
          NSFctPressSep,
          NSFEFctIndexPressSep, NSFEMultiIndexPressSep,
          NSN_ParamsPressSep, NSBeginParamPressSep);

        // assemble the right hand side
        N_FESpaces = 2;
        N_SquareMatrices = 0;
        N_RectMatrices = 0;
        N_Rhs = 3;
        RHSs[0] = rhs;
        RHSs[1] = rhs + N_U;
        RHSs[2] = rhs + 2*N_U;
        memset(rhs, 0, (3*N_U+N_P)*SizeOfDouble);
        ferhs[0] = USpaces[mg_level];
        ferhs[1] = USpaces[mg_level];
        ferhs[2] = USpaces[mg_level];
        DiscreteForm = DiscreteFormPressSep;

        Assemble3D(N_FESpaces, fesp,
          N_SquareMatrices, SQMATRICES,
          N_RectMatrices, MATRICES,
          N_Rhs, RHSs, ferhs,
          DiscreteForm,
          BoundaryConditions,
          BoundValues,
          aux);

        // initialize solution array for separated pressure
        memset(sol+3*N_U,0, N_P*SizeOfDouble);
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
        fefct[2] = U3Array[mg_level];
        aux =  new TAuxParam3D(NSN_FESpacesVelo_GradVelo, NSN_FctVelo_GradVelo,
          NSN_ParamFctVelo_GradVelo,
          NSN_FEValuesVelo_GradVelo,
          fesp+1, fefct,
          NSFctVelo_GradVelo,
          NSFEFctIndexVelo_GradVelo, NSFEMultiIndexVelo_GradVelo,
          NSN_ParamsVelo_GradVelo, NSBeginParamVelo_GradVelo);
        RHSs[0] = rhsPressSep;
        memset(rhsPressSep, 0, N_P_sep*SizeOfDouble);

        // assemble
        Assemble3D(N_FESpaces, fesp,
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
        //   this has to be said to the assembling routine by an aux object
        fesp[0] = USpaces[mg_level];
        fesp[1] = pressure_separation_space;

        fefct[0] = separated_pressure_fe_funct;

        aux =  new TAuxParam3D(NSN_FESpacesPressSep, NSN_FctPressSep, NSN_ParamFctPressSep,
          NSN_FEValuesPressSep,
          fesp+1, fefct,
          NSFctPressSep,
          NSFEFctIndexPressSep, NSFEMultiIndexPressSep,
          NSN_ParamsPressSep, NSBeginParamPressSep);

        // assemble the right hand side
        N_FESpaces = 2;
        N_SquareMatrices = 0;
        N_RectMatrices = 0;
        N_Rhs = 3;
        RHSs[0] = rhs;
        RHSs[1] = rhs + N_U;
        RHSs[2] = rhs + 2*N_U;
        memset(rhs, 0, (3*N_U+N_P)*SizeOfDouble);
        ferhs[0] = USpaces[mg_level];
        ferhs[1] = USpaces[mg_level];
        ferhs[2] = USpaces[mg_level];
        DiscreteForm = DiscreteFormPressSep;

        Assemble3D(N_FESpaces, fesp,
          N_SquareMatrices, SQMATRICES,
          N_RectMatrices, MATRICES,
          N_Rhs, RHSs, ferhs,
          DiscreteForm,
          BoundaryConditions,
          BoundValues,
          aux);

        // initialize solution array for separated pressure
        memset(sol+3*N_U,0, N_P*SizeOfDouble);
      }

      // end of assembling

      // reset MATRICES for solver
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          SQMATRICES[0] = MatricesA[mg_level];
          MATRICES[0] = MatricesB1[mg_level];
          MATRICES[1] = MatricesB2[mg_level];
          MATRICES[2] = MatricesB3[mg_level];
          break;
        case 2:
          SQMATRICES[0] = MatricesA[mg_level];
          MATRICES[0] = MatricesB1[mg_level];
          MATRICES[1] = MatricesB2[mg_level];
          MATRICES[2] = MatricesB3[mg_level];
          MATRICES[3] = MatricesB1T[mg_level];
          MATRICES[4] = MatricesB2T[mg_level];
          MATRICES[5] = MatricesB3T[mg_level];
          break;
        case 3:
          SQMATRICES[0] = MatricesA11[mg_level];
          SQMATRICES[1] = MatricesA12[mg_level];
          SQMATRICES[2] = MatricesA13[mg_level];
          SQMATRICES[3] = MatricesA21[mg_level];
          SQMATRICES[4] = MatricesA22[mg_level];
          SQMATRICES[5] = MatricesA23[mg_level];
          SQMATRICES[6] = MatricesA31[mg_level];
          SQMATRICES[7] = MatricesA32[mg_level];
          SQMATRICES[8] = MatricesA33[mg_level];
          MATRICES[0] = MatricesB1[mg_level];
          MATRICES[1] = MatricesB2[mg_level];
          MATRICES[2] = MatricesB3[mg_level];
          break;
        case 4:
          SQMATRICES[0] = MatricesA11[mg_level];
          SQMATRICES[1] = MatricesA12[mg_level];
          SQMATRICES[2] = MatricesA13[mg_level];
          SQMATRICES[3] = MatricesA21[mg_level];
          SQMATRICES[4] = MatricesA22[mg_level];
          SQMATRICES[5] = MatricesA23[mg_level];
          SQMATRICES[6] = MatricesA31[mg_level];
          SQMATRICES[7] = MatricesA32[mg_level];
          SQMATRICES[8] = MatricesA33[mg_level];
          MATRICES[0] = MatricesB1[mg_level];
          MATRICES[1] = MatricesB2[mg_level];
          MATRICES[2] = MatricesB3[mg_level];
          MATRICES[3] = MatricesB1T[mg_level];
          MATRICES[4] = MatricesB2T[mg_level];
          MATRICES[5] = MatricesB3T[mg_level];
          break;
      }

      // set rhs for Dirichlet nodes
      memcpy(sol+N_Active, rhs_high+N_Active, N_NonActive*SizeOfDouble);
      memcpy(sol+N_U+N_Active, rhs_high+N_U+N_Active, N_NonActive*SizeOfDouble);
      memcpy(sol+2*N_U+N_Active, rhs_high+2*N_U+N_Active, N_NonActive*SizeOfDouble);

      // compute defect
      memset(defect,0,N_Unknowns*SizeOfDouble);
      Defect(sqmatrices,matrices,sol,rhs_high,defect);
      if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
        IntoL20Vector3D(defect+3*N_U, N_P,pressure_space_code);
      residual =  Ddot(N_Unknowns,defect,defect);
      impuls_residual = Ddot(3*N_U,defect,defect);
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
        {
          total_time2 = GetTime();
          OutPut("time on geometry level " << i << " : " << total_time2-total_time1<< " sec" << endl);
          break;
        }

        if (((TDatabase::ParamDB->PRESSURE_SEPARATION==1)||
          (TDatabase::ParamDB->PRESSURE_SEPARATION==2))&&
          (pressure_separation == 1))
        {
          // compute pressure
          memcpy(sol+3*N_U, nosep_p, N_P*SizeOfDouble);
          //Daxpy(N_P,1.0,nosep_p,sol+2*N_U);
          break;
        }
        if ((TDatabase::ParamDB->PRESSURE_SEPARATION==3)||
          (TDatabase::ParamDB->PRESSURE_SEPARATION==4)||
          (TDatabase::ParamDB->PRESSURE_SEPARATION==5))
        {
          Daxpy(N_P,1.0,nosep_p,sol+3*N_U);
          break;
        }
        // save original pressure
        memcpy(nosep_p,sol+3*N_U, N_P*SizeOfDouble);

        // assemble separated pressure entry on right hand side
        // first : interpolate discrete normal pressure to the
        //         pressure separation space
        Prolongate(pressure_space, pressure_separation_space,
          sol+3*N_U, separated_pressure_array,
          separated_pressure_aux_array);

        // TESTING: USE INTERPOLATION
        // separated_pressure_fe_funct->Interpolate(ExactP);
        // second : assemble
        // the gradient of the separated pressure is needed for assembling
        // this has to be said to the assembling routine by an aux object
        fesp[0] = USpaces[mg_level];
        fesp[1] = pressure_separation_space;

        fefct[0] = separated_pressure_fe_funct;

        aux =  new TAuxParam3D(NSN_FESpacesPressSep, NSN_FctPressSep, NSN_ParamFctPressSep,
          NSN_FEValuesPressSep,
          fesp+1, fefct,
          NSFctPressSep,
          NSFEFctIndexPressSep, NSFEMultiIndexPressSep,
          NSN_ParamsPressSep, NSBeginParamPressSep);

        // assemble the right hand side
        N_FESpaces = 2;
        N_SquareMatrices = 0;
        N_RectMatrices = 0;
        N_Rhs = 3;
        RHSs[0] = rhs;
        RHSs[1] = rhs + N_U;
        RHSs[2] = rhs + 2*N_U;
        memset(rhs, 0, (3*N_U+N_P)*SizeOfDouble);
        ferhs[0] = USpaces[mg_level];
        ferhs[1] = USpaces[mg_level];
        ferhs[2] = USpaces[mg_level];
        DiscreteForm = DiscreteFormPressSep;

        Assemble3D(N_FESpaces, fesp,
          N_SquareMatrices, SQMATRICES,
          N_RectMatrices, MATRICES,
          N_Rhs, RHSs, ferhs,
          DiscreteForm,
          BoundaryConditions,
          BoundValues,
          aux);

        // initialize solution array for separated pressure
        memset(sol+3*N_U,0, N_P*SizeOfDouble);
        // compute defect
        memset(defect,0,N_Unknowns*SizeOfDouble);
        Defect(sqmatrices,matrices,sol,rhs_high,defect);
        if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
          IntoL20Vector3D(defect+3*N_U, N_P,pressure_space_code);
        residual =  Ddot(N_Unknowns,defect,defect);
        impuls_residual = Ddot(3*N_U,defect,defect);
        // reset counter
        j = 0;
        pressure_separation = 1;
        OutPut("p sep ");
        OutPut("nonlinear iteration step " << setw(3) << j);
        OutPut(setw(14) << impuls_residual);
        OutPut(setw(14) << residual-impuls_residual);
        OutPut(setw(14) << sqrt(residual) << endl);
      }

      switch(TDatabase::ParamDB->SOLVER_TYPE)
      {
        case AMG:
          TDatabase::ParamDB->SC_VERBOSE_AMG=0;
          TDatabase::ParamDB->CC_VERBOSE=0;
          t1 = GetTime();
          switch(TDatabase::ParamDB->NSTYPE)
          {
            case 1:
              Solver(sqmatrixA, matrixB1, matrixB2, matrixB3,
                rhs_high, sol);
              break;
          }
          /*
                  case 2:
                    Solver(sqmatrixA, matrixB1T, matrixB2T,
                           matrixB1, matrixB2, rhs, sol);
                    break;

                  case 3:
                    Solver(sqmatrixA11, sqmatrixA12, sqmatrixA21,
                           sqmatrixA22, matrixB1, matrixB2, rhs, sol);
                    break;

                  case 4:
                    Solver(sqmatrixA11, sqmatrixA12, sqmatrixA21,
                           sqmatrixA22, matrixB1T, matrixB2T,
                           matrixB1, matrixB2, rhs, sol);
                    break;
                    }*/
          t2 = GetTime();
          solver_time += (t2-t1);
          break;

        case DIRECT:
	  t1 = GetTime();
        switch(TDatabase::ParamDB->NSTYPE)
        {
          case 2:
	      DirectSolver(sqmatrixA, matrixB1T, matrixB2T, matrixB3T,
		     matrixB1, matrixB2, matrixB3, rhs_high, sol);
	      break;
          case 4:
	      /*DirectSolver(MatricesA11[mg_level], MatricesA12[mg_level], MatricesA13[mg_level],
			   MatricesA21[mg_level], MatricesA22[mg_level], MatricesA23[mg_level],
			   MatricesA31[mg_level], MatricesA32[mg_level], MatricesA33[mg_level],
			   MatricesB1T[mg_level], MatricesB2T[mg_level], MatricesB3T[mg_level],
			   MatricesB1[mg_level], MatricesB2[mg_level], MatricesB3[mg_level],
			   rhs_high, sol);*/
	      OutPut("Da stimmt noch etwas nicht !!!"<<endl);
	      exit(1);
	      break;
          default:
            OutPut("Direct solver not implemented !!!"<< endl);
            exit(4711);
            break;
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
          break;
      }
      for(l=0;l<N_Unknowns;l++)
      {
	  p2 = sol[l]-oldsol[l];
	  sol[l] = oldsol[l] + omega * p2;
      }


/*      #ifdef __BENCH__
      // compute characteristic values (deltaP, Cd, Cl)
      //   GetCdCl(U1Array[mg_level], U2Array[mg_level], U3Array[mg_level],
      //  PArray[mg_level], Cd, Cl);

      PArray[mg_level]->FindGradient(0.45, 0.2, 0.205, dP1);
      PArray[mg_level]->FindGradient(0.55, 0.2, 0.205, dP2);

      OutPut( "Cdrag = " << setprecision(16) <<Cd );
      OutPut( " C_lift = " << setprecision(16) << Cl);
      OutPut( " deltaP = " << setprecision(16) << dP1[0] - dP2[0] << endl);
      OutPut( setprecision(7) << endl);

      if(TDatabase::ParamDB->WRITE_GRAPE)
      {
        os.seekp(std::ios::beg);
        os << GrapeBaseName << i << ".dat" << ends;
        Output->WriteGrape(os.str().c_str());
      }
      if(TDatabase::ParamDB->WRITE_GMV)
      {
        os.seekp(std::ios::beg);
        os << GmvBaseName << i << ".gmv" << ends;
        Output->WriteGMV(os.str().c_str());
	}
      #endif
*/
    }                            // endfor

    if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
      IntoL20FEFunction3D(sol+3*N_U, N_P,PSpaces[mg_level]);
    psi = new double[N_V];
    //StreamFunction(velocity_space, sol, sol+N_U, streamfunction_space, psi);
    //StreamFct = new TFEFunction3D(streamfunction_space, PsiString, PsiString, psi, N_V);
    // DivU(u1, u2, StreamFct);
    // Output->AddFEFunction(StreamFct);

    if(TDatabase::ParamDB->WRITE_GRAPE)
    {
      OutPut("write u " << Ddot(3*N_U,sol,sol)<< endl);
      os.seekp(std::ios::beg);
      os << GrapeBaseName << i << ".dat" << ends;
      Output->WriteGrape(os.str().c_str());
    }
    if(TDatabase::ParamDB->WRITE_GMV)
    {
      os.seekp(std::ios::beg);
      os << GmvBaseName << i << ".gmv" << ends;
      Output->WriteGMV(os.str().c_str());
    }

    if(TDatabase::ParamDB->WRITE_VTK)
    {
      os.seekp(std::ios::beg);
      os << VTKBaseName << i << ".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
    }

    if (TDatabase::ParamDB->SAVE_DATA)
    {
	save_sol[0] = sol;
	save_N_Unknowns[0] = N_Unknowns;
	os << SaveDataFileName << i << ".save" << ends;
	SaveData(SaveDataFileName,1,save_sol,save_N_Unknowns);
    }	    

    // measure errors to known solution
    if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
      fesp[0] = USpaces[mg_level];
      fefct[0] = U1Array[mg_level];
      fefct[1] = U2Array[mg_level];
      fefct[2] = U3Array[mg_level];
      aux =  new TAuxParam3D(NSN_FESpacesVelo, NSN_FctVelo,
        NSN_ParamFctVelo,
        NSN_FEValuesVelo,
        fesp, fefct,
        NSFctVelo,
        NSFEFctIndexVelo, NSFEMultiIndexVelo,
        NSN_ParamsVelo, NSBeginParamVelo);

      // errors in first velocity component
      U1Array[mg_level]->GetErrors(ExactU1, 4, NSAllDerivatives, 2,
        L2H1Errors,
        NULL, aux, 1, USpaces+mg_level, errors);
      l2u1[i] = errors[0];
      h1u1[i] = errors[1];

      // errors in second velocity component
      U2Array[mg_level]->GetErrors(ExactU2, 4, NSAllDerivatives, 2,
        L2H1Errors,
        NULL, aux, 1, USpaces+mg_level, errors);
      l2u2[i] = errors[0];
      h1u2[i] = errors[1];

      // errors in third velocity component
      U3Array[mg_level]->GetErrors(ExactU3, 4, NSAllDerivatives, 2,
        L2H1Errors,
        NULL, aux, 1, USpaces+mg_level, errors);
      l2u3[i] = errors[0];
      h1u3[i] = errors[1];

      // errors in pressure
      PArray[mg_level]->GetErrors(ExactP, 4, NSAllDerivatives, 2,
        L2H1Errors,
        NULL, aux, 1, PSpaces+mg_level, errors);
      l2p[i] = errors[0];
      h1p[i] = errors[1];

      // output of errors
      // coarsest level
      OutPut(l2u1[i] << " AAA " << l2u2[i] << " " << l2u3[i] << endl);
      if(i<=FirstSolve)
      {
        p1 = sqrt(l2u1[i]*l2u1[i]+l2u2[i]*l2u2[i]+l2u3[i]*l2u3[i]);
        OutPut("L2(u): " <<  p1 << endl);
        p2 = sqrt(h1u1[i]*h1u1[i]+h1u2[i]*h1u2[i]+h1u3[i]*h1u3[i]);
        OutPut("H1-semi(u): " <<  p2 << endl);
        OutPut("L2(p): " <<  l2p[i] << endl);
        OutPut("H1-semi(p): " <<  h1p[i] << endl);
      }
      // not coarsest level
      else
      {
        errors[0] = sqrt(l2u1[i]*l2u1[i]+l2u2[i]*l2u2[i]+l2u3[i]*l2u3[i]) ;
        errors[1] = sqrt(l2u1[i-1]*l2u1[i-1]+l2u2[i-1]*l2u2[i-1]+l2u3[i-1]*l2u3[i-1]);
        OutPut("L2(u): " <<  errors[0] << " order ");
        OutPut(log(errors[1]/errors[0])/ln2 << endl);

        errors[0] = sqrt(h1u1[i]*h1u1[i]+h1u2[i]*h1u2[i]+h1u3[i]*h1u3[i]) ;
        errors[1] = sqrt(h1u1[i-1]*h1u1[i-1]+h1u2[i-1]*h1u2[i-1]+h1u3[i-1]*h1u3[i-1]);
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
          (ExactNull, ExactNull, ExactNull,
	   3, NSAllDerivatives,
	   2, DivergenceError,
	   NULL, aux, 1, USpaces+mg_level, errors);
      
      OutPut( "divergence error (L1/L2) : " << errors[0]*errors[0] << " " << errors[1] << endl);
      delete aux;
    }                            // endif MEASURE_ERRORS

    #ifdef __BENCH__
    // compute characteristic values (deltaP, Cd, Cl)
    GetCdCl(U1Array[mg_level], U2Array[mg_level], U3Array[mg_level],
	    PArray[mg_level], Cd, Cl);

    PArray[mg_level]->FindGradient(0.45, 0.2, 0.205, dP1);
    PArray[mg_level]->FindGradient(0.55, 0.2, 0.205, dP2);

    OutPut( "C_drag = " << setprecision(16) <<Cd );
    OutPut( " C_lift = " << setprecision(16) << Cl);
    OutPut( " deltaP = " << setprecision(16) << dP1[0] - dP2[0] << endl);
    OutPut( setprecision(7) << endl);
    #endif
    #ifdef  __CHANNEL_OBSTACLE__
    // compute characteristic values (deltaP, Cd, Cl)
    GetCdCl(U1Array[mg_level], U2Array[mg_level], U3Array[mg_level],
      PArray[mg_level], Cd, Cl);

    PArray[mg_level]->FindGradient(2.45, 0.2, 0.205, dP1);
    PArray[mg_level]->FindGradient(2.55, 0.2, 0.205, dP2);

    OutPut( "C_drag = " << setprecision(16) <<Cd );
    OutPut( " C_lift = " << setprecision(16) << Cl);
    OutPut( " deltaP = " << setprecision(16) << dP1[0] - dP2[0] << endl);
    OutPut( setprecision(7) << endl);
    #endif

    #ifdef __CHANNELSTEP__
    GetReattachmentLine(U1Array[mg_level], reatt_pt);
    //OutPut( "reattachment: " << reatt_pt<< endl);
    #endif

    // in case of convergence, set solution to zero
    // to prevent interpolation of useless data
    /*if ((slow_conv) || sqrt(residual)>limit)
      {
      memset(sol, 0, N_Unknowns*SizeOfDouble);
      OutPut("solution on level " << i <<
        " set to zero since method did not converge" << endl);
	}*/

    // remove data which will not be used later
    delete oldsol;
    delete defect;
    delete psi;

    if ((mg_type==1)||(TDatabase::ParamDB->SOLVER_TYPE == AMG)||
	(TDatabase::ParamDB->SOLVER_TYPE == DIRECT))
    {
      switch(TDatabase::ParamDB->NSTYPE)
      {
        case 1:
          delete matrixB1;
          delete matrixB2;
          delete matrixB3;
          delete sqmatrixA;
 	  if (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION)
	  {
	      delete sqmatrixL;
	      delete  matrix_tilde_G11;
	      delete  matrix_tilde_G22;
	      delete  matrix_tilde_G33;
	      delete  matrix_G11;
	      delete  matrix_G22;
	      delete  matrix_G33;
	      delete  sqstructureL;
	      delete  structure_tilde_G;
	      delete  structure_G;
	      delete  projection_space;
	  }
         break;
        case 2:
          delete matrixB1;
          delete matrixB2;
          delete matrixB3;
          delete matrixB1T;
          delete matrixB2T;
          delete matrixB3T;
          delete sqmatrixA;
          break;
        case 3:
          delete matrixB1;
          delete matrixB2;
          delete matrixB3;
          delete sqmatrixA11;
          delete sqmatrixA12;
          delete sqmatrixA13;
          delete sqmatrixA21;
          delete sqmatrixA22;
          delete sqmatrixA23;
          delete sqmatrixA31;
          delete sqmatrixA32;
          delete sqmatrixA33;
          break;
        case 4:
          delete matrixB1;
          delete matrixB2;
          delete matrixB3;
          delete matrixB1T;
          delete matrixB2T;
          delete matrixB3T;
          delete sqmatrixA11;
          delete sqmatrixA12;
          delete sqmatrixA13;
          delete sqmatrixA21;
          delete sqmatrixA22;
          delete sqmatrixA23;
          delete sqmatrixA31;
          delete sqmatrixA32;
          delete sqmatrixA33;
          break;
      }
      delete structureB;
      delete structureBT;
      delete sqstructureA;
      delete rhs_high;
      delete MGLevel;
    }                            // end if (mg_type==1)

    old_sol = sol;
    old_u = u;
    old_p = p;
    old_u_space = velocity_space;
    old_p_space = pressure_space;

    if (TDatabase::ParamDB->SOLVER_TYPE == GMG)
    {
      delete prec;
      delete itmethod;
      if (TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE==5)
      {
        delete itmethod_sol;
        delete itmethod_rhs;
      }
    }

    OutPut("memory after: " << setw(10) << GetMemory() << endl);

  }                              // endfor i

  OutPut("used time: " << GetTime() << endl);
  OutPut("used bytes: " << GetMemory() << endl);
  CloseFiles();

  return 0;
}
