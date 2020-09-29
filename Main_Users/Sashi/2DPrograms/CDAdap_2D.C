// =======================================================================
//
// Purpose:     main program
//
// Author:      Volker John   August 2000
// 
// =======================================================================

#include <Assemble2D.h>
#include <Assemble2D_edge_Oseen.h>
#include <Domain.h>
#include <Database.h>
#include <DiscreteForm2D.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <AuxParam2D.h>
#include <QuadAffin.h>
#include <Solver.h>
#include <DirectSolver.h>
#include <Constants.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>
#include <LocalProjection.h>
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
#include <FEM_TVD_FCT.h>
#include <MultiGrid2D.h>
#include <MGLevel2D.h>
#include <MainRoutines2D.h>
#include <ItMethod.h>

// =======================================================================
// include current example
// =======================================================================
//#include "../Examples/CD_2D/SineLaplace.h"
//#include "../Examples/CD_2D/TwoBoundaryLayers.h"
//#include "../Examples/CD_2D/TwoBoundaryLayersLaplace.h"
// #include "../Examples/CD_2D/Sphere.h"
// #include "../Examples/CD_2D/Constant1.h"
// #include "../Examples/CD_2D/Cylinder.h"
// #include "../Examples/CD_2D/Plane.h"
// #include "../Examples/CD_2D/MF1.h"
// #include "../Examples/CD_2D/MF2.h"
// #include "../Examples/CD_2D/MF3.h"
// #include "../Examples/CD_2D/BspLutz.h"
//#include "../Examples/CD_2D/Hexagon.h"
// #include "../Examples/CD_2D/Bw_Facing_Step_1_3.h"
// #include "../Examples/CD_2D/LaplaceComputerPraktikum.h"
//#include "../Examples/CD_2D/Robin.h"
// #include "../Examples/CD_2D/Robin_indian_stat.h"
//#include "../Examples/CD_2D/TwoInteriorLayers.h"
//#include "../Examples/CD_2D/PetrSkew.h"
//#include "../Examples/CD_2D/BurmanErn2002.h"
//#include "../Examples/CD_2D/BurmanHansbo2004.h"
//#include "../Examples/CD_2D/HMM1986.h"
//#include "../Examples/CD_2D/JohnMaubachTobiska1997.h"
#include "../Examples/CD_2D/Hemker1996.h"
//#include "../Examples/CD_2D/ParabolicLayers.h"
//#include "../Examples/CD_2D/MG_Vorlesung_rhs1.h"
//#include "../Examples/CD_2D/pw_linear_rhs.h"
//#include "../Examples/CD_2D/TestReaction.h"
//#include "../Examples/CD_2D/SinSin.h"
//#include "../Examples/CD_2D/Smooth.h"
//#include "../Examples/CD_2D/Poisson_Const.h"
//#include "../Examples/CD_2D/Burman2004.h"
//#include "../Examples/CD_2D/Leveque.h"
int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *coll=NULL, *coll_coarse;
  TFESpace2D *velocity_space;
  TFESpace2D *old_u_space, *shock_capturing_space, *pw_const_param_space;
  TFESpace2D **USpaces;
  TOutput2D *Output;

  double *rhs, *sol, *oldsol, *defect, *sol_adjoint, *pw_const_param;
  double *pw_const_param_deriv, *pw_const_param_old;
  double *old_sol, *update, *matrix_D_Entries;
  int i,j,k,l,low, m, k0, iter, n, type, loop_count, inner_loop_count, meas_errors, succ_step = 0;
  int N_U, N_Unknowns, N_Cells, N_RHSs, N_SqMat, N_Fesp, solve_adj = 1;
  int N_Active, N_NonActive;
  double *l2, *h1, *sd, *l_inf, *lp, *energy;
  double  eps, val_x[10], val_y[10], val_r[10], help, damp_adjoint;
  double damp_adjoint_start = 1e-6;
  double observed_functional, observed_functional_old = -1;
  char *PRM, *GEO;
  int LEVELS, BASELEVEL;
  int ret, ORDER;
  double errors[6],score[4];
  double t1, t2, res, res2, oldres, solver_time, hmin, hmax, oldres_stepm1;
  double nonlin_min_res;
  int N_LinIter, first_damp,linite, compute_matrix_D;

  std::ostringstream os;
  char *PsBaseName, *ReadGrapeBaseName;
  int current_estimator = TDatabase::ParamDB->ADAPTIVE_REFINEMENT_CRITERION;
  double maximal_local_error,estimated_global_error[5],*eta_K;
  TCD2DErrorEstimator *CDErrorEstimator;

  // TFEFunction2D **AuxFEFunctArray, **AuxFEVectFunctArray;
  TFEFunction2D *u, **UArray, *old_u, *u_adjoint, *pw_const_param_fe, *pw_const_param_deriv_fe;
  TFESpace2D *fesp[2], *ferhs[3];

  TAuxParam2D *aux;

  TSquareStructure2D *sqstructureA;
  TSquareMatrix2D *sqmatrixA, *SQMATRICES[1], *sqmatrixAadjoint;
  TSquareMatrix2D **MatricesA;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  double **RhsArray;

  TMGLevel2D *MGLevel;
  TMultiGrid2D *MG;

  double *RHSs[1];
  int *N_Uarray;
  char UString[] = "u";
  char CdString[] = "Conv-Diff";
  char GalString[] = "Galerkin";
  char SDFEMString[] = "SDFEM";
  char UpwString[] = "Upwind";
  char SoldString[] = "SOLD";
  char Readin[] = "readin.dat";
  char Name[] = "name";
  char Description[] = "description";

  int FirstSolve;
  int N_Parameters=2,n_aux;
  double Parameters[2];
  int mg_level,mg_type;
  int *neum_to_diri, N_neum_to_diri = 0, *neum_to_diri_bdry;
  double *neum_to_diri_param;

  double *sc_params, residual, omega, omega_max = 1.0;
  double omega_min = 0.01, *rhs_edge;
  int N_sc, sc_N_params, max_it, sold_parameter_type;
  double lin_red, res_norm_min;
  TFEFunction2D *sc_params_fe[2], *fefct[3];
  TFEVectFunct2D *sc_params_vect;
  DefectProc *Defect;

  double values[4];

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
    OutPut("NO READIN FOUND !!!"<<endl);
    exit(-1);
  }

  ManipulateFct2D *manipulate;
  if (TDatabase::ParamDB->SDFEM_NORM_B==0)
    manipulate = linfb;
  else
    manipulate = ave_l2b_quad_points;

  TDiscreteForm2D *DiscreteForm, *DiscreteForms[10];
  TDiscreteForm2D *DiscreteFormGalerkin = new TDiscreteForm2D
    (CdString, GalString, N_Terms, Derivatives, SpacesNumbers,
    N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace,
    BilinearAssemble, BilinearCoeffs, NULL);

  TDiscreteForm2D *DiscreteFormSDFEM = new TDiscreteForm2D
    (CdString, SDFEMString, N_Terms_SD, Derivatives_SD, SpacesNumbers_SD,
    N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace,
    BilinearAssemble_SD, BilinearCoeffs,  manipulate);

  TDiscreteForm2D *DiscreteFormUpwind = new TDiscreteForm2D
    (CdString, UpwString , N_Terms, Derivatives, SpacesNumbers,
    N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace,
    BilinearAssemble_UPW2, BilinearCoeffs, NULL);

  TDiscreteForm2D *DiscreteFormSOLD = new TDiscreteForm2D
    (CdString, SoldString, N_Terms_SD, Derivatives_SD, SpacesNumbers_SD,
    N_Matrices, 0 , RowSpace, ColumnSpace, NULL,
    BilinearAssemble_SOLD, BilinearCoeffs,  manipulate);

  TDiscreteForm2D *DiscreteFormSOLD_Orthogonal = new TDiscreteForm2D
    (CdString, SoldString, N_Terms_SD, Derivatives_SD, SpacesNumbers_SD,
    N_Matrices, 0 , RowSpace, ColumnSpace, NULL,
    BilinearAssemble_SOLD_Orthogonal, BilinearCoeffs,  manipulate);

  TDiscreteForm2D *DiscreteFormRhsLP96 = new TDiscreteForm2D
    (CdString, SDFEMString, N_Terms, Derivatives, SpacesNumbers,
    0, 1 , NULL, NULL,  RhsSpace,
    RhsAssemble_LP96, BilinearCoeffs,  manipulate);

  TDiscreteForm2D *DiscreteFormMH_Kno06 = new TDiscreteForm2D
    (CdString, SoldString, N_Terms_SD, Derivatives_SD, SpacesNumbers_SD,
    N_Matrices, 0 , RowSpace, ColumnSpace, NULL,
    BilinearAssemble_MH_Kno06, BilinearCoeffs,  manipulate);

  TDiscreteForm2D *DiscreteFormRhsAdjointEnergyErrorEstimate = new TDiscreteForm2D
    (CdString, SDFEMString, N_Terms_SD,
    Derivatives_SD, SpacesNumbers_SD,
    0, 1 , NULL, NULL,  RhsSpace,
    RhsAssemble_RhsAdjointEnergyEstimate, BilinearCoeffs,  manipulate);

  TDiscreteForm2D *DiscreteFormRhsAdjointL2Error = new TDiscreteForm2D
    (CdString, SDFEMString, N_Terms_rhs,
    Derivatives_rhs, SpacesNumbers_rhs,
    0, 1 , NULL, NULL,  RhsSpace,
    RhsAssemble_RhsAdjointL2Error, BilinearCoeffs,  manipulate);

  TDiscreteForm2D *DiscreteFormRhsAdjointH1Error = new TDiscreteForm2D
    (CdString, SDFEMString, N_Terms_SD,
    Derivatives_SD, SpacesNumbers_SD,
    0, 1 , NULL, NULL,  RhsSpace,
    RhsAssemble_RhsAdjointH1Error, BilinearCoeffs,  manipulate);

  DiscreteForms[0] = DiscreteFormGalerkin;
  DiscreteForms[1] = DiscreteFormSDFEM;
  DiscreteForms[2] = DiscreteFormUpwind;
  DiscreteForms[3] = DiscreteFormSOLD;
  DiscreteForms[4] = DiscreteFormSOLD_Orthogonal;
  DiscreteForms[5] = DiscreteFormRhsLP96;
  DiscreteForms[6] = DiscreteFormMH_Kno06;
  DiscreteForms[7] = DiscreteFormRhsAdjointEnergyErrorEstimate;
  DiscreteForms[8] = DiscreteFormRhsAdjointL2Error;
  DiscreteForms[9] = DiscreteFormRhsAdjointH1Error;

  BoundCondFunct2D *BoundaryConditions[3] =
  {
    BoundCondition, BoundCondition_FEM_FCT,
    BoundConditionAdjoint
  };
  BoundValueFunct2D *BoundaryValues[3] =
  {
    BoundValue, BoundValue_FEM_FCT,
    BoundValueAdjoint
  };
  CoeffFct2D *Coefficients[1];
  CheckWrongNeumannNodesFunct2D *CheckWrongNeumannNodesFct[1] ={CheckWrongNeumannNodes};

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  // sets parameters of the data base
  SetParametersCDAdapt2D();
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  // write data base
  Database->WriteParamDB(argv[0]);
  ExampleFile();

  //======================================================================
  // copy read parameters into local variables
  //======================================================================

  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;

  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  ReadGrapeBaseName = TDatabase::ParamDB->READGRAPEBASENAME;

  sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;
  Coefficients[0] = BilinearCoeffs;

  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SCALAR;
  if (TDatabase::ParamDB->SOLVER_TYPE==AMG_SOLVE||
    TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
    mg_type = 0;
  if (mg_type)
    mg_level = 0;
  else
    mg_level = -1;
  LEVELS = TDatabase::ParamDB->LEVELS;
  BASELEVEL = TDatabase::ParamDB->UNIFORM_STEPS;
  l2 = new double[LEVELS+1];
  h1 = new double[LEVELS+1];
  energy = new double[LEVELS+1];
  sd = new double[LEVELS+1];
  l_inf = new double[LEVELS+1];
  lp = new double[LEVELS+1];

  UArray = new TFEFunction2D*[LEVELS+1];
  RhsArray = new double* [LEVELS+1];
  N_Uarray = new int[LEVELS+1];

  USpaces = new TFESpace2D*[LEVELS+1];
  MatricesA = new TSquareMatrix2D*[LEVELS+1];

  //======================================================================
  // read boundary parameterization and initialize coarse grid
  //======================================================================
  Domain->Init(PRM, GEO);

  // refine up to user defined coarsest level

  for(i=0;i<TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR;i++)
    Domain->RegRefineAll();

  // initialize solver parameters
  Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
  Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;

  if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    i=1;
    MG = new TMultiGrid2D(i, N_Parameters, Parameters);
  }

  FirstSolve = TDatabase::ParamDB->SC_FIRST_SOLUTION_LEVEL_SCALAR;
  meas_errors = TDatabase::ParamDB->MEASURE_ERRORS;
  TDatabase::ParamDB->MEASURE_ERRORS = 0;
  //======================================================================
  // loop over all levels
  //======================================================================

  for(i=0;i<LEVELS;i++)
  {
    mg_level++;
    OutPut("*******************************************************" << endl);
    OutPut("******           GEOMETRY  LEVEL ");
    OutPut(i << "              *******" << endl);
    OutPut("******           MULTIGRID LEVEL ");
    OutPut(mg_level << "              *******" << endl);
    OutPut("*******************************************************" << endl);
    solver_time =  GetTime();
    N_LinIter = 0;
    OutPut("memory before: " << setw(10) << GetMemory() << endl);

    // refine grid if level is greater than 0
    if ((i)&&(i<LEVELS))
    {
      if ((i<= BASELEVEL)||(! TDatabase::ParamDB->ESTIMATE_ERRORS))
      {
        Domain->RegRefineAll();
      }
      else if (i < LEVELS)
      {
        if (TDatabase::ParamDB->GRID_TYPE)
          Domain->RefineByErrorEstimator(coll,eta_K,maximal_local_error,
            estimated_global_error[current_estimator],
            TRUE);
        else
          Domain->RefineByErrorEstimator(coll,eta_K,maximal_local_error,
            estimated_global_error[current_estimator],
            FALSE);
        delete eta_K;
      }
    }
    // loop over all levels, build new fe spaces for changed
    // adaptive discretization
    if ( TDatabase::ParamDB->ESTIMATE_ERRORS)
      k0 = 0;
    else
      k0 = mg_level-mg_type;
    for(k=k0;k<=mg_level;k++)
    {
      if (coll!=NULL)
        coll_coarse = coll;
      coll=Domain->GetCollection(It_LE, k+TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR);
      if (k==mg_level)
      {
        N_Cells = coll->GetN_Cells();
        OutPut( "number of cells: " << N_Cells << endl);
      }
      // for adaptive methods, all levels are created new
      // but picture only on the finest level
      if((k==mg_level)&&(TDatabase::ParamDB->WRITE_PS))
      {
        // write grid into an Postscript file
        os.seekp(std::ios::beg);
        os << PsBaseName << i << ".ps" << ends;
        Domain->PS(os.str().c_str(),It_Finest,0);
      }

      // fe space without boundary conditions
      if ((TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD)
        || (TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM))
      {
        BoundaryConditions[0] = BoundCondition_FEM_FCT;
        BoundaryValues[0] = BoundValue_FEM_FCT;
      }
      // get spaces for low order disc on finest geo grid
      if ((mg_type==1)&&(k<mg_level))
      {
        velocity_space =  new TFESpace2D(coll, Name,
          Description,
          BoundaryConditions[0], ORDER, NULL);
      }
      // get spaces of high order disc on grid
      else
        velocity_space =  new TFESpace2D(coll, Name,
          Description,
          BoundaryConditions[0], ORDER, NULL);

      if ((TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD)
        || (TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM))
      {
        BoundaryConditions[0] = BoundCondition;
        BoundaryValues[0] = BoundValue;
      }

      // build fespace hierarchy
      // set values and pointers for fe space
      USpaces[k] = velocity_space;
      N_U = velocity_space->GetN_DegreesOfFreedom();
      N_Uarray[k] = N_U;
      N_Active = velocity_space->GetActiveBound();
      N_NonActive = N_U - N_Active;

      // SOLD schemes
      if (TDatabase::ParamDB->SOLD_TYPE)
      {
        // allocate piecewise constant finite element space
        shock_capturing_space = new TFESpace2D(coll, Name, Description,
          BoundCondition, 0, NULL);
        sc_N_params = 2;
        N_sc = shock_capturing_space->GetN_DegreesOfFreedom();
        if (sc_N_params)
        {
          sc_params = new double[sc_N_params* N_sc];
          sc_params_vect = new TFEVectFunct2D(shock_capturing_space, UString, UString, sc_params,
            N_sc, sc_N_params);
          for (m=0;m<sc_N_params;m++)
            sc_params_fe[m] = sc_params_vect->GetComponent(m);
        }
        if ((sold_parameter_type >= BH04)&& (sold_parameter_type <= BE05_2))
          rhs_edge = new double[N_U];
        if (sold_parameter_type == MH_Kno06)
          rhs_edge = new double[N_U];
      }

      // build matrices
      sqstructureA = new TSquareStructure2D(velocity_space);
      sqstructureA->Sort();
      // allocate matrices
      sqmatrixA = new TSquareMatrix2D(sqstructureA);
      MatricesA[k] = sqmatrixA;
      N_Unknowns = N_U;

      if (TDatabase::ParamDB->SOLD_TYPE)
      {
        if (sold_parameter_type ==  FEM_TVD)
        {
          matrix_D_Entries = new double[sqmatrixA->GetN_Entries()];
          memset(matrix_D_Entries, 0, sqmatrixA->GetN_Entries()*SizeOfDouble);
          rhs_edge = new double[N_U];
          memset(rhs_edge, 0, N_U * SizeOfDouble);
        }
      }

      if (((k>=FirstSolve)||(mg_type==0))&&(k==mg_level))
        OutPut("dof all      : "<< setw(10) << N_Unknowns  << endl);
      if ((mg_type==1)&&(k==mg_level-1))
      {
        OutPut("dof low order disc     : "<<  setw(10) << N_Unknowns  << endl);
      }

      coll->GetHminHmax(&hmin,&hmax);
      OutPut("h_min : " << hmin << " h_max : " << hmax << endl);

      // initialize arrays
      rhs = new double[N_Unknowns];
      memset(rhs, 0, N_Unknowns*SizeOfDouble);
      RhsArray[k] = rhs;
      sol = new double[N_Unknowns];
      oldsol = new double[N_Unknowns];
      update = new double[N_Unknowns];
      memset(sol, 0, N_Unknowns*SizeOfDouble);
      memset(oldsol, 0, N_Unknowns*SizeOfDouble);
      memset(update, 0, N_Unknowns*SizeOfDouble);

      // build multigrid level(s)
      // ( A B' )
      // ( B 0  )
      switch(TDatabase::ParamDB->SOLVER_TYPE)
      {
        case AMG_SOLVE:
        case DIRECT:
          low = mg_level;
          break;

        case GMG:
          // coarsest grid number
          low = 0;
          // determine number of auxiliary arrays
          if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR)
            || (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR))
            n_aux= 4;
          else
            n_aux= 2;
          // build fe multigrid levels
          MGLevel = new TMGLevel2D(k, sqmatrixA,
            rhs, sol,
            n_aux, NULL);
          if ((k==mg_level)||((k==0)&&(mg_level==1)&&(mg_type==1)))
            MG->AddLevel(MGLevel);
          else
            MG->ReplaceLevel(k,MGLevel);
          break;
      }

      // build new fe functions
      u = new TFEFunction2D(velocity_space, UString, UString, sol, N_U);
      UArray[k] = u;
      if ((TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)&&((k==mg_level)))
        PrepareAdjointProblem2D(coll,
          u_adjoint, sqmatrixAadjoint,
          pw_const_param_space, pw_const_param_fe,
          pw_const_param_deriv_fe,
          velocity_space, sqstructureA,
          BoundCondition, BilinearCoeffs,
          sol_adjoint, rhs_edge, pw_const_param,
          pw_const_param_deriv, pw_const_param_old,
          N_U, N_Cells);

      // prepare data output
      if (k==mg_level)
      {
        if (!TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)
        {
          Output = new TOutput2D(1, 1, 1, 1,Domain);
          Output->AddFEFunction(u);
        }
        else
        {
          Output = new TOutput2D(1, 4, 1, 1,Domain);
          Output->AddFEFunction(u);
          Output->AddFEFunction(u_adjoint);
          Output->AddFEFunction(pw_const_param_fe);
          Output->AddFEFunction(pw_const_param_deriv_fe);
        }
      }

      // prolongation, to get a good starting iterate
      // HAS TO BE TESTED FOR ADAPTIVE GRIDS
      if ((k==mg_level) && k>FirstSolve+mg_type)
      {
        // for assessment of LCD vs. GMRES, can be cancelled after the assessment
        //if (!((TDatabase::ParamDB->SOLVER_TYPE==AMG_SOLVE))&&
        //((TDatabase::ParamDB->SC_SOLVER_SCALAR==14)||(TDatabase::ParamDB->SC_SOLVER_SCALAR==19)))
        //
        OutPut("prolongate solution from next coarser level"<<endl);
        Prolongate(old_u_space, USpaces[mg_level],
          old_u->GetValues(), UArray[mg_level]->GetValues(),
          oldsol);
        //OutPut("Starting iterate is zero !"<<endl);
        // copy current solution for assembling the nonlinear terms
        // memcpy(oldsol,sol, N_Unknowns*SizeOfDouble);
        if (mg_type==1)
        {
          delete old_sol;
          delete old_u;
          delete old_u_space;
        }
      }                                           // end of prolongate
    }

    // restrict solution to all grids
    if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
      MG->RestrictToAllGrids();

    // if no solution on this grid, continue
    if(FirstSolve>i)
      continue;
    // TEST
    /*TBaseCell *cell;
    coll->MarkBoundaryVertices();
    for (j=0;j<N_Cells;j++)
    {
    cell = coll->GetCell(j);
    pw_const_param[j] = ((TGridCell *)cell)->GetN_BoundaryEdges();
    }
    OutputData2D(os, Output,4711);
    exit(1);*/

    loop_count = 0;
    // **************************************************************************
    // loop over the update of stabilization parameters
    inner_loop_count=0;
    observed_functional_old = -1;
    while(loop_count<TDatabase::ParamDB->SC_NONLIN_ITE_ADJOINT)
    {
      // assembling of the primal problem
      Assemble_CD_2D(coll, USpaces, UArray,
        RhsArray, MatricesA,
        Coefficients[0],
        BoundaryConditions,
        BoundaryValues,
        DiscreteForms,
        sqmatrixAadjoint,
        CheckWrongNeumannNodesFct,
        N_Uarray,
        low, mg_level, mg_type, i,
        N_neum_to_diri, neum_to_diri,
        neum_to_diri_bdry,
        neum_to_diri_param);
      // set rhs for Dirichlet nodes
      memcpy(sol+N_Active, rhs+N_Active, N_NonActive*SizeOfDouble);
      // solve linear system
      SQMATRICES[0] = MatricesA[mg_level];
      Solver(sqmatrices, NULL,
        RhsArray[mg_level], sol,
        MatVect_Scalar, Defect_Scalar,
        MG,  N_Unknowns, 0);

      if (TDatabase::ParamDB->SC_VERBOSE>1)
        OutPut("time for linear problem: " << GetTime() - solver_time << endl);

      // the solution of the adjoint problem is required
      if (TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)
      {
        // compute error estimate
        switch(TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)
        {
          case 1:
            ComputeErrorEstimate(coll, u,
              Coefficients[0], BoundaryConditions,
              BoundaryValues, eta_K,
              &maximal_local_error,
              estimated_global_error, l2[i], h1[i], N_Unknowns);
            observed_functional = estimated_global_error[3];
            break;
          case 2:
            ComputeErrorEstimate(coll, u,
              Coefficients[0], BoundaryConditions,
              BoundaryValues, eta_K,
              &maximal_local_error,
              estimated_global_error, l2[i], h1[i], N_Unknowns);
            observed_functional = estimated_global_error[4];
            break;
          case 4:
            // mark vertices on Dirichlet boundaries
            coll->MarkBoundaryVertices();
            ComputeErrorEstimate(coll, u,
              Coefficients[0], BoundaryConditions,
              BoundaryValues, eta_K,
              &maximal_local_error,
              estimated_global_error, l2[i], h1[i], N_Unknowns);
            observed_functional = estimated_global_error[4];
            break;
          case 100:
            fesp[0] = USpaces[i];
            aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
            u->GetErrors(Exact, 3, AllDerivatives, 4, SDFEMErrors,
              BilinearCoeffs, aux, 1, fesp, errors);
            delete aux;
            observed_functional = errors[0];
            break;
          case 101:
            fesp[0] = USpaces[i];
            aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
            u->GetErrors(Exact, 3, AllDerivatives, 4, SDFEMErrors,
              BilinearCoeffs, aux, 1, fesp, errors);
            delete aux;
            observed_functional = errors[1];
            break;
        }
        inner_loop_count++;
        if ((loop_count == 0) && (solve_adj))
        {
          OutPut("initial observed functional " << observed_functional << endl);
          // write data for pictures if wished
          //OutputData2D(os, Output,loop_count);
        }

        // adjoint problem has to be assembled and solved
        if (solve_adj)
        {
          // copy current solution
          Daxpy(N_Unknowns,-1,sol,oldsol);
          // copy current solution
          memcpy(oldsol,sol, N_Unknowns*SizeOfDouble);
          // assign DiscreteForm
          switch(TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)
          {
            case 1:
            case 2:
            case 4:
              DiscreteForm = DiscreteFormRhsAdjointEnergyErrorEstimate;
              break;
            case 100: DiscreteForm = DiscreteFormRhsAdjointL2Error;
            break;
            case 101: DiscreteForm = DiscreteFormRhsAdjointH1Error;
            break;
            default: OutPut("SOLVE_ADJOINT_PROBLEM no defined !!!" << endl);
            exit(4711);
          }

          // assemble and solve adjoint problem
          // computes derivatives of error estimator wrt stabilization parameters
          SolveAdjointProblem2D(coll,
            DiscreteForm,
            USpaces[mg_level],
            UArray[mg_level],
            u_adjoint,
            sqmatrixAadjoint,
            BilinearCoeffs,
            BoundaryConditions,
            BoundaryValues,
            rhs, rhs_edge,
            sol_adjoint,
            pw_const_param_deriv,
            N_Unknowns, N_Active, N_neum_to_diri,
            neum_to_diri, neum_to_diri_bdry,
            neum_to_diri_param);
          // damping factor for update of stabilization parameters
          damp_adjoint = damp_adjoint_start;
          // save parameter array of previous step
          memcpy(pw_const_param_old, pw_const_param, N_Cells*SizeOfDouble);
          //  update stabilization parameters
          Daxpy(N_Cells,-damp_adjoint,pw_const_param_deriv,pw_const_param);
	  // restrict parameters
	  RestrictParametersToAdmissibleValues(coll,pw_const_param_space,Coefficients[0],pw_const_param);
	  // solution of adjoint problem done for this iteration
          solve_adj = 0;
          // no successful step
          succ_step = 0;
          // current value of functional (for standard parameter choice)
          observed_functional_old = observed_functional;
        }
        else
        {
          // improvement of error estimator, increase damp_adjoint
          if ((observed_functional_old<0)||(observed_functional<observed_functional_old))
          {
            // old error estimate
            observed_functional_old = observed_functional;
            // new damping factor
            damp_adjoint *= 2;
            OutPut("new damp " << damp_adjoint << endl);
            // parameter array of previous step
            memcpy(pw_const_param, pw_const_param_old, N_Cells*SizeOfDouble);
            // new stabilization parameters
            Daxpy(N_Cells,-damp_adjoint,pw_const_param_deriv,pw_const_param);
	    // restrict parameters
	    RestrictParametersToAdmissibleValues(coll,pw_const_param_space,Coefficients[0],pw_const_param);
            // there was a damping parameter with an improvement of the error estimate
            succ_step = 1;
            // now: solve problem with new stabilization parameter
            continue;
          }
          else
          {
            // the error estimate became worse with the new stabilization parameters
            // there are two situations:
            //   1. there was a previous damping parameter with an improvement of the
            //      error estimate, then take this damping parameter
            //   2. there was not such a damping parameter, then decrease the current
            //      damping factor
            // case 1
            if ((succ_step)||( inner_loop_count==100)||(damp_adjoint<1e-12))
            {
              // there were already successful step length
              // take last successful step length
              // parameter array of previous step
              memcpy(pw_const_param, pw_const_param_old, N_Cells*SizeOfDouble);
              // last successful step length
              if (succ_step)
                damp_adjoint /= 2.0;
              // last successful stabilization parameters
              Daxpy(N_Cells,-damp_adjoint,pw_const_param_deriv,pw_const_param);
              // restrict parameters
	      RestrictParametersToAdmissibleValues(coll,pw_const_param_space,Coefficients[0],pw_const_param);

              if ((loop_count == TDatabase::ParamDB->SC_NONLIN_ITE_ADJOINT-1) || (damp_adjoint<1e-12))
              {
                OutputData2D(os, Output,i);
                OutPut(" ACCEPTED: ");
              }
              else
              {
                OutPut(" accepted: ");
              }
              OutPut(loop_count << " functional " << observed_functional_old <<
                " damp " << damp_adjoint << " para "
                << sqrt(Ddot(N_Cells,pw_const_param,pw_const_param)) <<
                " sol "  << sqrt(Ddot(N_Unknowns,sol,sol)) <<
                " adj " << sqrt(Ddot(N_Unknowns,sol_adjoint,sol_adjoint)) <<
                " deriv " << sqrt(Ddot(N_Cells,pw_const_param_deriv,pw_const_param_deriv)) << endl);
              // write data for pictures if wished
              //OutputData2D(os, Output,loop_count+1);
              // set parameter for necessity of solving the adjoint problem in the next loop
              solve_adj = 1;
              // increase loop counter
              loop_count++;
              inner_loop_count=0;
              if (damp_adjoint<1e-12)
                break;
            }
            else
            {
              // no successful step length
              // reduce damping
              // recover solution of last successful step
              // memcpy(sol,oldsol, N_Unknowns*SizeOfDouble);
              // parameter array of previous step
              memcpy(pw_const_param, pw_const_param_old, N_Cells*SizeOfDouble);
              // recover stabilization parameters of last succesful step
              damp_adjoint /= 2;
              OutPut("red damp " << damp_adjoint << endl);
              // new stabilization parameters
              Daxpy(N_Cells,-damp_adjoint,pw_const_param_deriv,pw_const_param);
              // cut negative damping parameters
              // restrict parameters
	      RestrictParametersToAdmissibleValues(coll,pw_const_param_space,Coefficients[0],pw_const_param);
            }
          }
        }
      }
      else
        // no adjoint problem to solve
        break;
    }                                             // end loop over update of parameters
    // **************************************************************************

    // ************************************************************* //
    // end of iteration for adjoint problem
    // ************************************************************* //
    m = 0;
    linite = 0;
    compute_matrix_D = 1;
    nonlin_min_res = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR;
#ifdef   __HEMKER1996__
    if ((i==6) && (TDatabase::ParamDB->RE_NR > 5e4))
    {
      nonlin_min_res = 1e-8;
      OutPut("stopping criteria changed to " << nonlin_min_res * sqrt(N_Unknowns) << endl);
    }
    nonlin_min_res *= sqrt(N_Unknowns);
#endif

    //memset(sol, 0, N_Unknowns*SizeOfDouble);

    // ************************************************************* //
    // SOLD
    // ************************************************************* //
    if ((TDatabase::ParamDB->DISCTYPE==2)&&(TDatabase::ParamDB->SOLD_TYPE))
    {
      lin_red =  TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR;
      TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR =  TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR_SOLD;
      max_it = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
      TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR_SOLD;
      res_norm_min = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
      TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR_SOLD;

      // fixed damping parameter
      if (TDatabase::ParamDB->P5==1234)
      {
        omega = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SCALAR-0.05*mg_level;
        omega_max = omega;
      }
      else
        //omega = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SCALAR;
        omega = omega_max;                        // from the next coarser level
      //omega_max = 1.0;
      oldres_stepm1 = 1.0;
      Defect = Defect_Scalar;
      fesp[0] = USpaces[mg_level];
      fesp[1] = shock_capturing_space;
      SQMATRICES[0] = MatricesA[mg_level];
      memcpy(oldsol,sol, N_Unknowns*SizeOfDouble);
      // defect
      defect = new double[N_Unknowns];
      memset(defect,0,N_Unknowns*SizeOfDouble);
      if (((sold_parameter_type >= BH04)&&(sold_parameter_type <= BE05_2))
        || (sold_parameter_type == MH_Kno06))
        memset(rhs_edge, 0,  N_Unknowns*SizeOfDouble);
      // nonlinear loop
      for (m=0;m< TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR;m++)
      {
        first_damp = 1;
        // loop on the damping parameter
        while (1)
        {
          // compute new iterate with current damping parameter
          Dsum(N_Unknowns,1,omega,oldsol,update,sol);

          if (TDatabase::ParamDB->SOLD_TYPE)
          {
            fefct[0]=UArray[mg_level];
            fefct[1]=sc_params_fe[0];
            fefct[2]=sc_params_fe[1];
            // compute parameters for assembling, e.g. norm of residual in mesh cells
            aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
            u->GetMeshCellParams(ExactNull, 5, Derivatives_SD, 2, Parameters_DC_CD,
              BilinearCoeffs, aux, 1, fesp, errors,
              sc_params);
            for (k=0;k<2*N_sc;k++)
              sc_params[k] = sqrt(sc_params[k]);
            delete aux;

            // set aux object
            aux =  new TAuxParam2D(SOLD_N_FESpaces, SOLD_N_Fct, SOLD_N_ParamFct,
              SOLD_N_FEValues, fesp, fefct,
              SOLD_Fct,
              SOLD_FEFctIndex, SOLD_FEMultiIndex,
              SOLD_N_Params, SOLD_BeginParam);

            // set discrete form
            if (TDatabase::ParamDB->SOLD_TYPE==1)
              DiscreteForm = DiscreteFormSOLD;
            else
              DiscreteForm = DiscreteFormSOLD_Orthogonal;
            if (sold_parameter_type == MH_Kno06)
              DiscreteForm = DiscreteFormMH_Kno06;
            N_RHSs = 0;
            RHSs[0] = NULL;
            ferhs[0] = NULL;
            N_SqMat = 1;
            // defintions for LP96
            if (sold_parameter_type!= LP96)
              SQMATRICES[0]->Reset();
            else
            {
              N_RHSs = 1;
              RHSs[0] = rhs;
              memset(rhs, 0, N_Unknowns*SizeOfDouble);
              ferhs[0] = fesp[0];
              N_SqMat = 0;
              DiscreteForm = DiscreteFormRhsLP96;
            }
          }

          // assemble
          Assemble2D(2, fesp,
            N_SqMat, SQMATRICES,
            0, NULL,
            N_RHSs, RHSs, ferhs,
            DiscreteForm,
            BoundaryConditions,
            BoundaryValues,
            aux);

          delete aux;

          if(TDatabase::ParamDB->DISCTYPE == LOCAL_PROJECTION)
          {
            if(TDatabase::ParamDB->LP_FULL_GRADIENT)
              ;                                   //UltraLocalProjection(MatricesA[mg_level], FALSE, Coefficients[0]);

            if(TDatabase::ParamDB->LP_STREAMLINE)
            {
              OutPut("local projection stabilisation in streamline direction ");
              OutPut("is currently not available." << endl);
            }
          }

          // prepare rhs for edge stabilizations
          if ((sold_parameter_type >= BH04)&&(sold_parameter_type <=  BE05_2))
          {
            // restore rhs
            Daxpy(N_Active, 1, rhs_edge, rhs);

            memset(rhs_edge, 0,  N_Unknowns*SizeOfDouble);
            EdgeStabilization(USpaces[mg_level], UArray[mg_level],
              Coefficients[0], rhs_edge, 0, NULL, NULL);
            // subtract rhs_edge from rhs
            Daxpy(N_Active, -1, rhs_edge, rhs);
          }

          if (sold_parameter_type == MH_Kno06)
          {
            // restore rhs
            Daxpy(N_Active, -1, rhs_edge, rhs);

            memset(rhs_edge, 0,  N_Unknowns*SizeOfDouble);

            MizukamiHughes(sqmatrixA, rhs_edge, USpaces[mg_level],
              UArray[mg_level], Coefficients[0],
              BoundCondition);

            // add rhs_edge to rhs
            Daxpy(N_Active, 1, rhs_edge, rhs);
          }

          if (sold_parameter_type ==  FEM_TVD)
          {
            if (!first_damp)
              memcpy(rhs,rhs_edge,N_U*SizeOfDouble);
            OutPut("RHS_EDGE " << sqrt(Ddot(N_Unknowns,rhs_edge,rhs_edge)) <<
              " SOL " << sqrt(Ddot(N_Unknowns,sol,sol)) << endl);
            memcpy(rhs_edge,rhs,N_U*SizeOfDouble);
            FEM_TVD_ForConvDiff(sqmatrixA, N_U, N_Active, matrix_D_Entries,
              sol, rhs, N_neum_to_diri, neum_to_diri,
              compute_matrix_D);
            compute_matrix_D = 0;
          }

          if (N_neum_to_diri)
            SetDirichletNodesFromNeumannNodes(SQMATRICES, rhs, UArray[mg_level]->GetValues(),
              N_neum_to_diri, neum_to_diri,
              neum_to_diri_bdry, neum_to_diri_param,
              BoundValue);

          // set rhs for Dirichlet nodes
          memcpy(sol+N_Active, rhs+N_Active, N_NonActive*SizeOfDouble);

          // compute defect of the proposed new iterate
          memset(defect,0,N_Unknowns*SizeOfDouble);
          Defect(sqmatrices,NULL,sol,rhs,defect);
          residual =  Ddot(N_Unknowns,defect,defect);
          // first loop
          if (m==0)
          {
            oldres = sqrt(residual);
            break;
          }
          // fixed damping parameter
          if (TDatabase::ParamDB->P5==1234)
          {
            OutPut("fixed damping " << omega << " " << omega_max << endl);
            oldres = sqrt(residual);
            break;
          }
          if ((sqrt(residual) < oldres)||(omega<=omega_min*1.001))
          {
            OutPut("damping " << omega << " " << omega_max << endl);
            // first omega is accepted, then increase omega
            if ((first_damp)&&(sqrt(residual) < oldres))
            {
              omega *= 1.1;
              omega_max *= 1.001;
              if (omega_max > 1)
                omega_max = 1;
              if (omega > omega_max)
                omega = omega_max;
            }
            oldres = sqrt(residual);
            break;
          }
          else
          {
            if (first_damp)
            {
              omega_max  = 0.9*omega_max;
              if (omega_max < omega_min)
                omega_max = omega_min;
            }
            // new residual too large, decrease omega
            omega /=2;
            if (omega<omega_min)
              omega = omega_min;
            first_damp = 0;
          }
        }                                         // end loop of the damping parameter

        memcpy(oldsol, sol, N_Unknowns*SizeOfDouble);

        OutPut("nonlinear iteration step " << m << " residual " <<
          oldres << " " << oldres/oldres_stepm1 << " res/rhs " << sqrt(residual)/sqrt(Ddot(N_Unknowns,rhs,rhs)) << endl);
        oldres_stepm1 = oldres;
        //OutPut("MEMORY: " << setw(10) << GetMemory()/(1048576.0));
        //OutPut(" MB" << endl);

        if (sqrt(residual)<=nonlin_min_res)
          //if (sqrt(residual)/sqrt(Ddot(N_Unknowns,rhs,rhs))<=TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR)
        {
          solver_time =  GetTime() - solver_time;
          OutPut("NONLINEAR ITERATIONS: " << m << " (" << linite << ")  residual " <<
            sqrt(residual)  << " res/rhs " << sqrt(residual)/sqrt(Ddot(N_Unknowns,rhs,rhs)) << " time " << solver_time << endl);
          break;
        }

        // solve system
        Solver(sqmatrices, NULL,
          RhsArray[mg_level], sol,
          MatVect_Scalar, Defect_Scalar,
          MG, N_Unknowns, 0);

        // the update is computed
        Dsum(N_Unknowns, 1, -1, sol, oldsol, update);
        //for (k=0;k<N_Unknowns;k++)
        //  sol[k] = omega *sol[k] + (1-omega)*oldsol[k];

        if (sold_parameter_type ==  FEM_TVD)
        {
          memcpy(rhs,rhs_edge,N_U*SizeOfDouble);
        }

        if (m==TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR-1)
        {
          m+=1;
          OutPut("NONLINEAR ITERATIONS: " << m << " (" << linite << ")  residual " <<
            sqrt(residual) << endl);
          break;
        }
      }

      TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR = lin_red;
      TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR = max_it;
      TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR = res_norm_min;
    }
    // ************************************************************* //
    // end of iteration for nonlinear methods
    // ************************************************************* //

    OutPut("MEMORY: " << setw(10) << GetMemory()/(1048576.0));
    OutPut(" MB" << endl);

    // output of data for plotting
    if (!TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)
      OutputData2D(os, Output, i);

    // measure errors to known solution
    TDatabase::ParamDB->MEASURE_ERRORS = meas_errors;
    if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
      fesp[0] = USpaces[mg_level];
      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      u->GetErrors(Exact, 3, AllDerivatives, 4, SDFEMErrors,
        BilinearCoeffs, aux, 1, fesp, errors);
      delete aux;

      l2[i] = errors[0];
      h1[i] = errors[1];
      sd[i] = errors[2];
      energy[i] = sqrt(errors[1]*errors[1]/TDatabase::ParamDB->RE_NR + errors[0]*errors[0]);
      l_inf[i] = errors[3];
      if(TDatabase::ParamDB->DISCTYPE == LOCAL_PROJECTION)
      {
        if(TDatabase::ParamDB->LP_FULL_GRADIENT)
        {
          lp[i] = sqrt(h1[i]*h1[i]/TDatabase::ParamDB->RE_NR + l2[i]*l2[i] +
            0);
          /*UltraLocalError(u, Exact,
                       TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF,
                       TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT));*/
          //TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE));
        }

        if(TDatabase::ParamDB->LP_STREAMLINE)
        {
          OutPut("local projection stabilisation in streamline direction ");
          OutPut("is currently not available." << endl);
        }
      }

      if (i>FirstSolve)
      {
        OutPut( "L2: " << errors[0]<< " order " <<  log(l2[i-1]/l2[i])/ln2);
        OutPut( " H1-semi: " << errors[1] << " order " << log(h1[i-1]/h1[i])/ln2);
        OutPut( " energy: " << energy[i] << " order " << log(energy[i-1]/energy[i])/ln2 << endl);
        OutPut( "SD: " << errors[2] << endl);
        OutPut( "L_inf: " << errors[3] << " order " << log(l_inf[i-1]/l_inf[i])/ln2 << endl);

        if(TDatabase::ParamDB->DISCTYPE == LOCAL_PROJECTION)
        {
          OutPut( "LP: " << lp[i] << " order " << log(lp[i-1]/lp[i])/ln2 << endl);
        }
      }
      else
      {
        OutPut( "L2: " << errors[0]);
        OutPut( " H1-semi: " << errors[1]);
        OutPut( " energy: " << energy[i] << endl);
        OutPut( "SD: " << errors[2] << endl);
        OutPut( "L_inf: " << errors[3] << endl);
        if(TDatabase::ParamDB->DISCTYPE == LOCAL_PROJECTION)
        {
          OutPut( "LP: " << lp[i] << endl);
        }
      }
#ifdef __SMOOTH__
      ErrorInGivenPoints(UArray[mg_level]);
#endif
    }                                             // endif MEASURE_ERRORS
    // measure L1 error to known solution (check wether this routine is correct
    if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      //u->Interpolate(Exact);
      u->GetErrors(ExactNull, 1, ZeroDerivative, 1, L1Error,
        BilinearCoeffs, aux, 1, fesp, errors);
      delete aux;
      //OutPut(setw(20) << setprecision(20) << "L1: " << errors[0] << endl);
      OutPut("L1: " << errors[0] << endl);
    }                                             // endif MEASURE_ERRORS

    if(TDatabase::ParamDB->ESTIMATE_ERRORS)
    {
      ComputeErrorEstimate(coll, u,
        Coefficients[0], BoundaryConditions,
        BoundaryValues, eta_K,
        &maximal_local_error,
        estimated_global_error, l2[i], h1[i], N_Unknowns);

    }                                             // endif(ESTIMATE_ERRORS)

#ifdef  __TWO_INTERIOR_LAYERS__
    ComputeExtremalValues(N_Unknowns,sol,errors);
    OutPut("undershoots " << errors[0] << " overshoots " << errors[1] << endl);
    //OutPut(setprecision(3) << "min " << errors[0] << " & max " << errors[1]
    //  << " & min2 " << errors[2] << " & max2 " << errors[3] );
    //ComputeReferenceOutflowValues(score);
    //ComputeOutflowBoundaryErrorNew(u,errors);
    //OutPut(setprecision(3) << " & mino " << errors[0] << " & maxo " << errors[1]
    //  << " & midv " << errors[3]);
    /*ComputeDiagonal(mg_level, u);
    ComputeOutflowBoundaryError(u,errors);
    OutPut(endl);
    CheckMaximumPrinciple(sqmatrixA,sol,N_Active,errors);
        ComputeCutLines(i, u, errors, hmax);
        score[0] = 0;
        score[1] = 0.0908781;
        score[2] = 0.3788512;
        score[1] = 0.114518;
        score[2] = 0.385697;

    //smear_ref = 0.114518
    //width_ref = 0.385697

    OutPut(setprecision(3) << " & smear " << fabs((errors[0]+errors[1]-(score[0]+score[1]))/(score[0]+score[1]))  <<
    " & width " << fabs(errors[2]-score[2])/score[2] << endl);
    OutPut(setprecision(7) <<"ref: smear " << score[0]+score[1] << " width " << score[2] << " comp " << errors[0]+errors[1]
    << " " << errors[2] << endl);
    ComputeOutflowBoundary(mg_level, u);*/
#endif

#ifdef __BURMAN_HANSBO_2004__
    ComputeExtremalValues(N_Unknowns, sol,errors);
    OutPut(setprecision(4) << "min " << errors[0] << "& max " << errors[1]);
#endif

#ifdef __HMM_1986__
    //ComputeExtremalValues(N_Unknowns, sol,errors);
    //OutPut(setprecision(4) << "& min " << errors[0] << "& max " << errors[1]
    //  << " & " << errors[2] );
    ComputeCutLines(i, u, errors, hmax);
    ComputeGlobalErrors(u,errors, N_neum_to_diri, neum_to_diri);
    //    ComputeScores(errors,score);
    OutPut(endl << N_U << " oscint " << setprecision(4) << errors[1] << " oscexp " << errors[2]
      << " smearint " << errors[0] << " smearexp " << errors[3]);
    OutPut(" oscint_sold2 " << setprecision(4) << errors[4] <<
      " oscexp_sold2 " << setprecision(4) << errors[5] <<  endl);
    /* OutPut(" & " << setprecision(3) << errors[1] << " & " << setprecision(1)<< (int) score[1] <<
     " & " << setprecision(3) << errors[2] << " & " << setprecision(1)<< (int) score[2] <<
     " & " << setprecision(3) << errors[0] << " & " << setprecision(1) << score[0] <<
     " & " << setprecision(3) << errors[3] << " & " << setprecision(1) << score[3] <<
     " & " <<  setw(3) << setprecision(1) << score[1]+score[2] <<
     " & " <<  setw(3) << setprecision(1) << score[0]+score[3] <<
     " & " <<  setw(3) << setprecision(1) << score[0]+score[1]+score[2]+score[3] <<
     " & " << setw(4) << m << setprecision(7) << endl);*/
    //OutPut("L2 " << errors[0] << "& H1 " << errors[1] << "& nonlin " << m << "&");
#endif

#ifdef __PARABOLIC_LAYERS__
    ComputeCutLines(i, u, errors, hmax);
    OutPut("level " << i << setprecision(5) << " mesh parallel to x-axes & max(x=0.5) " << errors[0] << " & smear(x=0.5) " << -errors[1]
      << " & center " << errors[2] <<  " & explay " << errors[3] << endl);
    ComputeGlobalErrors(u,errors);
    OutPut(setprecision(5) << " & linfty(node) " << errors[0] << " & l2(node) " << errors[1]
      << " & l2(interior) " << errors[2] << endl);
    //ComputeGlobalErrorsIrr(u,errors);
    //OutPut(setprecision(5) << "IRR GRID  &para_val  " << errors[0] << " & para_der " << errors[1]
    //  << " &exp_der " << errors[2] << " &smear " << errors[3] << endl);
#endif

#ifdef __LEVEQUE__
    errors[0] = ComputeSherwoodNumber(u);
    // Peclet number
    errors[1] = 4*64*TDatabase::ParamDB->RE_NR;
    // Sherwood number
    errors[1] = 0.8075491 * pow( errors[1], 1.0/3.0);
    OutPut("level: " << i <<  " Sherwood number numerical " <<  -errors[0] << " theoretical " <<  errors[1] <<
      " error " << fabs(errors[0]+errors[1]) << endl);
#endif

#ifdef   __HEMKER1996__
    ComputeExtremalValues(N_Unknowns, sol,errors);
    OutPut(setprecision(4) << "ext min " << errors[0] << "& max " << errors[1]<<endl);
    //ComputeOutflowBoundary(mg_level, u);
    ComputeLocalExtrema(u,errors);
    OutPut(setprecision(3) << "circ_min " << errors[0] << "& circ_max " << errors[1]
      << " down_min " << errors[2] << "& down_max " << errors[3]<< endl);
    /* ComputeExtremalValues(N_Unknowns, sol,errors);
    OutPut(setprecision(4) << "min " << errors[0] << "& max " << errors[1]<<endl);
    ComputeOutflowBoundary(mg_level, u);*/
    ComputeLocalExtrema(u,errors);
    OutPut(setprecision(3) << "circ_min " << errors[0] << "& circ_max " << errors[1]
      << " down_min " << errors[2] << "& down_max " << errors[3]<< endl);
    // u->FindGradient(j*0.01,0.5,values);
    //  OutPut(j*0.01 << " " << 0.5 << " value: " << values[0]);
    // OutPut(" grad: " << values[1] << " " << values[2] << endl);

    eps = 1/TDatabase::ParamDB->RE_NR;
    val_x[0]=-1; val_x[1]=0; val_x[2]=1;
    val_x[3]=2; val_x[4]= 4; val_x[5]=8;
    val_y[0]=1+pow(eps,0.5); val_y[1]=1; val_y[2]=1-pow(eps,0.5); val_y[3]=0;

    for (iter=0;iter<4;iter++)
    {
      if (iter==0)
      {
        OutPut(" tab1: lev " << i  << "\\begin{tabular}{|l|l|l|l|l|l|l|}" << "\\hline " << endl);
        OutPut(" u(x,y;$\\varepsilon)$ & x=-1&x=0&x=1&x=2&x=4&x=8 \\\\ \\hline" << endl);
        OutPut(" y=1+$\\varepsilon^{\\frac{1}{2}}$  & ");
      }
      if (iter==1)  OutPut(" y=1 & ");
      if (iter==2)  OutPut("y=1-$\\varepsilon^{\\frac{1}{2}}$ & ");
      if (iter==3)  OutPut(" y=0 & ");
      for (n=0;n<6;n++)
      {
        u->FindGradient(val_x[n], val_y[iter], values);
        OutPut( setprecision(7) << values[0] );
        if(n!=5) OutPut( " & " );
      }
      OutPut(" \\\\ \\hline" << endl);

    }
    OutPut("\\end{tabular} \\\\[1em]" << endl);
    val_r[0]=1+eps; val_r[1]=1+pow(eps,2.0/3.0); val_r[2]=1+pow(eps,0.5);
    val_r[3]=1+pow(eps,1.0/3.0);
    val_y[0]=1; val_y[1]=1-pow(eps,0.5); val_y[2]=0;
    val_x[0]=-1; val_x[1]=0; val_x[2]=1;
    // Tabelle2
    for (iter=0;iter<3;iter++)
    {
      if (iter==0)
      {
        OutPut(" tab2: lev " << i );
        OutPut("\\begin{tabular}{|l|l|l|l|l|}" << "\\hline  " << endl);
        OutPut("u(x,y;$\\varepsilon)$ & r=1+$\\varepsilon$ & r=1+$\\varepsilon^{\\frac{2}{3}}$&r=1+$\\varepsilon^{\\frac{1}{2}}$&r=1+$\\varepsilon^{\\frac{1}{3}}$ \\\\ \\hline" << endl);
        OutPut(" y=1 & ");
      }
      if (iter==1) OutPut("y=1-$\\varepsilon^{\\frac{1}{2}}$ & ")
          if (iter==2) OutPut("y=0 & ");

      for (n=0;n<4;n++)
      {
        help=sqrt(val_r[n]*val_r[n]-val_y[iter]*val_y[iter]);
        u->FindGradient(help, val_y[iter], values);
        OutPut( setprecision(7) << values[0] );
        if(n!=3) OutPut( " & " );
        //OutPut("Level: " << i << setprecision(10) << " r: " << val_r[iter] << " x2.1: " << help << " y2.1: " << val_y[n] << " u2.1: " << values[0]<<endl);
      }
      OutPut(" \\\\ \\hline" << endl);

    }

    for (iter=0;iter<3;iter++)
    {
      if (iter==0) OutPut("x=-1 & ")
          if (iter==1) OutPut("x=0 & ")
            if (iter==2) OutPut("x=1 & ");

      for (n=0;n<4;n++)
      {
        help=sqrt(val_r[n]*val_r[n]-val_x[iter]*val_x[iter]);
        u->FindGradient(val_x[iter], help, values);
        OutPut( setprecision(7) << values[0] );
        if(n!=3) OutPut( " & " );
        // OutPut("Level: " << i << " r: " << val_r[iter] << " x2.2: " << val_x[n] << " y2.2: " << help << " u2.2: " << values[0]<<endl);
      }
      OutPut(" \\\\ \\hline" << endl);

    }
    OutPut("\\end{tabular}"<< endl);
    //if (mg_level > 0)
    //  {
    //    OutPut(mg_level<<endl);
    //   ComputeDifferenceToCoarseLevel(coll, coll_coarse, UArray[mg_level], UArray[mg_level-1]);
    // }
    coll_coarse = coll;
    ComputeCutLines_X(coll, UArray[mg_level],i);
    ComputeCutLines_Y(coll, UArray[mg_level],i);
    //ComputeCutLines_epsY(coll, UArray[mg_level],i);
    //ComputeCutLines_eps_radial(coll, UArray[mg_level],i);
    ComputeExtremalValues(N_Unknowns,sol,errors);
    OutPut("undershoots " << errors[0] << " overshoots " << errors[1] << endl);
#endif
#ifdef   __PW_LINEAR_RHS__
    /* ComputeExtremalValues(N_Unknowns, sol,errors);
    OutPut(setprecision(4) << "min " << errors[0] << "& max " << errors[1]<<endl);
    ComputeOutflowBoundary(mg_level, u);*/
    //ComputeCutLines(i, u, errors, hmax);
    ComputeLocalExtrema(u,errors);
    OutPut(setprecision(3) << " min_int " << errors[0] << "& osc_out " << errors[1] << endl);
#endif

    // remove data which will not be used later
    delete oldsol;
    if ((TDatabase::ParamDB->DISCTYPE==2)&&(TDatabase::ParamDB->SOLD_TYPE))
      delete defect;

    if ((mg_type==1)||(TDatabase::ParamDB->SOLVER_TYPE == AMG_SOLVE)
      ||(TDatabase::ParamDB->SOLVER_TYPE == DIRECT))
    {
      delete sqmatrixA;
      delete sqstructureA;
      delete rhs;
      if (mg_type==1)
        delete MGLevel;
    }                                             // end if (mg_type==1)

    old_sol = sol;
    old_u = u;
    old_u_space = velocity_space;

    OutPut("memory after: " << setw(10) << GetMemory() << endl);

  }                                               // endfor i

#ifdef   __POISSON_CONST__
  OutPut("here " << endl);
  ComputeValues(u);
#endif

  OutPut("used time: " << GetTime() << endl);
  OutPut("used bytes: " << GetMemory() << endl);
  CloseFiles();

  return 0;
}
