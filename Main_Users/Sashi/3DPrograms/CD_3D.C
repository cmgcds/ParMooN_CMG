// =======================================================================
//
// Purpose:     main program
//
// Author:      Volker John   August 2000
//
// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <DiscreteForm3D.h>
#include <FEDatabase3D.h>
#include <FESpace3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <AuxParam3D.h>
#include <QuadAffin.h>
#include <Solver.h>
#include <DirectSolver.h>
#include <Assemble3D.h>
#include <Output3D.h>
#include <LinAlg.h>

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
#include <Upwind3D.h>
#include <FEM_TVD_FCT.h>
#include <MultiGrid3D.h>
#include <MGLevel3D.h>
#include <ItMethod.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <JacobiIte.h>
#include <MultiGridScaIte.h>

#define AMG 0
#define GMG 1
#define DIRECT 2

// =======================================================================
// include current example
// =======================================================================
// #include "../Examples/CD_3D/Laplace.h"
// #include "../Examples/CD_3D/LaplaceBsp1_2.h"
//#include "../Examples/CD_3D/ThreeBoundaryLayers.h"
// #include "../Examples/CD_3D/Sphere.h"
// #include "../Examples/CD_3D/Constant1.h"
// #include "../Examples/CD_3D/Cylinder.h"
// #include "../Examples/CD_3D/Plane.h"
//#include "../Examples/CD_3D/CircularLayer.h"
//#include "../Examples/CD_3D/SkewConv.h"
//#include "../Examples/CD_3D/ParabolicLayers3D.h"
#include "../Examples/CD_3D/Hemker1996_3D.h"
//#include "../Examples/TCD_3D/bulk_compare.h"

int main(int argc, char* argv[])
{
  TDomain *Domain = new TDomain();
  TDatabase *Database = new TDatabase();
  TFEDatabase3D *FEDatabase = new TFEDatabase3D();
  TCollection *coll;
  TFESpace3D *velocity_space;
  TFESpace3D *velocity_space_low;
  TFESpace3D *old_u_space, *shock_capturing_space;
  TFESpace3D **USpaces;
  TOutput3D *Output;

  double *rhs, *sol, *oldsol, *defect, *itmethod_sol, *itmethod_rhs;
  double *rhs_low, *sol_low, *old_sol, *update;
  int i,j,k,l,m, low;
  int N_U, N_Unknowns, N_RHSs, N_SqMat;
  int N_Active, N_NonActive, N_U_low, N_Unknowns_low;
  double *l2, *h1, *sd, *l_inf, *l2smooth, *h1smooth;
  char *PRM, *GEO;
  int LEVELS;
  int ret, ORDER;
  double errors[6];
  double t1, t2, res, res2, oldres, solver_time, hmin, hmax, oldres_stepm1;
  double nonlin_min_res;
  int N_LinIter, first_damp,linite, compute_matrix_D;

  std::ostringstream os;
  char *PsBaseName, *GrapeBaseName, *GnuBaseName, *ReadGrapeBaseName, *GMVBaseName;
  char *VTKBaseName;

  // TFEFunction3D **AuxFEFunctArray, **AuxFEVectFunctArray;
  TFEFunction3D *u, **UArray, *u_low, *old_u;
  TFESpace3D *fesp[2], *ferhs[3];

  TAuxParam3D *aux;

  TSquareStructure3D *sqstructureA;
  TSquareMatrix3D *sqmatrixA, *SQMATRICES[1];
  TSquareMatrix3D **MatricesA;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;

  TSquareStructure3D *sqstructureA_low;
  TSquareMatrix3D *sqmatrixA_low;
  double **RhsArray;

  TMGLevel3D *MGLevel, *MGLevel_low;
  TMultiGrid3D *MG;
  TItMethod *itmethod, *prec;

  double *RHSs[1];
  int *N_Uarray;

  char UString[] = "u";
  char CdString[] = "Conv-Diff";
  char GalString[] = "Galerkin";
  char SDFEMString[] = "SDFEM";
  char UpwString[] = "Upwind";
  char SCDFEMString[] = "SCDFEM";
  char SC_DC_CDString[] = "SC_DC_CD";
  char SC_MBEString[] = "SC_MBE";
  char SC_1String[] = "SC_1";
  char SC_2String[] = "SC_2";
  char Name[] = "name";
  char Description[] = "description";
  char Readin[] = "readin.dat";

  TDiscreteForm3D *DiscreteForm;

  int FirstSolve;
  int N_Paramters=2,n_aux;
  double Parameters[2];
  int mg_level,mg_type;
  int *neum_to_diri, N_neum_to_diri = 0;
  double *neum_to_diri_x, *neum_to_diri_y, *neum_to_diri_z;

  double *sc_params, residual, omega, omega_max = 1.0, res_norm_min, lin_red;
  double  *lump_mass, *matrix_D_Entries;
  double omega_min = 0.1, *rhs_edge;
  int N_sc, sc_N_params, max_it, sold_parameter_type;

  TFEFunction3D *sc_params_fe[2], *fefct[3];
  TFEVectFunct3D *sc_params_vect;
  DefectProc *Defect;

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

  ManipulateFct3D *manipulate;
  if (TDatabase::ParamDB->SDFEM_NORM_B==0)
    manipulate = linfb;
  else
    manipulate = ave_l2b_quad_points;

  TDiscreteForm3D *DiscreteFormGalerkin = new TDiscreteForm3D
    (CdString, GalString, N_Terms, Derivatives, SpacesNumbers,
    N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace,
    BilinearAssemble, BilinearCoeffs, NULL);

  TDiscreteForm3D *DiscreteFormSDFEM = new TDiscreteForm3D
    (CdString, SDFEMString, N_Terms_SD, Derivatives_SD, SpacesNumbers_SD,
    N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace,
    BilinearAssemble_SD, BilinearCoeffs, manipulate);

  TDiscreteForm3D *DiscreteFormUpwind = new TDiscreteForm3D
    (CdString, UpwString, N_Terms, Derivatives, SpacesNumbers,
    N_Matrices, N_Rhs, RowSpace, ColumnSpace, RhsSpace,
    BilinearAssemble_UPW2, BilinearCoeffs, NULL);

  TDiscreteForm3D *DiscreteFormSOLD = new TDiscreteForm3D
    (CdString, SC_2String, N_Terms_SD, Derivatives_SD, SpacesNumbers_SD,
    N_Matrices, 0 , RowSpace, ColumnSpace, NULL,
    BilinearAssemble_SOLD, BilinearCoeffs, manipulate);

  TDiscreteForm3D *DiscreteFormSOLD_Orthogonal = new TDiscreteForm3D
    (CdString, SC_2String, N_Terms_SD, Derivatives_SD, SpacesNumbers_SD,
    N_Matrices, 0 , RowSpace, ColumnSpace, NULL,
    BilinearAssemble_SOLD_Orthogonal, BilinearCoeffs, manipulate);

  BoundCondFunct3D *BoundaryConditions[1] = { BoundCondition };
  BoundValueFunct3D *BoundaryValues[2] = { BoundValue, BoundValue };

  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->WriteParamDB(argv[0]);
  ExampleFile();

  //======================================================================
  // copy read parameters into local variables
  //======================================================================

  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;

  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  GrapeBaseName = TDatabase::ParamDB->GRAPEBASENAME;
  GnuBaseName = TDatabase::ParamDB->GNUBASENAME;
  GMVBaseName = TDatabase::ParamDB->GMVBASENAME;
  VTKBaseName = TDatabase::ParamDB->VTKBASENAME;
  ReadGrapeBaseName = TDatabase::ParamDB->READGRAPEBASENAME;

  sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;

  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SCALAR;
  if (TDatabase::ParamDB->SOLVER_TYPE==AMG||
    TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
    mg_type = 0;
  if (mg_type)
    mg_level = 0;
  else
    mg_level = -1;
  LEVELS = TDatabase::ParamDB->LEVELS;
  l2 = new double[LEVELS+1];
  h1 = new double[LEVELS+1];
  sd = new double[LEVELS+1];
  l2smooth = new double[LEVELS+1];
  h1smooth = new double[LEVELS+1];
  l_inf = new double[LEVELS+1];

  UArray = new TFEFunction3D*[LEVELS+1];
  RhsArray = new double* [LEVELS+1];
  N_Uarray = new int[LEVELS+1];

  USpaces = new TFESpace3D*[LEVELS+1];
  MatricesA = new TSquareMatrix3D*[LEVELS+1];

  //======================================================================
  // read boundary parameterization and initialize coarse grid
  //======================================================================
  Domain->Init(PRM, GEO);

  // refine up to user defined coarsest level
  for(i=0;i<TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR;i++)
  {
    Domain->RegRefineAll();
    coll=Domain->GetCollection(It_Finest, 0);
    Domain->MakeBdParamsConsistent(coll);
    delete coll;
  }

  // initialize solver parameters
  Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
  Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;

  if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    i=1;
    MG = new TMultiGrid3D(i, N_Paramters, Parameters);
  }

  FirstSolve = TDatabase::ParamDB->SC_FIRST_SOLUTION_LEVEL_SCALAR;

  if(TDatabase::ParamDB->P7 == 13579)
  {
    Domain->RegRefineAll();
    coll=Domain->GetCollection(It_Finest, 0);
    coll->GetCell(0)->GetVertex(6)->SetCoords(
      TDatabase::ParamDB->P4,
      TDatabase::ParamDB->P5,
      TDatabase::ParamDB->P6);
    OutPut("midpoint moved" << endl);
  }

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
    solver_time = 0.0;
    N_LinIter = 0;
    OutPut("memory before: " << setw(10) << GetMemory() << endl);

    // refine grid if level is greater than 0
    if (i)
      Domain->RegRefineAll();

    coll=Domain->GetCollection(It_Finest, 0);
#ifdef __HEMKER1996__
    SetPeriodicFaceJoints(coll);
#endif
    Domain->MakeBdParamsConsistent(coll);
    Output = new TOutput3D(1, 1, 1, 1,Domain);

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
        -1);
    }
    // get spaces of high order disc on finest geo grid
    if ((i>=FirstSolve)||(mg_type==0))
      velocity_space =  new TFESpace3D(coll, Name,
        Description,
        BoundCondition, ORDER);

    // build fespace hierarchy
    // set values and pointers for low order fe space
    if (mg_type==1)
    {
      USpaces[i] = velocity_space_low;
      N_U_low = velocity_space_low->GetN_DegreesOfFreedom();
      N_Uarray[i] = velocity_space_low->GetN_DegreesOfFreedom();
    }
    // set values and pointers for high order fe space
    if ((i>=FirstSolve)||(mg_type==0))
    {
      USpaces[mg_level] = velocity_space;
      N_U = velocity_space->GetN_DegreesOfFreedom();
      N_Uarray[mg_level] = N_U;
      N_Active = velocity_space->GetActiveBound();
      N_NonActive = N_U - N_Active;
    }

    // SOLD schemes

    if (TDatabase::ParamDB->SOLD_TYPE)
    {
      // allocate piecewise constant finite element space
      shock_capturing_space = new TFESpace3D(coll, Name, Description,
        BoundCondition, 0);
      sc_N_params = 2;
      N_sc = shock_capturing_space->GetN_DegreesOfFreedom();
       if (sc_N_params)
      {
        sc_params = new double[sc_N_params* N_sc];
        sc_params_vect = new TFEVectFunct3D(shock_capturing_space, UString, UString, sc_params,
          N_sc, sc_N_params);
        for (m=0;m<sc_N_params;m++)
          sc_params_fe[m] = sc_params_vect->GetComponent(m);
      }
    }

    // build matrices for high order disc
    if ((i>=FirstSolve)||(mg_type==0))
    {
      // matrix structures
      sqstructureA = new TSquareStructure3D(velocity_space);
      sqstructureA->Sort();
      // allocate matrices
      sqmatrixA = new TSquareMatrix3D(sqstructureA);
      OutPut("level "<< i << " number of matrix entries " << sqmatrixA->GetN_Entries() << endl);
      MatricesA[mg_level] = sqmatrixA;
    }

    // FEM TVD methos
    if (TDatabase::ParamDB->SOLD_TYPE)
    {
	if ((sold_parameter_type ==  FEM_TVD)&&((i>=FirstSolve)||(mg_type==0)))
      {
        matrix_D_Entries = new double[sqmatrixA->GetN_Entries()];
        memset(matrix_D_Entries, 0, sqmatrixA->GetN_Entries()*SizeOfDouble);
        rhs_edge = new double[N_U];
        memset(rhs_edge, 0, N_U * SizeOfDouble);
      }
    }

    // build matrices for low order disc
    if (mg_type==1)
    {
      // matrix structures
      sqstructureA_low = new TSquareStructure3D(velocity_space_low);
      sqstructureA_low->Sort();
      sqmatrixA_low = new TSquareMatrix3D(sqstructureA_low);
      OutPut("level "<< i << " number of matrix entries " << sqmatrixA_low->GetN_Entries() << endl);
      MatricesA[i] = sqmatrixA_low;
    }                                             // end if (mg_type==1)

    N_Unknowns = N_U;
    if (mg_type==1)
      N_Unknowns_low = N_U_low;

    OutPut("dof all      : "<< setw(10) << N_Unknowns  << endl);
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

    coll->GetHminHmax(&hmin,&hmax);
    OutPut("h_min : " << hmin << " h_max : " << hmax << endl);

    // high order disc
    if ((i>=FirstSolve)||(mg_type==0))
    {
      rhs = new double[N_Unknowns];
      memset(rhs, 0, N_Unknowns*SizeOfDouble);
      RhsArray[mg_level] = rhs;
      sol = new double[N_Unknowns];
      oldsol = new double[N_Unknowns];
      update = new double[N_Unknowns];
      memset(sol, 0, N_Unknowns*SizeOfDouble);
      memset(oldsol, 0, N_Unknowns*SizeOfDouble);
      memset(update, 0, N_Unknowns*SizeOfDouble);
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
        if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR)
          || (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR))
          n_aux= 4;
        else
          n_aux= 2;
        // build fe multigrid levels
        // low order disc
        if (mg_type==1)
        {
          MGLevel_low = new TMGLevel3D(i, sqmatrixA_low,
            rhs_low,
            sol_low,
            n_aux, NULL);
          if (i==0)
            MG->AddLevel(MGLevel_low);
          else
            MG->ReplaceLevel(i,MGLevel_low);
        }
        // high order disc
        if ((i>=FirstSolve)||(mg_type==0))
        {
          MGLevel = new TMGLevel3D(mg_level, sqmatrixA,
            rhs, sol,
            n_aux, NULL);
          MG->AddLevel(MGLevel);
        }
        break;
    }

    // build new fe functions
    // high order fe space
    if ((i>=FirstSolve)||(mg_type==0))
    {
      u = new TFEFunction3D(velocity_space, UString, UString, sol, N_U);
      UArray[mg_level] = u;
    }

    // low order fe space
    if (mg_type==1)
    {
      u_low = new TFEFunction3D(velocity_space_low, UString, UString, sol_low, N_U_low);
      UArray[i] = u_low;
    }

    // CHECK THIS !!!old_u_space
    if ((i>=FirstSolve)||(mg_type==0))
    {
      Output->AddFEFunction(u);
    }

    // prolongation, to get a good starting iterate
    if(i && i>FirstSolve)
    {
      if (!((TDatabase::ParamDB->SOLVER_TYPE==AMG)&&
        // for assessment of LCD vs. GMRES, can be cancelled after the assessment
        ((TDatabase::ParamDB->SC_SOLVER_SCALAR==14)||(TDatabase::ParamDB->SC_SOLVER_SCALAR==19))))
        Prolongate(old_u_space, USpaces[mg_level],
          old_u->GetValues(), UArray[mg_level]->GetValues(),
          oldsol);

      if ((TDatabase::ParamDB->SOLVER_TYPE==AMG)||(TDatabase::ParamDB->SOLVER_TYPE==DIRECT))
      {
        delete USpaces[i-1];
        delete UArray[i-1]->GetValues();
      }
      // copy current solution for assembling the nonlinear terms
      memcpy(oldsol,sol, N_Unknowns*SizeOfDouble);

      if (mg_type==1)
      {
        delete old_sol;
        delete old_u;
        delete old_u_space;
      }
    }                                             // end of prolongate
    // restrict solution to all grids
    if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
      MG->RestrictToAllGrids();

    // if no solution on this grid, continue
    if(FirstSolve>i)
      continue;

    /*if(TDatabase::ParamDB->READ_GRAPE_FILE)
    {
      AuxFEFunctArray = new TFEFunction3D*[1];
      AuxFEFunctArray[0] = PArray[mg_level];
      AuxFEVectFunctArray = new TFEVectFunct3D*[1];
      AuxFEVectFunctArray[0] = UArray[mg_level];
      ReadGrapeFile3D(ReadGrapeBaseName, 1 , 1 , AuxFEFunctArray,AuxFEVectFunctArray);
      TDatabase::ParamDB->READ_GRAPE_FILE = 0;
      OutPut("u " << Ddot(N_U,sol,sol)<< endl);
      memcpy(oldsol,sol, N_Unknowns*SizeOfDouble);
      if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
    MG->RestrictToAllGrids();
    }*/

    // build the discretizations
    for(k=low;k<=mg_level;k++)
    {
      rhs = RhsArray[k];
      N_U = N_Uarray[k];
      N_Active = USpaces[k]->GetActiveBound();
      N_NonActive = N_U - N_Active;

      RHSs[0] = rhs;
      memset(rhs, 0, N_U*SizeOfDouble);

      fesp[0] = USpaces[k];
      ferhs[0] = USpaces[k];
     // find discrete form
      if ((mg_type==1) && (k<i+1))
      {
        DiscreteForm = DiscreteFormUpwind;
	OutPut("UPWIND"<<endl);
      }
      else
        switch(TDatabase::ParamDB->DISCTYPE)
        {
          case GALERKIN:
            DiscreteForm = DiscreteFormGalerkin;
            break;

          case SDFEM:
            DiscreteForm = DiscreteFormSDFEM;
            break;

          case UPWIND:
	      DiscreteForm = DiscreteFormUpwind;
	      break;

        default:
          Error("Unknown DISCTYPE" << endl);
          return -1;
      }

      // initialize matrices
      SQMATRICES[0] = MatricesA[k];
      SQMATRICES[0]->Reset();
      aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

      SetPolynomialDegree();
      // assemble
       // homogeneous Neumann boundary conditions  
      if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD)
      {
	  BoundaryValues[0] = BoundValue_FEM_FCT;
      }
     Assemble3D(1, fesp,
        1, SQMATRICES,
        0, NULL,
        1, RHSs, ferhs,
        DiscreteForm,
        BoundaryConditions,
        BoundaryValues,
        aux);

      // reset
      if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD)
      {
	  BoundaryValues[0] = BoundaryValues[1];
      }

      if (DiscreteForm == DiscreteFormUpwind)
      {
        UpwindForConvDiff(SQMATRICES[0],RHSs[0],
          fesp[0],DiscreteFormUpwind);
      }                                           // endif
      delete aux;

      if ((TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD)&&(k==mg_level))
      {
#ifdef __HEMKER1996__
        CheckWrongNeumannNodes(coll, USpaces[k], N_neum_to_diri, neum_to_diri,
          neum_to_diri_x, neum_to_diri_y, neum_to_diri_z);
#endif
        if (N_neum_to_diri)
          SetDirichletNodesFromNeumannNodes(SQMATRICES, RHSs[0], UArray[k]->GetValues(),
            N_neum_to_diri, neum_to_diri,
            neum_to_diri_x, neum_to_diri_y, neum_to_diri_z,
            BoundValue);
      }
    }                                             // endfor, assembling done

    // set rhs for Dirichlet nodes
    memcpy(sol+N_Active, rhs+N_Active, N_NonActive*SizeOfDouble);

    // compute defect
    defect = new double[N_Unknowns];
    memset(defect,0,N_Unknowns*SizeOfDouble);
    // solve system
    switch(TDatabase::ParamDB->SOLVER_TYPE)
    {
      case DIRECT:
        t1 = GetTime();
        DirectSolver(MatricesA[mg_level], RhsArray[mg_level], sol);
        t2 = GetTime();
        OutPut( "time for direct solving: " << t2-t1 << endl);
        OutPut("solution " << sqrt(Ddot(N_Unknowns,sol,sol)) << endl);
        break;

      case AMG:
        t1 = GetTime();
        Solver(MatricesA[mg_level], RhsArray[mg_level], sol);
        t2 = GetTime();
        OutPut( "time for AMG solving: " << t2-t1 << endl);
        OutPut("solution " << sqrt(Ddot(N_Unknowns,sol,sol)) << endl);
        break;

      case GMG:
        t1 = GetTime();
        // build preconditioner
        switch (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR)
        {
          case 1:
            prec = new TJacobiIte(MatVect_Scalar, Defect_Scalar, NULL,
              0, N_Unknowns, 1);
            break;
          case 5:
            prec = new TMultiGridScaIte(MatVect_Scalar, Defect_Scalar, NULL,
              0, N_Unknowns, MG, 1);
            break;
          default:
            OutPut("Unknown preconditioner !!!" << endl);
            exit(4711);
        }
        switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
        {
          case 11:
            itmethod = new TFixedPointIte(MatVect_Scalar, Defect_Scalar, prec,
              0, N_Unknowns, 1);
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
            {
              itmethod_sol = new double[N_Unknowns];
              itmethod_rhs = new double[N_Unknowns];
              memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
              memcpy(itmethod_rhs, rhs, N_Unknowns*SizeOfDouble);
            }
            else
            {
              itmethod_sol = sol;
              itmethod_rhs = rhs;
            }
            break;
          case 16:
            itmethod = new TFgmresIte(MatVect_Scalar, Defect_Scalar, prec,
              0, N_Unknowns, 1);
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
            {
              itmethod_sol = new double[N_Unknowns];
              itmethod_rhs = new double[N_Unknowns];
              memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
              memcpy(itmethod_rhs, rhs, N_Unknowns*SizeOfDouble);
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
        // solve linear system
        itmethod->Iterate(sqmatrices,NULL,itmethod_sol,itmethod_rhs);
        switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
        {
          case 11:
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
            {
              memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
              memcpy(rhs, itmethod_rhs, N_Unknowns*SizeOfDouble);
            }
            break;
          case 16:
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
            {
              memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
              memcpy(rhs, itmethod_rhs, N_Unknowns*SizeOfDouble);
            }
            break;
        }
    }

    // ************************************************************* //
    // end of iteration
    // ************************************************************* //

    m = 0;
    linite = 0;
    compute_matrix_D = 1;
    nonlin_min_res = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR;

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
	  omega = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SCALAR;
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
            aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
            u->GetMeshCellParams(ExactNull, 5, Derivatives_SD, 2, Parameters_DC_CD,
              BilinearCoeffs, aux, 1, fesp, errors,
              sc_params);
            for (k=0;k<2*N_sc;k++)
              sc_params[k] = sqrt(sc_params[k]);
            delete aux;

            // set aux object
            aux =  new TAuxParam3D(SOLD_N_FESpaces, SOLD_N_Fct, SOLD_N_ParamFct,
              SOLD_N_FEValues, fesp, fefct,
              SOLD_Fct,
              SOLD_FEFctIndex, SOLD_FEMultiIndex,
              SOLD_N_Params, SOLD_BeginParam);

            // set discrete form
            if (TDatabase::ParamDB->SOLD_TYPE==1)
              DiscreteForm = DiscreteFormSOLD;
            else
              DiscreteForm = DiscreteFormSOLD_Orthogonal;
            //if (sold_parameter_type == MH_Kno06)
            //  DiscreteForm = DiscreteFormMH_Kno06;
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
              exit(1);
              //DiscreteForm = DiscreteFormRhsLP96;
            }
          }

          // assemble
          Assemble3D(2, fesp,
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
            exit(4711);
            //if(TDatabase::ParamDB->LP_FULL_GRADIENT)
            //	  UltraLocalProjection(MatricesA[mg_level], FALSE, Coefficients[0]);

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
            exit(1);
            Daxpy(N_Active, 1, rhs_edge, rhs);

            memset(rhs_edge, 0,  N_Unknowns*SizeOfDouble);
            //EdgeStabilization(USpaces[mg_level], UArray[mg_level],
            //  Coefficients[0], rhs_edge, 0, NULL, NULL);
            // subtract rhs_edge from rhs
            Daxpy(N_Active, -1, rhs_edge, rhs);
          }

          if (sold_parameter_type == MH_Kno06)
          {
            // restore rhs
            exit(1);
            Daxpy(N_Active, -1, rhs_edge, rhs);

            memset(rhs_edge, 0,  N_Unknowns*SizeOfDouble);

            //MizukamiHughes(sqmatrixA, rhs_edge, USpaces[mg_level],
            //  UArray[mg_level], Coefficients[0],
            //  BoundCondition);

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
              neum_to_diri_x, neum_to_diri_y, neum_to_diri_z,
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
        switch(TDatabase::ParamDB->SOLVER_TYPE)
        {
          case DIRECT:
            t1 = GetTime();
            DirectSolver(MatricesA[mg_level], RhsArray[mg_level], sol);
            t2 = GetTime();
            OutPut( "time for direct solving: " << t2-t1 << endl);
            OutPut("solution " << sqrt(Ddot(N_Unknowns,sol,sol)) << endl);
            break;

          case AMG:
            t1 = GetTime();
            Solver(MatricesA[mg_level], RhsArray[mg_level], sol);
            t2 = GetTime();
            OutPut( "time for AMG solving: " << t2-t1 << endl);
            OutPut("solution " << sqrt(Ddot(N_Unknowns,sol,sol)) << endl);
            break;

          case GMG:
            t1 = GetTime();
            // build preconditioner
            /* switch (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR)
            {
              case 1:
                prec = new TJacobiIte(MatVect_Scalar, Defect_Scalar, NULL,
                  0, N_Unknowns, 1);
                break;
              case 5:
                prec = new TMultiGridScaIte(MatVect_Scalar, Defect_Scalar, NULL,
                  0, N_Unknowns, MG, 0);
                break;
              default:
            OutPut("Unknown preconditioner !!!" << endl);
            exit(4711);
            }*/
            delete itmethod;
            switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
            {
              case 11:
                itmethod = new TFixedPointIte(MatVect_Scalar, Defect_Scalar, prec,
                  0, N_Unknowns, 1);
                if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
                {
                  //itmethod_sol = new double[N_Unknowns];
                  //itmethod_rhs = new double[N_Unknowns];
                  memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
                  memcpy(itmethod_rhs, rhs, N_Unknowns*SizeOfDouble);
                }
                else
                {
                  itmethod_sol = sol;
                  itmethod_rhs = rhs;
                }
                break;
              case 16:
                itmethod = new TFgmresIte(MatVect_Scalar, Defect_Scalar, prec,
                  0, N_Unknowns, 1);
                if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
                {
                  //itmethod_sol = new double[N_Unknowns];
                  //itmethod_rhs = new double[N_Unknowns];
                  memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
                  memcpy(itmethod_rhs, rhs, N_Unknowns*SizeOfDouble);
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
            // solve linear system
            linite += itmethod->Iterate(sqmatrices,NULL,itmethod_sol,itmethod_rhs);
            switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
            {
              case 11:
                if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
                {
                  memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
                  memcpy(rhs, itmethod_rhs, N_Unknowns*SizeOfDouble);
                }
                break;
              case 16:
                if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
                {
                  memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
                  memcpy(rhs, itmethod_rhs, N_Unknowns*SizeOfDouble);
                }
                break;
            }
        }                                         // end solver

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

    OutPut("MEMORY: " << setw(10) << GetMemory()/(1048576.0));
    OutPut(" MB" << endl);

    if(TDatabase::ParamDB->WRITE_GRAPE)
    {
      os.seekp(std::ios::beg);
      os << GrapeBaseName << i << ".dat" << ends;
      Output->WriteGrape(os.str().c_str());
    }

    if(TDatabase::ParamDB->WRITE_GMV)
    {
      os.seekp(std::ios::beg);
      os << GMVBaseName << i << ".gmv" << ends;
      Output->WriteGMV(os.str().c_str());
    }

    if (TDatabase::ParamDB->WRITE_VTK)
    {
      os.seekp(std::ios::beg);
      os << VTKBaseName << i << ".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
    }

    if(TDatabase::ParamDB->WRITE_GNU)
    {
      os.seekp(std::ios::beg);
      os << GnuBaseName << i << ".gnu" << ends;
      //Output->WriteGnuplot(os.str().c_str());
    }

    // measure errors to known solution
    if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
      aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      // for(j=0;j<N_Unknowns;j++)
      //   cout << setw(3) << j << setw(30) << sol[j] << endl;
      // u->Interpolate(Exact);
      // k = 0;
      // l = velocity_space->GetN_ActiveDegrees();
      // for(j=0;j<N_Unknowns;j++)
      // {
      //   if(fabs(sol[j])>1e-10 && j>=l) k++;
      //   cout << setw(3) << j << setw(30) << sol[j] << endl;
      // }
      // if(k)
      // {
      //   cout << "Number of non-zero Dirichlet values: " << k << endl;
      // }
      // memset(sol+l, 0, (N_Unknowns-l)*SizeOfDouble);
      u->GetErrors(Exact, 4, AllDerivatives, 2, L2H1Errors,
        BilinearCoeffs, aux, 1, fesp, errors);
      delete aux;

      l2[i] = errors[0];
      h1[i] = errors[1];
      // sd[i] = errors[2];
      if (i>FirstSolve)
      {
        OutPut( "L2: " << errors[0]<< " order " <<  log(l2[i-1]/l2[i])/ln2 << endl);
        OutPut( "H1-semi: " << errors[1] << " order " << log(h1[i-1]/h1[i])/ln2 << endl);
        // OutPut( "SD: " << errors[2] << endl);
      }
      else
      {
        OutPut( "L2: " << errors[0] << endl);
        OutPut( "H1-semi: " << errors[1] << endl);
        // OutPut( "SD: " << errors[2] << endl);
      }
    }                                             // endif MEASURE_ERRORS

    if ((TDatabase::ParamDB->MEASURE_ERRORS)&&(TDatabase::ParamDB->P5==123456789))
    {
      aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      u->GetErrors(Exact, 4, AllDerivatives, 2, L2H1ErrorsSmooth,
        BilinearCoeffs, aux, 1, fesp, errors);

      l2smooth[i] = errors[0];
      h1smooth[i] = errors[1];
      // sd[i] = errors[2];
      if (i>FirstSolve)
      {
        OutPut( "smooth part: L2: " << errors[0]<< " order " <<  log(l2smooth[i-1]/l2smooth[i])/ln2 << endl);
        OutPut( "smooth part: H1-semi: " << errors[1] << " order " << log(h1smooth[i-1]/h1smooth[i])/ln2 << endl);
        // OutPut( "SD: " << errors[2] << endl);
      }
      else
      {
        OutPut( "smooth part: L2: " << errors[0] << endl);
        OutPut( "smooth part: H1-semi: " << errors[1] << endl);
        // OutPut( "SD: " << errors[2] << endl);
      }

      delete aux;
    }                                             // endif MEASURE_ERRORS

    /*
    #ifdef __CIRCULAR_LAYER__
        ComputeExtremalValues(N_Unknowns, sol,errors);
        OutPut(setprecision(3) << errors[0] << "& " << errors[1] << "&bdry ");
        ComputeOutflowBoundary(mg_level, u, errors);
        OutPut(setprecision(3) << errors[0] << "& " << errors[1]);
        CheckMaximumPrinciple(sqmatrixA,sol,N_Active,errors);
        OutPut(endl);
    #endif
    */
#ifdef __SKEW_CONV__
    ComputeExtremalValues(N_Unknowns, sol,errors);
    OutPut("minmax " << setprecision(3) << errors[0] << "& " << errors[1] << " ");
    ComputeCutLines(i, u, errors, hmax);
    OutPut(" smear " << setprecision(3) << errors[0] << endl);
#endif

#ifdef   __HEMKER1996__
    ComputeExtremalValues(N_Unknowns, sol,errors);
    OutPut(setprecision(4) << "ext min " << errors[0] << "& max " << errors[1]<<endl);
#endif

    // remove data which will not be used later
    delete oldsol;
    delete defect;

    if ((mg_type==1)||(TDatabase::ParamDB->SOLVER_TYPE == AMG)
      ||(TDatabase::ParamDB->SOLVER_TYPE == DIRECT))
    {
      delete sqmatrixA;
      delete sqstructureA;
      delete rhs;
      if (mg_type==1)
        delete MGLevel;
    }                                             // end if (mg_type==1)

    if (TDatabase::ParamDB->SOLVER_TYPE == GMG)
    {
      delete prec;
      delete itmethod;
      if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR==5)
      {
        delete itmethod_sol;
        delete itmethod_rhs;
      }
    }

    old_sol = sol;
    old_u = u;
    old_u_space = velocity_space;

    OutPut("memory after: " << setw(10) << GetMemory() << endl);

  }                                               // endfor i

  OutPut("used time: " << GetTime() << endl);
  OutPut("used bytes: " << GetMemory() << endl);
  CloseFiles();

  return 0;
}
