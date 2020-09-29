// =======================================================================
// @(#)MainRoutines.C
//
// Purpose: contains routines which are called from the main program 
//
// Author: Volker John 
//
// History: start of implementation 22.09.2009
//
// =======================================================================

#include <Constants.h>
#include <Database.h>
#include <AuxParam2D.h>
#include <DirectSolver.h>
#include <Solver.h>
#include <ItMethod.h>
#include <FEDatabase2D.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <JacobiIte.h>
#include <MultiGridScaIte.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>
#include <Output2D.h>
#include <ConvDiff2D_Routines.h>
#include <Assemble2D.h>
#include <MainUtilities.h>
#include <MainRoutines2D.h>
#include <Upwind.h>
#include <DiscreteForm2D.h>
#include <LocalProjection.h>

#include <stdlib.h>
#include <string.h>
#include <sstream>
#ifndef __MAC64__
#include <malloc.h>
#endif

/******************************************************************************/
// SetParametersCDAdapt2D()
// sets parameters of the data base for the main program CDAdapt2D.C
/******************************************************************************/

void SetParametersCDAdapt2D()
{
  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE== KLR02_3)
    TDatabase::ParamDB->SOLD_S = 0;
  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE== LP96)
  {
    OutPut("SOLD_PARAMETER_TYPE == LP96 should be used with higher quadrature rule,"<<endl);
    OutPut("since right hand side is in general not linear !!!"<<endl);
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
  

  if ((TDatabase::ParamDB->SDFEM_TYPE == 100)&&(!TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)&&
      (TDatabase::ParamDB->DISCTYPE==SDFEM))
  {
      TDatabase::ParamDB->SDFEM_TYPE = 2;
      OutPut("Changed Database::ParamDB->SDFEM_TYPE to " << TDatabase::ParamDB->SDFEM_TYPE
	     << " since no adjoint problem is solved !!! "<<endl);
  } 
  if ((TDatabase::ParamDB->SDFEM_TYPE != 100)&&(TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)&&
      (TDatabase::ParamDB->DISCTYPE==SDFEM))
  {
      TDatabase::ParamDB->SDFEM_TYPE = 100;
      OutPut("Changed Database::ParamDB->SDFEM_TYPE to " << TDatabase::ParamDB->SDFEM_TYPE
	     << " since adjoint problem is solved !!! "<<endl);
  } 
  if ((TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM==4))
    TDatabase::ParamDB->INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS = 1;
  // SUPG 
  if ((TDatabase::ParamDB->DISCTYPE==SDFEM)&&(TDatabase::ParamDB->SOLD_TYPE==0))
  {
      // this excludes some not wished side effects
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 0;
  } 
  if  (!(TDatabase::ParamDB->DISCTYPE==SDFEM))
    {
      TDatabase::ParamDB->SOLD_TYPE = 0;
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE =0;
    }
  if  ((TDatabase::ParamDB->DISCTYPE==SDFEM)&&(TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD))
    {
      TDatabase::ParamDB->SDFEM_TYPE = 0;
      TDatabase::ParamDB->DELTA0 =  TDatabase::ParamDB->DELTA1 = 0;
      OutPut("FEM-TVD: switched stabilization off!" << endl);
    }

  TDatabase::ParamDB->NSTYPE = 0;

  if (TDatabase::ParamDB->DISCTYPE==CIP)
    {
      TDatabase::ParamDB->DISCTYPE=GALERKIN;
      TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS = 1;
    }
  if (TDatabase::ParamDB->DISCTYPE==DG)
    {
      TDatabase::ParamDB->DISCTYPE=GALERKIN;
      TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS = 2;
      if ( TDatabase::ParamDB->ANSATZ_ORDER < 10)
	TDatabase::ParamDB->ANSATZ_ORDER = -TDatabase::ParamDB->ANSATZ_ORDER-10;
      else 
	// P elements on quads
	TDatabase::ParamDB->ANSATZ_ORDER = -10*TDatabase::ParamDB->ANSATZ_ORDER;
      if (TDatabase::ParamDB->ESTIMATE_ERRORS)
	{
	  TDatabase::ParamDB->ESTIMATE_ERRORS = 0;
	  OutPut("Error estimation does not work for DG !!!"<< endl);
	}
    }
}

/******************************************************************************/
// Solver
// solves linear system
// output in sol
// scalar problems: ns_type == 0
/******************************************************************************/
void Solver(TSquareMatrix **sqmatrices, TMatrix **matrices,
	    double *rhs, double *sol, 
	    MatVecProc *MatVect, DefectProc *Defect,
	    TMultiGrid2D *MG, 
	    int N_Unknowns, int ns_type)
{
  int solver_type, prec_type;
  double *itmethod_sol, *itmethod_rhs, t1, t2;
  TItMethod *itmethod, *prec;
  TSquareMatrix2D *A;

  t1 = GetTime();
  
  if (!ns_type)
  {
      prec_type = TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR;
      solver_type = TDatabase::ParamDB->SC_SOLVER_SCALAR;
  }
  else
  {
      prec_type = TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE;
      solver_type = TDatabase::ParamDB->SC_SOLVER_SADDLE;
  }

  switch(ns_type)
  {
      case 0:
	  A = (TSquareMatrix2D*) sqmatrices[0];
	  break;
  }
  /* if (TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE == 0)
    {
      TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE =  TDatabase::ParamDB->SOLVER_TYPE;
      TDatabase::ParamDB->SOLVER_TYPE = 2;
      }*/
  switch(TDatabase::ParamDB->SOLVER_TYPE)
  {
      case DIRECT:
	  switch(ns_type)
	  {
	      case 0:
		  DirectSolver(A, rhs, sol);
		  break;
	  }
	  //TDatabase::ParamDB->SOLVER_TYPE = 1;
	  break;
	  
      case AMG_SOLVE:
	  switch(ns_type)
	  {
	      case 0:
		  Solver(A, rhs, sol);
		  break;
	  }
	  break;
	  
      case GMG:
	  switch (prec_type)
	  {
	      case 1:
		  prec = new TJacobiIte(MatVect, Defect, NULL,
					0, N_Unknowns, 1);
		  break;
	      case 5:
		  prec = new TMultiGridScaIte(MatVect, Defect, NULL,
					      0, N_Unknowns, MG, 0);
		  break;
	      default:
		  OutPut("Unknown preconditioner !!!" << endl);
		  exit(4711);
	  }
	  switch (solver_type)
	  {
	      case 11:
		  itmethod = new TFixedPointIte(MatVect, Defect, prec,
						0, N_Unknowns, 1);
		  if (prec_type == 5)
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
		  itmethod = new TFgmresIte(MatVect, Defect, prec,
					    0, N_Unknowns, 1);
		  if (prec_type == 5)
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
	  
	  delete prec;
	  delete itmethod;
	  
	  switch (solver_type)
	  {
	      case 11:
		  if (prec_type == 5)
		  {
		      memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
		      memcpy(rhs, itmethod_rhs, N_Unknowns*SizeOfDouble);
		      delete itmethod_sol;
		      delete itmethod_rhs;
		  }
		  break;
	      case 16:
		  if (prec_type == 5)
		  {
		      memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
		      memcpy(rhs, itmethod_rhs, N_Unknowns*SizeOfDouble);
		      delete itmethod_sol;
		      delete itmethod_rhs;
		  }
		  break;
	  }
	  break;
  }
  t2 = GetTime();
  if (TDatabase::ParamDB->SC_VERBOSE>1)
    OutPut("time for solving: " << t2-t1 << endl);
}

/******************************************************************************/
// ComputeErrorEstimate
// computes residual based error estimates
// output in eta_K and estimated_global_error
/******************************************************************************/

void ComputeErrorEstimate(TCollection *coll, TFEFunction2D *u,
			  CoeffFct2D *Coeffs, BoundCondFunct2D **BoundaryConditions,
			  BoundValueFunct2D **BoundaryValues,
			  double* &eta_K,
			  double *maximal_local_error,
			  double *estimated_global_error,
			  double l2, double h1, double sd, int N_Unknowns)
{
    int N_Cells, bdry_tvd = 0, bdry_ad = 0;
    double t1, t2;
    TAuxParam2D *aux;
    TCD2DErrorEstimator *CDErrorEstimator;
    TFESpace2D *fesp[2];
    MultiIndex2D Derivatives_All[5] = { D10, D01, D00, D20, D02 };
 
    // set correct boundary conditions
    if (TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)
    {
	bdry_ad = TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM;
	TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM = 0;
    }
    if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD);
    {
	bdry_tvd = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;
	TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 0;
    }

    N_Cells = coll->GetN_Cells();
    // allocate arrays for local estimate
    eta_K = new double[N_Cells];

    CDErrorEstimator = new TCD2DErrorEstimator(TDatabase::ParamDB->ADAPTIVE_REFINEMENT_CRITERION,
					     u,
					     TDatabase::ParamDB->ERROR_CONTROL);
    t1 = GetTime();
    aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
    fesp[0] = u->GetFESpace2D();
    CDErrorEstimator->GetErrorEstimate(5, Derivatives_All,
				       Coeffs, BoundaryConditions,
				       BoundaryValues, aux, 1, fesp,
				       eta_K, &maximal_local_error[0], &estimated_global_error[0]);
    delete aux;
    delete CDErrorEstimator;
    t2 = GetTime();

    if (TDatabase::ParamDB->SC_VERBOSE>0)
      {
	OutPut("time for error estimation: " << t2-t1 << endl);
	
	OutPut("gradient estimate " << estimated_global_error[0] << endl);
	OutPut("L2 estimate " << estimated_global_error[2]);
	if ((l2>0)&&(TDatabase::ParamDB->MEASURE_ERRORS))
	  {
	    OutPut("   effectivity index " << estimated_global_error[2]/l2 << endl);
	  }
	else
	  {
	    OutPut(endl);
	  }
	OutPut("H1 estimate " << estimated_global_error[1]);
	if ((h1>0)&&(TDatabase::ParamDB->MEASURE_ERRORS))
	  {
	    OutPut("   effectivity index " << estimated_global_error[1]/h1 << endl);
	  }
	else
	  {
	    OutPut(endl);
	  }
	OutPut("energy estimate " << estimated_global_error[3]);
	if ((l2*l2+h1*h1>0)   &&(TDatabase::ParamDB->MEASURE_ERRORS))
	  {
	    OutPut("   effectivity index " << estimated_global_error[3]/
		   sqrt(l2*l2+h1*h1/TDatabase::ParamDB->RE_NR));// << endl);
	  }
	else
	  {
	      ;//OutPut(endl);
	  }
	OutPut(" (without edge res.) " << estimated_global_error[4]);
	OutPut(endl);
  
  if (TDatabase::ParamDB->DISCTYPE == SUPG)
  {
    OutPut(TDatabase::TimeDB->CURRENTTIME << " SUPG     estimate " << setw(12) << estimated_global_error[5]);
    if ((sd>0)&&(TDatabase::ParamDB->MEASURE_ERRORS))
    {
      OutPut(" effectivity index " << estimated_global_error[5]/sd << endl);
    }
    else
      OutPut(endl);
    //for (i=0;i<N_Cells;i++)
    //{
    //TDatabase::ParamDB->INTERNAL_P1_Array[i+N_Cells] = sqrt(eta_K[i]);
    //}
  }

  if (TDatabase::ParamDB->DISCTYPE == SDFEM)
  {
    OutPut(TDatabase::TimeDB->CURRENTTIME << " L2 rob   estimate " << setw(12) << estimated_global_error[6]);
    if ((sd>0)&&(TDatabase::ParamDB->MEASURE_ERRORS))
    {
      OutPut(" effectivity index " << estimated_global_error[6]/l2 << endl);
    }
    else
      OutPut(endl);
  }
  if (TDatabase::ParamDB->DISCTYPE == SDFEM)
  {
    OutPut("SD all " << sd << " " << estimated_global_error[5] << " " << estimated_global_error[6] << endl);
  }

	/*
	if (TDatabase::ParamDB->MEASURE_ERRORS)
	  OutPut("eff_ind " << N_Unknowns << " " << estimated_global_error[2]/l2
		 << " " << estimated_global_error[1]/h1 <<
		 " " << estimated_global_error[3]/
		 sqrt(l2*l2+h1*h1/TDatabase::ParamDB->RE_NR) << endl);
		 */
      }
    // reset parameters
      if (bdry_ad)
        TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM = bdry_ad;

      if (bdry_tvd)
        TDatabase::ParamDB->SOLD_PARAMETER_TYPE = bdry_tvd;
}

/******************************************************************************/
// ComputeAdjointMatrix
// computes the adjoint matrix of a matrix A
// it is assumed that both matrices have the same structure
/******************************************************************************/

void ComputeAdjointMatrix(TSquareMatrix2D *A,TSquareMatrix2D *AT)
{
    int i, j, k, n, m, ii, k1, k2, index, *row, *col;
    double *a, *at, val;

    m = A->GetN_Rows();
    n = A->GetN_Columns();
    row = A->GetRowPtr();
    col = A->GetKCol();
    a = A->GetEntries();
    at = AT->GetEntries();
  
    // loop over A
    j = row[0];
    for (i=0;i<m;i++)
    {
        ii = row[i+1];
	for (;j<ii;j++)
	{
	    index = col[j];
	    val = a[j];
	    // row in AT
	    k1 = row[index];
	    k2 = row[index+1];
	    for (k=k1;k<k2;k++)
	    {
		// compare column of AT with row of A
		if (col[k] == i)
		{
		    at[k] = val;
		    break;
		}
	    }
	}
    }
}

/******************************************************************************/
// OutputData2D
// writes the outputs
/******************************************************************************/
void OutputData2D(std::ostringstream& os, TOutput2D *Output,
		  int counter)
{
    if(TDatabase::ParamDB->WRITE_GRAPE)
    {
      os.seekp(std::ios::beg);
      os << TDatabase::ParamDB->GRAPEBASENAME << counter << ".dat" << ends;
      Output->WriteGrape(os.str().c_str());
    }
    if(TDatabase::ParamDB->WRITE_GMV)
    {
      os.seekp(std::ios::beg);
      os << TDatabase::ParamDB->GMVBASENAME  << counter << ".gmv" << ends;
      Output->WriteGMV(os.str().c_str());
    }

    if(TDatabase::ParamDB->WRITE_GNU)
    {
      os.seekp(std::ios::beg);
      os <<  TDatabase::ParamDB->GNUBASENAME  << counter << ".gnu" << ends;
      Output->WriteGnuplot(os.str().c_str());
    }

    if(TDatabase::ParamDB->WRITE_VTK)
    {
      os.seekp(std::ios::beg);
      os <<  TDatabase::ParamDB->VTKBASENAME << counter << ".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
    }
}

/******************************************************************************/
// PrepareAdjointProblem2D
// preparation of assembling and solution of adjoint problem
// the array for the piecewise constant stabilization parameters is filled
//    with standard values
/******************************************************************************/
void PrepareAdjointProblem2D(TCollection *coll,
			     TFEFunction2D* &u_adjoint,
			     TSquareMatrix2D* &sqmatrixAadjoint,
			     TFESpace2D* &pw_const_param_space,
			     TFEFunction2D* &pw_const_param_fe,
			     TFEFunction2D* &pw_const_param_deriv_fe,
			     TFESpace2D *velocity_space,
			     TSquareStructure2D *sqstructureA,
			     BoundCondFunct2D *BoundCondition,
			     CoeffFct2D *Coeffs,
			     double* &sol_adjoint,
			     double* &rhs_edge,
			     double* &pw_const_param,
			     double* &pw_const_param_deriv,
			     double* &pw_const_param_old,
			     int N_U, int N_Cells)
{
  int i, j, N_Edges, type;
  double *coeff, x, y, xs, ys;
  char UadjointString[] = "u_adjoint";
  char PwString[] = "pw_const_supg_param";
  char PwdString[] = "pw_const_supg_param_deriv";
  char Name[] = "name";
  char Description[] = "space for pw constant parameters";
  TBaseCell *cell;

  // coefficients of the equation
  coeff = new double[13];
  // allocate objects and arrays for solution and right hand side of adjoint problem
  sol_adjoint = new double[N_U];
  memset(sol_adjoint, 0, N_U*SizeOfDouble);
  u_adjoint = new TFEFunction2D(velocity_space, UadjointString, UadjointString, sol_adjoint, N_U);
  rhs_edge = new double[N_U];
  // allocate matrix for adjoint problem
  sqmatrixAadjoint = new TSquareMatrix2D(sqstructureA);
  // allocate objects and array that contains the piecewise constant stabilization parameters
  pw_const_param = new double[N_Cells];
  memset(pw_const_param, 0, N_Cells*SizeOfDouble);
  TDatabase::ParamDB->INTERNAL_P1_Array = pw_const_param;
  // piecewise constant finite element space
  pw_const_param_space = new TFESpace2D(coll, Name, Description, BoundCondition, 0, NULL);
  pw_const_param_fe = new TFEFunction2D(pw_const_param_space, PwString, PwString, pw_const_param, N_Cells);
  // allocate objects and array that contains the piecewise constant derivatives of the stabilization parameters
  pw_const_param_deriv = new double[N_Cells];
  pw_const_param_deriv_fe = new TFEFunction2D(pw_const_param_space, 
					      PwdString, PwdString, 
					      pw_const_param_deriv, N_Cells);
  // array for the previous stabilization parameters
  pw_const_param_old = new double[N_Cells];
  
  // initialization: fill pw_const_param_fe with the standard stabilization parameters
  type = TDatabase::ParamDB->SDFEM_TYPE;
  TDatabase::ParamDB->SDFEM_TYPE = 2;
  // loop over all cells
  for (i=0;i<N_Cells;i++)
  {
      cell = coll->GetCell(i);
      N_Edges=cell->GetN_Edges();
      // prepare data for computing mesh size in convection direction
      for (j=0;j<N_Edges;j++)
      {
	  TDatabase::ParamDB->INTERNAL_VERTEX_X[j] = cell->GetVertex(j)->GetX();
	  TDatabase::ParamDB->INTERNAL_VERTEX_Y[j] = cell->GetVertex(j)->GetY();
      }
      if (N_Edges==3)
	  TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
      TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;
      xs = ys = 0;
      for (j=0;j<N_Edges;j++)
      {
	  // compute coordinates
	  x = cell->GetVertex(j)->GetX();
	  y = cell->GetVertex(j)->GetY();
	  xs += x;
	  ys += y;
      }
      xs /= N_Edges;
      ys /= N_Edges;
      Coeffs(1, &xs, &ys, NULL, &coeff);
      pw_const_param[i] = 
	  Compute_SDFEM_delta(1, coeff[0], coeff[1], coeff[2], coeff[3], 0);
  }
  TDatabase::ParamDB->SDFEM_TYPE = type;
  delete coeff;
}

/******************************************************************************/
// SolveAdjointProblem2D
//  - assembles right hand side for adjoint problem
//  - solves adjoint problem
//  - computes derivatives of error estimator wrt stabilization parameters
/******************************************************************************/
void SolveAdjointProblem2D(TCollection *coll,
			   TDiscreteForm2D *DiscreteForm,
			   TFESpace2D *velocity_space,
			   TFEFunction2D *velo,
			   TFEFunction2D *u_adjoint,
			   TSquareMatrix2D* &sqmatrixAadjoint,
			   CoeffFct2D *Coeffs,
			   BoundCondFunct2D **BoundaryConditions,
			   BoundValueFunct2D **BoundaryValues,
			   double *rhs, double *rhs_edge,
			   double *sol_adjoint,
			   double *pw_const_param_deriv,
			   int N_U, int N_Active, int N_neum_to_diri,
			   int *neum_to_diri, int *neum_to_diri_bdry,
			   double *neum_to_diri_param)
{
    int solver_type, i;
    TAuxParam2D *aux;
    TSquareMatrix2D *SQMATRICES[1];
    TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
    MultiIndex2D FEMultiIndex_All_Deriv[5] = { D00, D10, D01, D20, D02 };
    MultiIndex2D FEMultiIndex_Sol[1] = { D00 };

    // assemble right hand side for adjoint problem
    // coercivity constant for rhs of adjoint problem
    TDatabase::ParamDB->INTERNAL_COERCIVITY = EstimateCoercivityConstant(coll, Coeffs);
    memset(rhs, 0, N_U*SizeOfDouble);
    // set aux object
    switch (TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)
    {
      case 100: 
	aux =  new TAuxParam2D(N_FESpaces_Sol, N_Fct_Sol, N_ParamFct_Sol,
			       N_FEValues_Sol, &velocity_space, &velo,
			       Fct_Sol,
			       FEFctIndex_Sol,  FEMultiIndex_Sol,
			       N_Params_Sol, BeginParam_Sol);
	break;
    default: 
      aux =  new TAuxParam2D(N_FESpaces_All_Deriv, N_Fct_All_Deriv, N_ParamFct_All_Deriv,
			     N_FEValues_All_Deriv, &velocity_space, &velo,
			     Fct_All_Deriv,
			     FEFctIndex_All_Deriv, FEMultiIndex_All_Deriv,
			     N_Params_All_Deriv, BeginParam_All_Deriv);
      break;
    }

    // this is just a dummy for the cell measure
    solver_type = TDatabase::ParamDB->CELL_MEASURE;
    // diameter
    TDatabase::ParamDB->CELL_MEASURE = 0;
    // assemble with homogeneous boundary conditions
    // IMPORTANT: bc. of primal problem to rule out Dirichlet boundaries, if set
    Assemble2D(1, &velocity_space,
	       0, NULL,
	       0, NULL,
	       1, &rhs, &velocity_space,
	       DiscreteForm,
	       BoundaryConditions,
	       BoundaryValues+2,
	       aux);
    // reset parameter
    TDatabase::ParamDB->CELL_MEASURE = solver_type;    
    // compute the contributions from the jumps across the edges
    // >= 100 -- adjoint problem with known error
    switch (TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)
      {
      case 1:
	memset(rhs_edge, 0, N_U*SizeOfDouble);
	JumpTermsForAdjointProblem(velocity_space, velo, Coeffs, BoundaryConditions[0], rhs_edge);
	// sum up both contributions
	Daxpy(N_U, 1, rhs, rhs_edge);
	break;
      default:
	memcpy(rhs_edge, rhs, N_U*SizeOfDouble);
	break;
      }

    // set Dirichlet boundary conditions
    SQMATRICES[0] = sqmatrixAadjoint;
    SetDirichletNodesFromNeumannNodes(SQMATRICES, rhs_edge, sol_adjoint,
				      N_neum_to_diri, neum_to_diri,
				      neum_to_diri_bdry, neum_to_diri_param,
				      BoundaryValues[2]);
    // this is just a dummy
    solver_type = TDatabase::ParamDB->SOLVER_TYPE;
    // change solver to direct solver
    TDatabase::ParamDB->SOLVER_TYPE = 2;
    Solver(sqmatrices, NULL,
	   rhs_edge, sol_adjoint,
	   NULL, NULL,
	   NULL, N_U, 0);
    // reset parameter
    TDatabase::ParamDB->SOLVER_TYPE = solver_type;

    // compute derivatives of error estimator wrt stabilization parameters
    ComputeDerivativeOfEstimator2D(coll, Coeffs, velo, u_adjoint, pw_const_param_deriv);
}

/******************************************************************************/
// ComputeDerivativeOfEstimator2D
/******************************************************************************/
void ComputeDerivativeOfEstimator2D(TCollection *coll,
				    CoeffFct2D *Coeffs,
				    TFEFunction2D *velo,
				    TFEFunction2D *u_adjoint,
				    double *pw_const_param_deriv)
{
  int i,j, N_Cells, N_Edges;
  double area, integral, xs, ys, x[4], y[4], val_velo[4], val_u_ad[4], *coeff;
  TBaseCell *cell;

  // coefficients of the equation
  coeff = new double[13];
  // number of mesh cells
  N_Cells = coll->GetN_Cells();
  // loop over all mesh cells
  for(i=0;i<N_Cells;i++)
  {
    // get cell no. i
    cell = coll->GetCell(i);
    // get number of edges
    N_Edges=cell->GetN_Edges();
    // get area of mesh cell
    area = cell->GetMeasure();
    // compute coordinates of vertices
    for (j=0;j<N_Edges;j++)
    {
      x[j] = cell->GetVertex(j)->GetX();
      y[j] = cell->GetVertex(j)->GetY();
    }
    // initialize integral
    integral = 0;
    // loop over the edges for edge mid point rule
    for (j=0;j<N_Edges;j++)
    {
	// barycenter of the edge
	xs = (x[j]+x[(j+1)%N_Edges])/2.0;
	ys = (y[j]+y[(j+1)%N_Edges])/2.0;
	// get value and derivatives of current solution of primal problem
	velo->FindGradientLocal(cell,i,xs,ys,val_velo);
	// get value and derivatives of solution of adjoint problem
	u_adjoint->FindGradientLocal(cell,i,xs,ys,val_u_ad);
	// get coefficients of the equation
	Coeffs(1, &xs, &ys, NULL, &coeff);
	// linearization of error estimator
	integral += (coeff[1]*val_velo[1] + coeff[2]*val_velo[2] + 
		     coeff[3]*val_velo[0] - coeff[4])*
	    (coeff[1]*val_u_ad[1] + coeff[2]*val_u_ad[2]);
    }
    // scaling for edge midpoint rule
    pw_const_param_deriv[i]  = -integral*area/3;
  }
  delete coeff;
}


void Assemble_CD_2D(TCollection *coll,
		    TFESpace2D **USpaces, TFEFunction2D **UArray,
		    double **RhsArray, TSquareMatrix2D **MatricesA,
		    CoeffFct2D *Coeffs,
		    BoundCondFunct2D **BoundaryConditions,
		    BoundValueFunct2D **BoundaryValues,
		    TDiscreteForm2D **DiscreteForms,
		    TSquareMatrix2D *sqmatrixAadjoint,
		    CheckWrongNeumannNodesFunct2D **CheckWrongNeumannNodesFct,
		    int *N_Uarray,
		    int low, int mg_level, int mg_type, int i,
		    int &N_neum_to_diri, int* &neum_to_diri,
		    int* &neum_to_diri_bdry, 
		    double* &neum_to_diri_param)
{
    int k, N_U, N_Active, N_NonActive, type;
    double *rhs, *RHSs[1];
    TFESpace2D *fesp[1], *ferhs[1];
    TSquareMatrix2D *SQMATRICES[1];
    TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
    TAuxParam2D *aux;
    TDiscreteForm2D *DiscreteForm;
    BoundCondFunct2D *BdCond;
    BoundValueFunct2D *BdVal;

    for(k=low;k<=mg_level;k++)
    {
      rhs = RhsArray[k];
      N_U = N_Uarray[k];
      N_Active = USpaces[k]->GetActiveBound();
      N_NonActive = N_U - N_Active;

      RHSs[0] = rhs;
      memset(rhs, 0, N_U*SizeOfDouble);

      // find discrete form
      if ((mg_type==1) && (k<i+1))
      {
        DiscreteForm = DiscreteForms[2];
	OutPut("Upwind ");
      }
      else
        switch(TDatabase::ParamDB->DISCTYPE)
        {
          case GALERKIN:
          case LOCAL_PROJECTION:
	      DiscreteForm = DiscreteForms[0];
	      OutPut("Galerkin ");
	      break;

          case SDFEM:
            DiscreteForm = DiscreteForms[1];
	    type = TDatabase::ParamDB->SDFEM_TYPE;
	    if ((TDatabase::ParamDB->SDFEM_TYPE == 100)&&(k<mg_level))
	    {
		TDatabase::ParamDB->SDFEM_TYPE = 2;
	    }
            break;

          case UPWIND:
	      DiscreteForm = DiscreteForms[2];
	      OutPut("Upwind ");
	      break;

        default:
          OutPut("Unknown DISCTYPE" << endl);
          exit(4711);;
      }

      fesp[0] = USpaces[k];
      ferhs[0] = USpaces[k];
      // initialize matrix
      SQMATRICES[0] = MatricesA[k];
      SQMATRICES[0]->Reset();

      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      SetPolynomialDegree();
      BdCond = BoundaryConditions[0];
      BdVal = BoundaryValues[0];
      // set homogeneous Neumann boundary conditions
      if ((TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD)
	  || (TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM))
      {
	  BdCond = BoundaryConditions[1];
	  BdVal = BoundaryValues[1];
      }
      // assemble
      Assemble2D(1, fesp,
		 1, SQMATRICES,
		 0, NULL,
		 1, RHSs, ferhs,
		 DiscreteForm,
		 &BdCond,
		 &BdVal,
		 aux);
      if (TDatabase::ParamDB->SC_VERBOSE>1)
	OutPut("Assembling done on mg-level: "<< mg_level<<endl);

      // DG and CIP
      switch (TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS)
	{
	case 1:  Assemble2D_CIP(Coeffs,1, fesp,
				1, SQMATRICES,
				0, NULL,
				1, RHSs, ferhs,
				BoundaryConditions,
				BoundaryValues,
				aux);
	  break;
	case 2:  Assemble2D_DG(Coeffs,1, fesp,
			       1, SQMATRICES,
			       0, NULL,
			       1, RHSs, ferhs,
			       BoundaryConditions,
			       BoundaryValues,
			       aux);
	  break;
	}

      // reset SUPG parameter
      if (TDatabase::ParamDB->DISCTYPE == SDFEM)
      {
	  TDatabase::ParamDB->SDFEM_TYPE = type;
      }	  
      // LPS ???
      if(TDatabase::ParamDB->DISCTYPE == LOCAL_PROJECTION)
      {
	  if(TDatabase::ParamDB->LP_FULL_GRADIENT)
	      ;//UltraLocalProjection(MatricesA[k], FALSE, Coefficients[0]);
	  
	  if(TDatabase::ParamDB->LP_STREAMLINE)
	  {
	      OutPut("local projection stabilisation in streamline direction ");
	      OutPut("is currently not available." << endl);
	      exit(4711);
	  }
      }
      // apply upwind
      if (DiscreteForm == DiscreteForms[2])
      {
	  UpwindForConvDiff(Coeffs, SQMATRICES[0],RHSs[0],
			    fesp[0],DiscreteForms[2],NULL,NULL,0);
	  cout << "UPWINDING DONE : level " << k << endl;
      }                                           // endif
      delete aux;

      // compute adjoint matrix, only on the fines level
      if ((TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)&&(k==mg_level))
      {
	  ComputeAdjointMatrix(MatricesA[k],sqmatrixAadjoint);
      }
      if ((TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD)
	  || (TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM))
      {
	  CheckWrongNeumannNodesFct[0](coll, USpaces[k], N_neum_to_diri, neum_to_diri,
				       neum_to_diri_bdry, neum_to_diri_param);
	  if (N_neum_to_diri)
	      SetDirichletNodesFromNeumannNodes(SQMATRICES, RHSs[0], UArray[k]->GetValues(),
						N_neum_to_diri, neum_to_diri,
						neum_to_diri_bdry, neum_to_diri_param,
						BoundaryValues[0]);
      }
    }	  
}

/******************************************************************************/
// restrict parameters to admitted values 
/******************************************************************************/
void  RestrictParametersToAdmissibleValues(TCollection *coll,
					   TFESpace2D *pw_const_param_space,
                                           CoeffFct2D *Coeffs,
					   double* pw_const_param)
{
    int i, j, N_Cells, N_Edges, index, order;
    int *GlobalNumbers, *BeginIndex, *DOF;
  double omega = TDatabase::ParamDB->INTERNAL_COERCIVITY, val1, val2, x, y;
  double *coeff, c_inv2 = 1, hK, lower_bound;
  TBaseCell *cell;
  FE2D CurrentElement;

  N_Cells = coll->GetN_Cells();
  GlobalNumbers = pw_const_param_space->GetGlobalNumbers();
  BeginIndex = pw_const_param_space->GetBeginIndex();
  order = TDatabase::ParamDB->ANSATZ_ORDER;

  coeff = new double[13];
  // loop over the cells
  for (i=0;i<N_Cells;i++)
  {
      cell = coll->GetCell(i);
      hK = cell->GetDiameter();
      lower_bound = 1e-12 * hK * hK;
      // dof of current mesh cell
      DOF = GlobalNumbers + BeginIndex[i];
      index = DOF[0];
      // apply lower bound
      if (pw_const_param[index]<=lower_bound)
      {
	  pw_const_param[index] = lower_bound;
	  continue;
      }
      // compute maximal value
      // compute first barycenter
      N_Edges = cell->GetN_Edges();
      x = y = 0;
      for (j=0; j<N_Edges; j++)
      {
	  x += cell->GetVertex(j)->GetX();
	  y += cell->GetVertex(j)->GetY();
      }
      x /= N_Edges;
      y /= N_Edges;
      // get coefficients
      Coeffs(1, &x, &y, NULL, &coeff);
      // finite element on the mesh cell
      CurrentElement = pw_const_param_space->GetFE2D(i, cell);
      switch(CurrentElement)
      {
	  // P_1, Q_1
	  case C_P1_2D_T_A:
	  case C_Q1_2D_Q_A:
	  case C_Q1_2D_Q_M:
	      c_inv2 = 1.0;
	      break;
	case C_P2_2D_T_A:
	    c_inv2 = 48.0;
	    break;
	case C_Q2_2D_Q_A:
	case C_Q2_2D_Q_M:
	    c_inv2 = 24.0;
	    break;
	case C_P3_2D_T_A:
	    c_inv2 = (435+sqrt(26025.0))/4.0;
	    break;
	case C_Q3_2D_Q_A:
	case C_Q3_2D_Q_M:
	    c_inv2 = (244+sqrt(9136.0))/3.0;
	    break;
	default:
	    c_inv2 = (435+sqrt(26025.0))/4.0;
	    break;
      }
	    
      // first parameter
      val1 = hK * hK;
      val1 /= (coeff[0]*c_inv2);
      // this bound not necessary for linear and bilinear fe (on parallelograms)
      if ((CurrentElement==C_P1_2D_T_A) || (CurrentElement==C_Q1_2D_Q_A) || (CurrentElement==C_Q1_2D_Q_M))
	  val1 = 1e36;
      // there is a second bound for omega > 0
      if (omega > 0)
      {
	  // second parameter
	  if (coeff[3]!=0)
	      val2 = omega/(coeff[3]*coeff[3]);
	  else
	      val2 = 1e36;
	  // put minimum on val1
	  if (val2 < val1)
	      val1 = val2/2.0;
	  else
	      val1 /= 2.0;
      }
      // apply upper bound
      if (pw_const_param[index] > val1)
      {
	  pw_const_param[index] = val1;
      }
  }
  delete coeff;
}

/******************************************************************************/
// computes the wrong Neumann nodes for problems on the unit square with
// Dirichlet boundary conditions
/******************************************************************************/

void CheckWrongNeumannNodesUnitSquareDiri(TCollection *Coll, TFESpace2D *fespace,
			    int &N_neum_to_diri, int* &neum_to_diri,
			    int* &neum_to_diri_bdry, 
			    double* &neum_to_diri_param)
{
    const int max_entries = 4096;  
    int i, j, N_, min_val, type;
  int N_Cells, N_V, diri_counter = 0, found, diri_counter_1 = 0;
  int *global_numbers, *begin_index, *dof;
  int boundary_vertices[4], tmp_diri[max_entries], tmp_bdry[max_entries];
  double x[4], y[4], eps = 1e-6, tmp_param[max_entries];
  TBaseCell *cell;
  TVertex *vertex;
  FE2D CurrentElement;

  // number of mesh cells
  N_Cells = Coll->GetN_Cells();
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();

  diri_counter = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    found = 0;
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      boundary_vertices[j] = 0;
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j]);
      // vertex on the upper lid
      if ((fabs(x[j])<eps)||(fabs(y[j])<eps)||(fabs(1-x[j])<eps)||(fabs(1-y[j])<eps))
      {
	   boundary_vertices[j] = 1;
	   found++;
      }
    }
    // no cell with edge with vertex on the boundary
    if (found<2) 
	continue;
    // finite element on the mesh cell
    CurrentElement = fespace->GetFE2D(i, cell);
    // number of basis functions (= number of d.o.f.)
    N_ = TFEDatabase2D::GetN_BaseFunctFromFE2D(CurrentElement);
    // the array which gives the mapping of the local to the global d.o.f.
    dof = global_numbers+begin_index[i];
    switch(CurrentElement)
    {
	// P_1, Q_1
	case C_P1_2D_T_A:
	case C_Q1_2D_Q_A:
	case C_Q1_2D_Q_M:
	    for (j=0;j<N_V;j++)
	    {
		// vertex on the boundary
		if (boundary_vertices[j])
		{
		    if (CurrentElement==C_P1_2D_T_A)
			tmp_diri[diri_counter] = dof[j];
		    else
		    {
			if (j<2)
			    tmp_diri[diri_counter] = dof[j];
			else
			{
			    if (j==2)
				tmp_diri[diri_counter] = dof[3];
			    else
				tmp_diri[diri_counter] = dof[2];
			}
		    }
		    if (diri_counter > max_entries)
		    {
			OutPut("tmp_diri too short !!!"<<endl);
			exit(4711);
		    }
		    if (fabs(y[j])<eps)
		    {
			tmp_bdry[diri_counter] = 0;
			tmp_param[diri_counter] = x[j];
		    }
		    if (fabs(1-y[j])<eps)
		    {
			tmp_bdry[diri_counter] = 2;
			tmp_param[diri_counter] = 1-x[j];
		    }
		    if (fabs(x[j])<eps)
		    {
			tmp_bdry[diri_counter] = 3;
			tmp_param[diri_counter] = 1-y[j];
		    }
		    if (fabs(1-x[j])<eps)
		    {
			tmp_bdry[diri_counter] = 1;
			tmp_param[diri_counter] = y[j];
		    }
		    diri_counter++;
		}
	    }
	    break;
	// P_2, Q_2
	case C_P2_2D_T_A:
	case C_Q2_2D_Q_A:
	case C_Q2_2D_Q_M:
            // loop over the edges
 	    for (j=0;j<N_V;j++)
	    {
              // check of edge j is on boundary  
              if (boundary_vertices[j] && boundary_vertices[(j+1)%N_V])
              {
		// check if this is a boundary edge
		type = cell->GetJoint(j)->GetType();
		if (!((type == BoundaryEdge)||(type == IsoBoundEdge)))
		  continue;
	        switch(j)
                {
                   case 0:
                     tmp_diri[diri_counter] = dof[0];
                     tmp_diri[diri_counter+1] = dof[1];
                     tmp_diri[diri_counter+2] = dof[2];
                   break;
                  case 1:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[2];
                       tmp_diri[diri_counter+1] = dof[4];
                       tmp_diri[diri_counter+2] = dof[5];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[2];
                       tmp_diri[diri_counter+1] = dof[5];
                       tmp_diri[diri_counter+2] = dof[8];
                     }
                   break;
                  case 2:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[5];
                       tmp_diri[diri_counter+1] = dof[3];
                       tmp_diri[diri_counter+2] = dof[0];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[8];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[6];
                     }
                   break;
                   case 3:
                     tmp_diri[diri_counter] = dof[6];
                     tmp_diri[diri_counter+1] = dof[3];
                     tmp_diri[diri_counter+2] = dof[0];
                   break;

                }
              
		if (diri_counter+2 > max_entries)
		{
			OutPut("tmp_diri too short !!!"<<endl);
			exit(4711);
		}
		// boundary at y=0
		if ((fabs(y[j])<eps)&&(fabs(y[(j+1)%N_V])<eps))
		  {
		    tmp_bdry[diri_counter] = 0;
		    tmp_bdry[diri_counter+1] = 0;
		    tmp_bdry[diri_counter+2] = 0;
		    tmp_param[diri_counter] = x[j];
		    tmp_param[diri_counter+1] = (x[j] + x[(j+1)%N_V])/2.0;
		    tmp_param[diri_counter+2] = x[(j+1)%N_V];
		  }
		// boundary at y = 1
		if ((fabs(1-y[j])<eps)&&(fabs(1-y[(j+1)%N_V])<eps))
		  {
		    tmp_bdry[diri_counter] = 2;
		    tmp_bdry[diri_counter+1] = 2;
		    tmp_bdry[diri_counter+2] = 2;
		    tmp_param[diri_counter] = 1-x[j];
		    tmp_param[diri_counter+1] = (1-x[j] + 1-x[(j+1)%N_V])/2.0;
		    tmp_param[diri_counter+2] = 1-x[(j+1)%N_V];
		  }
		// boundary at x = 0
		if ((fabs(x[j])<eps)&&(fabs(x[(j+1)%N_V])<eps))
		  {
		    tmp_bdry[diri_counter] = 3;
		    tmp_bdry[diri_counter+1] = 3;
		    tmp_bdry[diri_counter+2] = 3;
		    tmp_param[diri_counter] = 1-y[j];
		    tmp_param[diri_counter+1] = (1-y[j] + 1-y[(j+1)%N_V])/2.0;
		    tmp_param[diri_counter+2] = 1-y[(j+1)%N_V];
		  }
		// boundary at x = 1
		if ((fabs(1-x[j])<eps)&&(fabs(1-x[(j+1)%N_V])<eps))
		  {
		    tmp_bdry[diri_counter] = 1;
		    tmp_bdry[diri_counter+1] = 1;
		    tmp_bdry[diri_counter+2] = 1;
		    tmp_param[diri_counter] = y[j];
		    tmp_param[diri_counter+1] = (y[j] + y[(j+1)%N_V])/2.0;
		    tmp_param[diri_counter+2] = y[(j+1)%N_V];
		  }
		diri_counter +=3;
	      }
	    }
	    break;
	// P_3, Q_3
	case C_P3_2D_T_A:
	case C_Q3_2D_Q_A:
	case C_Q3_2D_Q_M:
            // loop over the edges
 	    for (j=0;j<N_V;j++)
	    {
              // check of edge j is on boundary  
              if (boundary_vertices[j] && boundary_vertices[(j+1)%N_V])
              {
		// check if this is a boundary edge
		type = cell->GetJoint(j)->GetType();
		if (!((type == BoundaryEdge)||(type == IsoBoundEdge)))
		  continue;

               // P3: local dof 0, 1, 2, 3 are on the boundary
               // Q3: local dof 0, 1, 2, 3 are on the boundary
	        switch(j)
                {
                   case 0:
                     tmp_diri[diri_counter] = dof[0];
                     tmp_diri[diri_counter+1] = dof[1];
                     tmp_diri[diri_counter+2] = dof[2];
		     tmp_diri[diri_counter+3] = dof[3];
                   break;
                  case 1:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[3];
                       tmp_diri[diri_counter+1] = dof[6];
                       tmp_diri[diri_counter+2] = dof[8];
		       tmp_diri[diri_counter+3] = dof[9];
                     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[3];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[11];
		       tmp_diri[diri_counter+3] = dof[15];
                     }
                   break;
                  case 2:
                     if (N_V==3)
                     {
                       tmp_diri[diri_counter] = dof[9];
                       tmp_diri[diri_counter+1] = dof[7];
                       tmp_diri[diri_counter+2] = dof[4];
                       tmp_diri[diri_counter+3] = dof[0];
		     }
                     else
                     {
                       tmp_diri[diri_counter] = dof[15];
                       tmp_diri[diri_counter+1] = dof[14];
                       tmp_diri[diri_counter+2] = dof[13];
			tmp_diri[diri_counter+3] = dof[12];
                     }
                   break;
                   case 3:
                     tmp_diri[diri_counter] = dof[12];
                     tmp_diri[diri_counter+1] = dof[8];
                     tmp_diri[diri_counter+2] = dof[4];
		     tmp_diri[diri_counter+3] = dof[0];
                   break;
                }
              
		if (diri_counter+3 > max_entries)
		{
			OutPut("tmp_diri too short !!!"<<endl);
			exit(4711);
		}

		if ((fabs(y[j])<eps)&&(fabs(y[(j+1)%N_V])<eps))
		{
		    tmp_bdry[diri_counter] = 0;
		    tmp_bdry[diri_counter+1] = 0;
		    tmp_bdry[diri_counter+2] = 0;
		    tmp_bdry[diri_counter+3] = 0;
		    tmp_param[diri_counter] = x[j];
		    tmp_param[diri_counter+1] = 2*x[j]/3.0 + x[(j+1)%N_V]/3.0;
		    tmp_param[diri_counter+2] = x[j]/3.0 + 2*x[(j+1)%N_V]/3.0; 
		    tmp_param[diri_counter+3]= x[(j+1)%N_V];
		}
		if ((fabs(1-y[j])<eps)&&(fabs(1-y[(j+1)%N_V])<eps))
		{
		    tmp_bdry[diri_counter] = 2;
		    tmp_bdry[diri_counter+1] = 2;
		    tmp_bdry[diri_counter+2] = 2;
		    tmp_bdry[diri_counter+3] = 2;
		    tmp_param[diri_counter] = 1-x[j];
		    tmp_param[diri_counter+1] = 2*(1-x[j])/3.0 + (1-x[(j+1)%N_V])/3.0;
		    tmp_param[diri_counter+2] = (1-x[j])/3.0 + 2*(1-x[(j+1)%N_V])/3.0;
		    tmp_param[diri_counter+3] = 1-x[(j+1)%N_V];
		}
		    if ((fabs(x[j])<eps)&&(fabs(x[(j+1)%N_V])<eps))
		    {
			tmp_bdry[diri_counter] = 3;
			tmp_bdry[diri_counter+1] = 3;
			tmp_bdry[diri_counter+2] = 3;
			tmp_bdry[diri_counter+3] = 3;
			tmp_param[diri_counter] = 1-y[j];
			tmp_param[diri_counter+1] = 2*(1-y[j])/3.0 + (1-y[(j+1)%N_V])/3.0;
			tmp_param[diri_counter+2] = (1-y[j])/3.0 + (1-y[(j+1)%N_V])/3.0;
			tmp_param[diri_counter+3] = 1-y[(j+1)%N_V];
		    }
		    if ((fabs(1-x[j])<eps)&&(fabs(1-x[(j+1)%N_V])<eps))
		    {
			tmp_bdry[diri_counter] = 1;
			tmp_bdry[diri_counter+1] = 1;
			tmp_bdry[diri_counter+2] = 1;
			tmp_bdry[diri_counter+3] = 1;
			tmp_param[diri_counter] = y[j];
			tmp_param[diri_counter+1] = 2*y[j]/3.0 + y[(j+1)%N_V]/3.0;
			tmp_param[diri_counter+2] = y[j]/3.0 + 2*y[(j+1)%N_V]/3.0;
			tmp_param[diri_counter+3] = y[(j+1)%N_V];
		    }
		    diri_counter +=4;
		}
	    }
	    break;
	default:
	    OutPut("CheckNeumannNodesForVelocity not implemented for element "
		   << CurrentElement << endl);
	    OutPut("code can be run without CheckNeumannNodesForVelocity, just delete the exit" << endl);
	    exit(4711);
    }	    
  }
 
  // condense
  for (i=0;i<diri_counter;i++)
  {
      if (tmp_diri[i] == -1)
	  continue;
      diri_counter_1++;
      for (j=i+1;j<diri_counter;j++)
      {
	  if (tmp_diri[i] == tmp_diri[j])
	  {
	      tmp_diri[j] = -1;
	  }
      }
  }

  //OutPut("CheckNeumannNodesForVelocity: N_neum_to_diri " << diri_counter_1 << endl);
  N_neum_to_diri = diri_counter_1;
  // allocate array for the indices
  neum_to_diri = new int[diri_counter_1];
  // allocate array for the corresponding boundary numbers
  neum_to_diri_bdry = new int[diri_counter_1];
  // allocate array for the corresponding boundary parameters
  neum_to_diri_param = new double[diri_counter_1];
  // fill array and sort
  for (i=0;i<diri_counter_1;i++)
  {
      min_val = tmp_diri[0];
      found = 0;
      for (j=1;j<diri_counter;j++)
      {
	  if ((tmp_diri[j]>0) && ((tmp_diri[j] < min_val) || 
				  (min_val == -1)))
	  {
	       min_val =  tmp_diri[j];
	       found = j;
	  }
      }
      neum_to_diri[i] = tmp_diri[found];
      neum_to_diri_bdry[i] = tmp_bdry[found];
      neum_to_diri_param[i] = tmp_param[found];
      tmp_diri[found] = -1;
  }
}


/******************************************************************************/
// AllocateArraysCD_2D
// allocates arrays for main program CDAdap_2D.C
/******************************************************************************/
void AllocateArraysCD_2D(int LEVELS, TFEFunction2D **&UArray,
TFEVectFunct2D  **&VelocityArray,
TFESpace2D **&USpaces,
TSquareMatrix2D **&MatricesA, TMultiGrid2D *&MG,
double **&RhsArray,
double *&l2, double *&h1, double *&energy,
double *&sd, double *&l_inf, double *&lp,
double *Parameters,
int *&N_Uarray,
int &mg_type, int &mg_level)
{
  int i;

  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SCALAR;
  if (TDatabase::ParamDB->SOLVER_TYPE==AMG_SOLVE||
    TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
    mg_type = 0;
  if (mg_type)
    mg_level = 0;
  else
    mg_level = -1;

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

  // initialize multigrid
  if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
    Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
    MG = new TMultiGrid2D(1, 2, Parameters);
  }
  // allocate array for velocity field
  VelocityArray = new TFEVectFunct2D*[LEVELS+1];
  for (i=0;i<LEVELS+1;i++)
    VelocityArray[i] = NULL;
}


/******************************************************************************/
// RefineGridCD_2D
// refines the grid for main program CDAdap_2D.C
/******************************************************************************/
void RefineGridCD_2D(TDomain *Domain, TCollection *coll,
int LEVELS, int curr_level,
double *&eta_K, double &maximal_local_error,
double *estimated_global_error,
int current_estimator)
{
  int BASELEVEL = TDatabase::ParamDB->UNIFORM_STEPS;

  // refine grid if level is greater than 0
  if ((curr_level)&&(curr_level<LEVELS))
  {
    // regular refinement if
    // adaptive procedure should not yet start
    // or no error estimation
    if ((curr_level<= BASELEVEL)||(! TDatabase::ParamDB->ESTIMATE_ERRORS))
    {
      Domain->RegRefineAll();
    }
    else if (curr_level < LEVELS)
    {
      // conforming closure
      if (TDatabase::ParamDB->GRID_TYPE)
        Domain->RefineByErrorEstimator(coll,eta_K,maximal_local_error,
          estimated_global_error[current_estimator],
          TRUE);
      else
        // hanging nodes
        Domain->RefineByErrorEstimator(coll,eta_K,maximal_local_error,
          estimated_global_error[current_estimator],
          FALSE);
      delete eta_K;
    }
  }
}


/******************************************************************************/
// GetCollectionCD_2D()
// gets the current collection
/******************************************************************************/
void GetCollectionCD_2D(TDomain *Domain, TCollection *&coll,
int i, int level, int mg_level, int &N_Cells)
{
  std::ostringstream os;

  coll=Domain->GetCollection(It_LE, level+TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR);
  if (level==mg_level)
  {
    N_Cells = coll->GetN_Cells();
    OutPut( "level " << level << " number of cells: " << N_Cells << endl);
  }

  // for adaptive methods, all levels are created new
  // but picture only on the finest level
  if((level==mg_level)&&(TDatabase::ParamDB->WRITE_PS))
  {
    // write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << TDatabase::ParamDB->PSBASENAME << i << ".grid.ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
  }
}

/******************************************************************************/
// PrepareDataOutputCD_2D()
// add (vector) functions to Output
/******************************************************************************/
void PrepareDataOutputCD_2D(
TOutput2D *Output,
TDomain *Domain,
TFEFunction2D *u,
TFEVectFunct2D *pw_const_proj_fefct)
{

  if (!TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)
  {
    /*  N_Parameters = 4;
        // allocate objects and array that contains the piecewise constant stabilization parameters
        pw_const_param = new double[N_Parameters * N_Cells];
        memset(pw_const_param, 0, N_Parameters * N_Cells*SizeOfDouble);
        TDatabase::ParamDB->INTERNAL_P1_Array = pw_const_param;
        // piecewise constant finite element space for parameters
        pw_const_param_space = new TFESpace2D(coll, Name, Description, BoundCondition, 0, NULL);
        // the corresponding finite element function
        pw_const_param_fe = new TFEVectFunct2D(pw_const_param_space, ErrString, ErrString, pw_const_param, N_Cells, 2);
        pw_const_param_deriv_fe = new TFEVectFunct2D(pw_const_param_space, SoldString, SoldString, pw_const_param+2*N_Cells, N_Cells, 2);

    Output = new TOutput2D(2, 2, 2, 1,Domain);
    Output->AddFEFunction(u);
    Output->AddFEVectFunct(pw_const_param_fe);
    Output->AddFEVectFunct(pw_const_param_deriv_fe);
    */
    //Output = new TOutput2D(1, 1, 0, 0,Domain);
    //Output->AddFEFunction(u);
    Output = new TOutput2D(2, 3, 0, 0,Domain);
    Output->AddFEFunction(u);
    if (pw_const_proj_fefct != NULL)
    {
      Output->AddFEFunction(pw_const_proj_fefct->GetComponent(0));
      Output->AddFEFunction(pw_const_proj_fefct->GetComponent(1));
    }
  }
  // code for 'else' is now in OptimizeSolutionAdjoint
}

/******************************************************************************/
// allocate spaces and arrays for CDAdap_2D.C
/****************F**************************************************************/
void AllocateSpacesArraysCD_2D(TCollection *coll, TFESpace2D **&USpaces,
TFESpace2D *&sold_space, TFESpace2D *&pw_const_proj_space,
TFEFunction2D **sc_params_fe, TFEVectFunct2D *&pw_const_proj_fefct,
TFEFunction2D **UArray,
BoundCondFunct2D **BoundaryConditions,
BoundValueFunct2D **BoundaryValues,
int level, int mg_type, int mg_level, int FirstSolve, int &N_Active,
int *N_Uarray,
double *&sol, double *&oldsol, double *&update,
double *&rhs,
double *&rhs_edge, double *&sc_params, double *&pw_const_proj,
double **RhsArray)
{
  TFESpace2D *fe_space;
  TFEFunction2D *u;
  TFEVectFunct2D *sc_params_vect;
  int bd, N_U, sc_N_params, N_Cells, N_sc,m;
  int sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;
  double hmin, hmax;
  char Description[] = "description";
  char Name[] = "name";
  char SoldString[] = "sold";
  char UString[] = "concentration";

  // fe space without boundary conditions
  if ((TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD)
    || (TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM))
    bd = 1;
  else
    bd = 0;

  // get spaces for low order disc on finest geometry grid
  if ((mg_type==1)&&(level<mg_level))
  {
    fe_space =  new TFESpace2D(coll, Name, Description,
      BoundaryConditions[bd], 1, NULL);
  }
  // get spaces of high order disc on grid
  else
  {
    fe_space =  new TFESpace2D(coll, Name, Description,
      BoundaryConditions[bd], TDatabase::ParamDB->ANSATZ_ORDER, NULL);
  }

  // build fespace hierarchy
  // set values and pointers for fe space
  USpaces[level] = fe_space;
  N_U = fe_space->GetN_DegreesOfFreedom();
  N_Uarray[level] = N_U;
  N_Active = fe_space->GetActiveBound();

  // SOLD schemes
  if (TDatabase::ParamDB->SOLD_TYPE)
  {
    // allocate piecewise constant finite element space
    sold_space = new TFESpace2D(coll, Name, Description, BoundaryConditions[0], 0, NULL);
    sc_N_params = 2;
    N_sc = sold_space->GetN_DegreesOfFreedom();
    if (sc_N_params)
    {
      sc_params = new double[sc_N_params* N_sc];
      sc_params_vect = new TFEVectFunct2D(sold_space, SoldString, SoldString, sc_params,
        N_sc, sc_N_params);
      for (m=0;m<sc_N_params;m++)
        sc_params_fe[m] = sc_params_vect->GetComponent(m);
    }
    if (((sold_parameter_type >= BH04)&& (sold_parameter_type <= BE05_2)) ||
      (sold_parameter_type == MH_Kno06) ||
      (sold_parameter_type ==  FEM_TVD))
    {
      rhs_edge = new double[N_U];
      memset(rhs_edge, 0, N_U*SizeOfDouble);
    }
  }

  // local projection scheme with two levels
  if (TDatabase::ParamDB->DISCTYPE == LOCAL_PROJECTION_2_LEVEL)
  {
    pw_const_proj_space = new TFESpace2D(coll, SoldString, SoldString, BoundaryConditions[0], 0, NULL);
    pw_const_proj = new double[2*N_Cells];
    memset(pw_const_proj, 0, 2*N_Cells*SizeOfDouble);
    pw_const_proj_fefct = new TFEVectFunct2D(pw_const_proj_space, SoldString, SoldString, pw_const_proj, N_Cells, 2);
  }
  else
    pw_const_proj_fefct = NULL;

  if (((level>=FirstSolve)||(mg_type==0))&&(level==mg_level))
    OutPut("dof scalar   : "<< setw(10) << N_U  << endl);

  if ((mg_type==1)&&(level==mg_level-1))
  {
    OutPut("dof low order disc     : "<<  setw(10) << N_U  << endl);
  }
  // for a posteriori paper with Julia Novo
  if (N_U>TDatabase::ParamDB->P2)
  {
    OutPut("Has to be included for a posteriori error estimator " <<endl);
    //TDatabase::ParamDB->UNIFORM_STEPS = level;
    //OutPut("Change parameter TDatabase::ParamDB->UNIFORM_STEPS = "<< TDatabase::ParamDB->UNIFORM_STEPS << endl);
    //TDatabase::ParamDB->P2 = 1e36;
  }
  // this is just to stop the simulations for the a posteriori error estimators
  if (N_U>TDatabase::ParamDB->P3)
  {
    OutPut("Has to be included for a posteriori error estimator " <<endl);
    //TDatabase::ParamDB->LEVELS=level+1;
  }

  coll->GetHminHmax(&hmin,&hmax);
  OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
  TDatabase::ParamDB->P4 = hmax;
  // initialize arrays
  rhs = new double[N_U];
  memset(rhs, 0, N_U*SizeOfDouble);
  RhsArray[level] = rhs;
  sol = new double[N_U];
  oldsol = new double[N_U];
  update = new double[N_U];
  memset(sol, 0, N_U*SizeOfDouble);
  memset(oldsol, 0, N_U*SizeOfDouble);
  memset(update, 0, N_U*SizeOfDouble);

  // build new fe functions
  u = new TFEFunction2D(fe_space, UString, UString, sol, N_U);
  UArray[level] = u;
  // read velocity field, this is done directly in the example file
}


/******************************************************************************/
// allocate spaces and arrays for CDAdap_2D.C
/****************F**************************************************************/
void AllocateMatricesCD_2D(TFESpace2D *fe_space,TSquareStructure2D *&sqstructureA,
TSquareMatrix2D *&sqmatrixA, TSquareMatrix2D **MatricesA,
int k, double *&matrix_D_Entries)
{
  // build matrices
  sqstructureA = new TSquareStructure2D(fe_space);
  sqstructureA->Sort();
  // allocate matrices
  sqmatrixA = new TSquareMatrix2D(sqstructureA);
  MatricesA[k] = sqmatrixA;

  if (TDatabase::ParamDB->SOLD_TYPE)
  {
    if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE ==  FEM_TVD)
    {
      matrix_D_Entries = new double[sqmatrixA->GetN_Entries()];
      memset(matrix_D_Entries, 0, sqmatrixA->GetN_Entries()*SizeOfDouble);
    }
  }
}


/******************************************************************************/
// prepares the solver for CDAdap_2D.C
/****************F**************************************************************/
void PrepareSolverCD_2D(TSquareMatrix2D *sqmatrixA, TMultiGrid2D *MG, TMGLevel2D *&MGLevel,
int level, int mg_type, int mg_level, int &low,
double *sol, double *rhs)
{
  int n_aux;

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
      MGLevel = new TMGLevel2D(level, sqmatrixA,
        rhs, sol, n_aux, NULL);
      if ((level==mg_level)||((level==0)&&(mg_level==1)&&(mg_type==1)))
        MG->AddLevel(MGLevel);
      else
        MG->ReplaceLevel(level,MGLevel);
      break;
  }
}


/******************************************************************************/
// prepares the output for CDAdap_2D.C
/****************F**************************************************************/
void PrepareOutPutCD_2D(TDomain *Domain, TFEFunction2D *conc, TFEVectFunct2D *velocity,
TFEFunction2D *conc_adjoint,
TFEVectFunct2D *pw_const_proj_fefct, TFEVectFunct2D *pw_const_param_fe,
TFEVectFunct2D *pw_const_param_deriv_fe, TOutput2D *&Output,
int level, int mg_level)
{
  // prepare data output
  if (level==mg_level)
  {
    if (!TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)
    {
      Output = new TOutput2D(3, 3, 1, 0, Domain);
      Output->AddFEFunction(conc);
      if (pw_const_proj_fefct != NULL)
      {
        Output->AddFEFunction(pw_const_proj_fefct->GetComponent(0));
        Output->AddFEFunction(pw_const_proj_fefct->GetComponent(1));
      }
    }
    else
    {
      Output = new TOutput2D(3, 5, 3, 1, Domain);
      Output->AddFEFunction(conc);
      Output->AddFEFunction(conc_adjoint);
      if (TDatabase::ParamDB->SOLD_ADJOINT)
      {
        Output->AddFEVectFunct(pw_const_param_fe);
        Output->AddFEVectFunct(pw_const_param_deriv_fe);
      }
      else
      {
        Output->AddFEFunction(pw_const_param_fe->GetComponent(0));
        //Output->AddFEFunction(pw_const_param_deriv_fe->GetComponent(0));
      }
    }
    if (TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
      Output->AddFEVectFunct(velocity);
  }
}


/******************************************************************************/
// gets the information from the coarse grid for CDAdap_2D.C
/****************F**************************************************************/
void GetInformationFromCoarseGridCD_2D(TFESpace2D **USpaces, TFESpace2D *old_u_space,
TFEFunction2D **UArray, TFEFunction2D *old_u,
int level, int mg_type, int mg_level, int FirstSolve,
double *sol, double *oldsol)
{
  // prolongation, to get a good starting iterate
  // HAS TO BE TESTED FOR ADAPTIVE GRIDS
  if ((level==mg_level) &&                        // only on finest level
    level>FirstSolve+mg_type &&                   // there is a coarser level
    (TDatabase::ParamDB->P5!=4711) &&
    (!TDatabase::ParamDB->SOLD_ADJOINT)&&         // not for algorithm with adjoint problem
    (level<=TDatabase::ParamDB->UNIFORM_STEPS))   // only for uniform grids
  {
    OutPut("prolongate solution from next coarser level" << endl);
    Prolongate(old_u_space, USpaces[mg_level],
      old_u->GetValues(), UArray[mg_level]->GetValues(), oldsol);
    // copy current solution for assembling the nonlinear terms
    // memcpy(oldsol,sol, N_Unknowns*SizeOfDouble);
    if (mg_type==1)
    {
      delete old_u->GetValues();
      delete old_u;
      delete old_u_space;
    }
  }                                               // end of prolongate
  else
  {
    if (level==mg_level)
      OutPut("Starting iterate is zero !"<<endl);
  }
}

void MeasureErrorsCD_2D(TFESpace2D *fe_space, TFEFunction2D *conc,
DoubleFunct2D *Exact,
BoundCondFunct2D **BoundaryConditions,
BoundValueFunct2D **BoundaryValues,
CoeffFct2D *BilinearCoeffs, int level, int mg_level, int FirstSolve,
int *N_Uarray, double *sol, double *rhs,
double *l2, double *h1, double *sd, double *energy, double *l_inf, double *lp)
{
  double errors[5], estimated_global_error[7];
  TFESpace2D *fesp[1];
  TFEFunction2D *fefct[2];
  TAuxParam2D *aux;
  ErrorMethod2D *ErrorMeth;

  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
  fesp[0] = fe_space;
  aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
  conc->GetErrors(Exact, 3, AllDerivatives, 5, SDFEMErrors,
    BilinearCoeffs, aux, 1, fesp, errors);
  delete aux;
  l2[level] = errors[0];
  h1[level] = errors[1];
  sd[level] = errors[2];
  energy[level] = sqrt(errors[1]*errors[1]/TDatabase::ParamDB->RE_NR + errors[0]*errors[0]);
  l_inf[level] = errors[4];
  if(TDatabase::ParamDB->DISCTYPE == LOCAL_PROJECTION)
  {
    if(TDatabase::ParamDB->LP_FULL_GRADIENT)
    {
      lp[level] = sqrt(h1[level]*h1[level]/TDatabase::ParamDB->RE_NR + l2[level]*l2[level] + 0);
      /*UltraLocalError(u, Exact,
                   TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF,
                   TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT));*/
      //TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE));
    }

    if(TDatabase::ParamDB->LP_STREAMLINE)
    {
      OutPut("error computation: local projection stabilisation in streamline direction ");
      OutPut("is currently not available." << endl);
    }
  }

  if (level>FirstSolve)
  {
    OutPut(level<< " L2: " << errors[0]<< " order " <<  log(l2[level-1]/l2[level])/ln2);
    OutPut( " H1-semi: " << errors[1] << " order " << log(h1[level-1]/h1[level])/ln2);
    OutPut( " energy: " << energy[level] << " order " << log(energy[level-1]/energy[level])/ln2 << endl);
    OutPut(level << " SD: " << errors[2] << " order " << log(sd[level-1]/sd[level])/ln2 <<  " cross wind: " << errors[3] << endl);
    OutPut(level<< " L_inf: " << errors[4] << " order " << log(l_inf[level-1]/l_inf[level])/ln2 << endl);

    if(TDatabase::ParamDB->DISCTYPE == LOCAL_PROJECTION)
    {
      OutPut( "LP: " << lp[level] << " order " << log(lp[level-1]/lp[level])/ln2 << endl);
    }
    OutPut("SD parts: H1 " << sqrt( 1.0/TDatabase::ParamDB->RE_NR) * h1[level]
      << " l2 " <<  sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY) * l2[level]
      << " sd " << sqrt(sd[level]*sd[level] - 1.0/TDatabase::ParamDB->RE_NR * h1[level] * h1[level] -  TDatabase::ParamDB->INTERNAL_COERCIVITY * l2[level] * l2[level])
      << " " << sqrt(sd[level]*sd[level] - 1.0/TDatabase::ParamDB->RE_NR * h1[level] * h1[level]
      -  TDatabase::ParamDB->INTERNAL_COERCIVITY * l2[level] * l2[level])/ (sqrt( 1.0/TDatabase::ParamDB->RE_NR) * h1[level]) <<  endl);
  }
  else
  {
    OutPut(level<< " L2: " << errors[0]);
    OutPut( " H1-semi: " << errors[1]);
    OutPut( " energy: " << energy[level] << endl);
    OutPut(level<< " SD: " << errors[2] << " cross wind: " << errors[3] << endl);
    OutPut(level<< " L_inf: " << errors[4] << endl);
  }
  if(TDatabase::ParamDB->DISCTYPE == LOCAL_PROJECTION)
  {
    OutPut( "LP: " << lp[level] << endl);
  }

  // compute error to interpolant, for paper with J. Novo,
  if (TDatabase::ParamDB->P9 ==  4711)
  {
    if(TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY != 120814)
    {
      // save solution
      memcpy(rhs,sol,N_Uarray[mg_level]*SizeOfDouble);
      // interpolate analytical solution
      conc->Interpolate(Exact);
    }

    // compute errors
    MultiIndex2D Derivatives_SD[5] = { D10, D01, D00, D20, D02 };
    fesp[0] = fe_space;
    aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
    conc->GetErrors(Exact, 5, Derivatives_SD, 4, SDFEMErrorsInterpolant,
      BilinearCoeffs, aux, 1, fesp, errors);

    errors[4] = ComputeErrorToInterpolantL2Edges(fesp[0], conc, Exact, 3, AllDerivatives,
      BilinearCoeffs, BoundaryConditions,
      BoundaryValues, aux, 1, fesp, &estimated_global_error[0]);
    delete aux;

    // recover solution
    if(TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY != 120814)
      memcpy(sol,rhs,N_Uarray[mg_level]*SizeOfDouble);

    OutPut(level << "Interpolant SUPG: " << errors[0] << " (" << 2*sd[level]*sd[level]/ (errors[0]*errors[0])<< ") "
      << errors[1] << " (" << 2*sd[level]*sd[level]/(errors[1]*errors[1]) << ") "
      << errors[4] << " (" << 2*sd[level]*sd[level]/(errors[4]*errors[4]) << ") " << endl);
    OutPut("additional terms SUPG: " << errors[2]  << " (" << errors[2]/ sd[level]<< ") "
      <<  errors[3]  << " (" << errors[3]/ sd[level]<< ") " << endl);
  }

  // measure L1 error to known solution (check whether this routine is correct
  //aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
  // conc->GetErrors(Exact, 1, ZeroDerivative, 1, L1Error, BilinearCoeffs, aux, 1, fesp, errors);
  // delete aux;
  //OutPut(setw(20) << setprecision(20) << "L1: " << errors[0] << endl);
  //OutPut("L1: " << errors[0] << endl);
}

/******************************************************************************/
// function restriction like Laplacian interpolation
// only for P1, Q1 on nested grids
/****************F**************************************************************/

void RestrictFunctionLikeLaplacianInterpolation(TFEFunction2D *fine, TFEFunction2D *coarse)
{
  int i, j, k, l, N_Cells, N_V, N_child, loc_dof_coarse, loc_dof_fine, no_fine;
  int *GlobalNumbers_coarse, *BeginIndex_coarse, *GlobalNumbers_fine, *BeginIndex_fine;
  int *dof, *dof_fine;
  double eps = 1e-6, x, y, x_c, y_c;
  TCollection *coarse_coll, *fine_coll;
  TBaseCell *cell, *child;
  TVertex *vertex;
  FE2D CurrentElement;

  // get collection
  coarse_coll = coarse->GetFESpace2D()->GetCollection();
  // get number of cells
  N_Cells = coarse_coll->GetN_Cells();
  // get information on dof for coarse fe function
  GlobalNumbers_coarse  = coarse->GetFESpace2D()->GetGlobalNumbers();
  BeginIndex_coarse  = coarse->GetFESpace2D()->GetBeginIndex();

  // get collection
  fine_coll = fine->GetFESpace2D()->GetCollection();
  // get number of cells
  k = fine_coll->GetN_Cells();
  // set clip board
  for (i=0;i<k;i++)
    fine_coll->GetCell(i)->SetClipBoard(i);

  // get information on dof for fine fe function
  GlobalNumbers_fine  = fine->GetFESpace2D()->GetGlobalNumbers();
  BeginIndex_fine  = fine->GetFESpace2D()->GetBeginIndex();

  // loop over coarse grid
  for(i=0;i<N_Cells;i++)
  {
    cell = coarse_coll->GetCell(i);
    // the array which gives the mapping of the local to the global d.o.f.
    dof = GlobalNumbers_coarse+BeginIndex_coarse[i];
    // finite element on the mesh cell
    CurrentElement = coarse->GetFESpace2D()->GetFE2D(i, cell);
    //OutPut("cell  " << i << endl);
    N_V = cell->GetN_Vertices();
    // loop over the vertices
    for (j=0;j<N_V;j++)
    {
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x, y);
      //OutPut(" vertex " << j)
      switch(CurrentElement)
      {
        // P_1, Q_1
        case C_P1_2D_T_A:
          loc_dof_coarse = dof[j];
          break;
        case C_Q1_2D_Q_A:
        case C_Q1_2D_Q_M:
          if (j<2)
            loc_dof_coarse = dof[j];
          else
          {
            if (j==2)
              loc_dof_coarse = dof[3];
            else
              loc_dof_coarse = dof[2];
          }
          break;
          break;
        default:
          OutPut("RestrictFunctionLikeLaplacianInterpolation not implemented" << endl);
          exit(4711);
      }
      // find this vertex on the finer grid
      N_child = cell->GetN_Children();
      // loop over the children
      for (k=0;k<N_child;k++)
      {
        // get child
        child = cell->GetChild(k);
        // loop over the vertices
        // ASSUMED THAT CHILDREN HAVE SAME NUMBER OF VERTICES
        for (l=0;l<N_V;l++)
        {
          vertex = child->GetVertex(l);
          vertex->GetCoords(x_c, y_c);
          if ((fabs(x-x_c)<eps) && (fabs(y-y_c)<eps ))
          {
            no_fine = child->GetClipBoard();
            //OutPut(" child vertex found on cell " << no_fine << endl);
            // the array which gives the mapping of the local to the global d.o.f.
            dof_fine = GlobalNumbers_fine+BeginIndex_fine[no_fine];
            // ASSUMED THAT ON FINE GRID SAME FE AS ON COARSE GRID
            switch(CurrentElement)
            {
              // P_1, Q_1
              case C_P1_2D_T_A:
                loc_dof_fine = dof_fine[l];
                break;
              case C_Q1_2D_Q_A:
              case C_Q1_2D_Q_M:
                if (j<2)
                  loc_dof_fine = dof_fine[l];
                else
                {
                  if (j==2)
                    loc_dof_fine = dof_fine[l];
                  else
                    loc_dof_fine = dof_fine[l];
                }
                break;
              default:
                loc_dof_fine = dof_fine[l];
                break;
            }
            //OutPut(" corresponding dofs found " << loc_dof_coarse << " " <<  loc_dof_fine<< endl);
            coarse->GetValues()[loc_dof_coarse] = fine->GetValues()[loc_dof_fine];
            break;
          }
        }
      }
    }
  }
}

void MeasureErrorsToFinestGridSolutionCD_2D(TFESpace2D **fe_spaces, TFEFunction2D **concs,
DoubleFunct2D *Exact,
BoundCondFunct2D **BoundaryConditions,
BoundValueFunct2D **BoundaryValues,
CoeffFct2D *BilinearCoeffs, int level, int mg_level, int FirstSolve,
int *N_Uarray, double *sol, double *rhs,
double *l2, double *h1, double *sd, double *energy, double *l_inf, double *lp,
double **&sol_array, int &first_comp)
{
  int i, j, len;
  double *current_solution, *aux_vec, *dummy_array, para;
  TFEFunction2D *conc_finest;
  char String[] = "concentration";

  if ((TDatabase::ParamDB->SOLVER_TYPE != GMG) || (TDatabase::ParamDB->SC_MG_TYPE_SCALAR !=0))
  {
    OutPut("MeasureErrorsFinestGridSolutionCD_2D not implemented for this type of solver" << endl);
    return;
  }

  if (first_comp)
  {
    sol_array = new double* [TDatabase::ParamDB->LEVELS+1];
    first_comp = 0;
  }

  // save current solution
  len = concs[mg_level]->GetLength();
  aux_vec = new double[len];
  current_solution = new double[len];
  memcpy(current_solution, concs[mg_level]->GetValues(), len * SizeOfDouble);
  sol_array[mg_level] = current_solution;
  conc_finest =  new TFEFunction2D(fe_spaces[mg_level], String, String, current_solution, len);

  // dummy array
  dummy_array = new double[TDatabase::ParamDB->LEVELS+1];

  // compute restriction
  for (i=mg_level;i>FirstSolve;i--)
  {
    if (TDatabase::ParamDB->ANSATZ_ORDER == 1)
      RestrictFunctionLikeLaplacianInterpolation(concs[i],concs[i-1]);
    else
      RestrictFunction(fe_spaces[i-1],fe_spaces[i],concs[i-1]->GetValues(),
        concs[i]->GetValues(),aux_vec);
    //for (j=0;j<N_Uarray[i-1];j++)
    //OutPut(j << " rest " << concs[i-1]->GetValues()[j]<< " sol " << sol_array[i-1][j] <<  endl);
  }

  // interpolation part
  for (i=mg_level;i>=FirstSolve;i--)
  {
    // copy coarse grid solutions to the fefunctions
    len =  concs[i]->GetLength();
    TDatabase::ParamDB->INTERNAL_LEVEL = i;
    // prolongate to finest grid
    for (j=i;j<mg_level;j++)
    {
      //memcpy(concs[i]->GetValues(), sol_array[i], len * SizeOfDouble);
      Prolongate(fe_spaces[j], fe_spaces[j+1], concs[j]->GetValues(), concs[j+1]->GetValues(),
        aux_vec);
    }
    // compute difference to solution on finest grid
    len  =  concs[mg_level]->GetLength();
    Daxpy(len, -1, current_solution, concs[mg_level]->GetValues());
    OutPut("interpolation_error of level " << i << " to level " << mg_level << endl);

    MeasureErrorsCD_2D(fe_spaces[mg_level], concs[mg_level],
      ExactNull, BoundaryConditions, BoundaryValues, BilinearCoeffs,
      mg_level, mg_level, FirstSolve, N_Uarray, NULL, NULL,
      dummy_array, dummy_array, dummy_array, dummy_array, dummy_array, dummy_array);

  }

  // no need to compute values for paper with J. Novo
  para = TDatabase::ParamDB->P9;
  TDatabase::ParamDB->P9 = 0;

  for (i=FirstSolve;i<=mg_level;i++)
  {
    // copy coarse grid solutions to the fefunctions
    len =  concs[i]->GetLength();
    memcpy(concs[i]->GetValues(), sol_array[i], len * SizeOfDouble);
    TDatabase::ParamDB->INTERNAL_LEVEL = i;

    // prolongate to finest grid
    for (j=i;j<mg_level;j++)
    {
      Prolongate(fe_spaces[j], fe_spaces[j+1], concs[j]->GetValues(), concs[j+1]->GetValues(),
        aux_vec);
    }
    // compute difference to solution on finest grid
    len  =  concs[mg_level]->GetLength();
    Daxpy(len, -1, current_solution, concs[mg_level]->GetValues());
    OutPut("error of level " << i << " to level " << mg_level << endl);

    MeasureErrorsCD_2D(fe_spaces[mg_level], concs[mg_level],
      ExactNull, BoundaryConditions, BoundaryValues, BilinearCoeffs,
      mg_level, mg_level, FirstSolve, N_Uarray, NULL, NULL,
      dummy_array, dummy_array, dummy_array, dummy_array, dummy_array, dummy_array);
  }
  // rewrite current solution
  memcpy(concs[mg_level]->GetValues(), current_solution, len * SizeOfDouble);
  TDatabase::ParamDB->P9 = para;

  delete aux_vec;
  delete conc_finest;
}


/******************************************************************************/
// restricts the solution on the currently finest grid and saves it
// this routine works only for
// SOLVER_TYPE: 1  (geometric multigrid)
// SC_MG_TYPE_SCALAR: 0
/****************F**************************************************************/

void SaveRestrictedFinestGridSolutionCD_2D(TFEFunction2D **concs, int mg_level,
int FirstSolve)
{
  int i, len, save_N_Unknowns[1];
  double *save_sol[1];
  TFEFunction2D *conc_finest;
  char *SaveDataFileName;
  std::ostringstream os;
  //string SaveDataFileName, filename;

  SaveDataFileName = TDatabase::ParamDB->SAVE_DATA_FILENAME;

  if ((TDatabase::ParamDB->SOLVER_TYPE != GMG) || (TDatabase::ParamDB->SC_MG_TYPE_SCALAR !=0))
  {
    OutPut("SaveRestrictedFinestGridSolutionCD_2D not implemented for this type of solver" << endl);
    return;
  }

  for (i=mg_level;i>FirstSolve;i--)
  {
    RestrictFunctionLikeLaplacianInterpolation(concs[i],concs[i-1]);
    // check if function should be saved
    len  = concs[i-1]->GetLength();
    // save only solutions on coarser grids
    if (len<2e4)
    {
      // construct the name for saving
      os.seekp(std::ios::beg);
      os << SaveDataFileName << "_fromlevel_" << mg_level <<  "_" << len <<ends;
      OutPut("save coarse grid solution in file "<< os.str().c_str() << endl);
      save_sol[0] = concs[i-1]->GetValues();
      save_N_Unknowns[0] = len;
      SaveData((char *)(os.str().c_str()),1,save_sol,save_N_Unknowns);
    }
  }
}

/******************************************************************************/
// reads the restricted solution from the finest grid that was saved
// this routine works only for
// SOLVER_TYPE: 1  (geometric multigrid)
// SC_MG_TYPE_SCALAR: 0
/****************F**************************************************************/
void ReadRestrictedFinestGridSolutionCD_2D(TFEFunction2D *concs)
{
  int i, len, save_N_Unknowns[1];
  double *save_sol[1];
  TFEFunction2D *conc_finest;
  char *ReadDataFileName;
  std::ostringstream os;

  if (!TDatabase::ParamDB->READ_DATA)
    return;
  ReadDataFileName = TDatabase::ParamDB->READ_DATA_FILENAME;

  // construct the name for reading
  len  = concs->GetLength();
  if (len>2e4)
    return;
  os.seekp(std::ios::beg);
  os << ReadDataFileName << len  <<ends;
  OutPut("read coarse grid solution from file "<< os.str().c_str() << endl);
  save_sol[0] = concs->GetValues();
  save_N_Unknowns[0] = len;
  ReadData((char *)(os.str().c_str()),1,save_sol,save_N_Unknowns);
}

/******************************************************************************/
// Assemble_CD_2D
// assembles the matrices
/******************************************************************************/
void Assemble_CD_2D(TDomain *Domain,
TFEFunction2D **UArray,
TFEVectFunct2D **VelocityArray,
TFEVectFunct2D *pw_const_proj_fefct,
double **RhsArray, TSquareMatrix2D **MatricesA,
CoeffFct2D *Coeffs,
BoundCondFunct2D **BoundaryConditions,
BoundValueFunct2D **BoundaryValues,
TDiscreteForm2D **DiscreteForms,
TSquareMatrix2D *sqmatrixAadjoint,
CheckWrongNeumannNodesFunct2D **CheckWrongNeumannNodesFct,
int low, int mg_level, int mg_type, int i,
int &N_neum_to_diri, int* &neum_to_diri,
int* &neum_to_diri_bdry,
double* &neum_to_diri_param, double *rhs_edge,
int problem_id, TFEFunction2D **param_func)
{
  int k, N_U, N_Active, N_NonActive, type, flag, N_fesp, N_fefct;
  int N_ParamFct, N_FEValues, N_Params;
  int FEFctIndex[10], BeginParam[1];
  double *rhs, *RHSs[1], tmp;
  TFESpace2D *fesp[3], *ferhs[1], *conc_space;
  TFEFunction2D  *fefct[5];
  TSquareMatrix2D *SQMATRICES[1];
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TAuxParam2D *aux;
  TDiscreteForm2D *DiscreteForm;
  BoundCondFunct2D *BdCond;
  BoundValueFunct2D *BdVal;
  MultiIndex2D FEMultiIndex_All_Deriv[5] = { D00, D10, D01, D20, D02 };
  MultiIndex2D FEMultiIndex_Fct[1] = {D00};
  MultiIndex2D FEMultiIndex[10];
  TCollection *coll;
  ParamFct *ParamFct[1];

  // loop over the levels
  for(k=low;k<=mg_level;k++)
  {
    TDatabase::ParamDB->INTERNAL_MOMENT = k;
    rhs = RhsArray[k];
    conc_space = UArray[k]->GetFESpace2D();
    N_U = UArray[k]->GetLength();
    N_Active = conc_space->GetActiveBound();
    N_NonActive = N_U - N_Active;

    // initialize rhs
    RHSs[0] = rhs;
    memset(rhs, 0, N_U*SizeOfDouble);
    ferhs[0] = conc_space;

    // initialize matrix
    SQMATRICES[0] = MatricesA[k];
    SQMATRICES[0]->Reset();

    // find discrete form
    if ((mg_type==1) && (k<i+1))
    {
      DiscreteForm = DiscreteForms[2];
      OutPut("Upwind ");
    }
    else
      switch(TDatabase::ParamDB->DISCTYPE)
      {
        case GALERKIN:
        case LOCAL_PROJECTION:
          DiscreteForm = DiscreteForms[0];
          break;

        case SDFEM:
          if (!((TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)&&(TDatabase::ParamDB->SOLD_ADJOINT)))
            DiscreteForm = DiscreteForms[1];
          else
          DiscreteForm = DiscreteForms[7];
        type = TDatabase::ParamDB->SDFEM_TYPE;
        // standard parameter choice on coarser levels
        if ((TDatabase::ParamDB->SDFEM_TYPE == 100)&&(k<mg_level))
        {
          TDatabase::ParamDB->SDFEM_TYPE = 2;
      }
      break;

      case UPWIND:
        DiscreteForm = DiscreteForms[2];
        OutPut("Upwind ");
        break;

      case LOCAL_PROJECTION_2_LEVEL:
        //DiscreteForm = DiscreteForms[1];
        //TDatabase::ParamDB->SDFEM_TYPE = 1;
        //TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 0;
        DiscreteForm = DiscreteForms[17];
        LocalProjectionCoarseGridQ0(UArray[k], pw_const_proj_fefct->GetComponent(0), Coeffs, 0);
        LocalProjectionCoarseGridQ0(UArray[k], pw_const_proj_fefct->GetComponent(1), Coeffs, 1);
        tmp = TDatabase::ParamDB->LP_CROSSWIND_COEFF;
        TDatabase::ParamDB->LP_CROSSWIND_COEFF = 0;
        break;

      default:
        OutPut("Unknown DISCTYPE" << endl);
        exit(4711);
    }

    // compute aux object
    // initialize parameters for aux object
    fesp[0] = conc_space;
    N_fesp = 1;
    fefct[0] = UArray[k];
    N_fefct = 1;
    N_ParamFct = 0;
    N_FEValues = 0;
    N_Params = 0;

    // different stabilizations that cannot be used together
    // SOLD
    if (DiscreteForm == DiscreteForms[7])
    {
      // always one parameter function
      ParamFct[N_ParamFct] = Params_All_Deriv;
      N_ParamFct = 1;
      N_FEValues += 5;
      FEFctIndex[N_Params] = 0;
      FEFctIndex[N_Params + 1] = 0;
      FEFctIndex[N_Params + 2] = 0;
      FEFctIndex[N_Params + 3] = 0;
      FEFctIndex[N_Params + 4] = 0;
      FEMultiIndex[N_Params] = D00;
      FEMultiIndex[N_Params + 1] = D10;
      FEMultiIndex[N_Params + 2] = D01;
      FEMultiIndex[N_Params + 3] = D20;
      FEMultiIndex[N_Params + 4] = D02;
      N_Params += 5;
      // always first parameter
      BeginParam[0] = 0;
    }
    else
    {
      // 2level LPS
      if (DiscreteForm == DiscreteForms[17])
      {
        fesp[N_fesp] = pw_const_proj_fefct->GetFESpace2D();
        N_fesp += 1;
        fefct[N_fefct] = pw_const_proj_fefct->GetComponent(0);
        fefct[N_fefct + 1] = pw_const_proj_fefct->GetComponent(1);
        N_fefct += 2;
        // always one parameter function
        ParamFct[N_ParamFct] = Params_Sol_And_Pw_Const_Proj;
        N_ParamFct = 1;
        N_FEValues += 4;
        FEFctIndex[N_Params] = 0;
        FEFctIndex[N_Params + 1] = 0;
        FEFctIndex[N_Params + 2] = N_fefct-2;
        FEFctIndex[N_Params + 3] = N_fefct-1;
        FEMultiIndex[N_Params] = D10;
        FEMultiIndex[N_Params + 1] = D01;
        FEMultiIndex[N_Params + 2] = D00;
        FEMultiIndex[N_Params + 3] = D00;
        N_Params += 4;
        // always first parameter
        BeginParam[0] = 0;
      }
    }

    // convection is velocity field or temperature equation for Boussinesq
    if ((TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)||(problem_id == 1))
    {
      if (k==mg_level)
      {
        if (TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD)
        {
          // add velocity space
          fesp[N_fesp] = VelocityArray[k]->GetFESpace2D();
          N_fesp += 1;
          // fefct for both components
          fefct[N_fefct] = VelocityArray[k]->GetComponent(0);
          fefct[N_fefct + 1] = VelocityArray[k]->GetComponent(1);
          N_fefct += 2;
        }
        else
        {
          // add velocity space
          fesp[N_fesp] = param_func[0]->GetFESpace2D();
          N_fesp += 1;
          // fefct for both components
          fefct[N_fefct] = param_func[0];
          fefct[N_fefct + 1] = param_func[1];
          N_fefct += 2;
        }

        // always one parameter function (more efficient than several functions)
        if (N_ParamFct==0)
          ParamFct[N_ParamFct] = Params_Velo;
        else
        {
          if (ParamFct[0] == Params_All_Deriv)
            ParamFct[0] =  Params_All_Deriv_And_Velo;
          else
          {
            if (ParamFct[0] == Params_Sol_And_Pw_Const_Proj)
              ParamFct[0] =  Params_Sol_And_Pw_Const_Proj_And_Velo;
            else
            {
              OutPut("Error in finding the parameter function for the aux object" << endl);
              exit(4711);
            }
          }
        }

        N_ParamFct = 1;
        N_FEValues += 2;
        FEFctIndex[N_Params] = N_fefct-2;
        FEFctIndex[N_Params + 1] = N_fefct-1;
        FEMultiIndex[N_Params] = D00;
        FEMultiIndex[N_Params + 1] = D00;
        // just that example knows where velocity starts
        TDatabase::ParamDB->INTERNAL_START_PARAM = N_Params;
        N_Params += 2;
        // always first parameter
        BeginParam[0] = 0;
      }
      else
      {
        OutPut("velocity field as convection not yet implemented on coarser grids"<<endl);
        OutPut(4711);
      }
    }

    if (N_ParamFct)
      aux =  new TAuxParam2D(N_fesp, N_fefct, N_ParamFct, N_FEValues, fesp, fefct,
        ParamFct, FEFctIndex, FEMultiIndex, N_Params, BeginParam);
    else
      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

    SetPolynomialDegree();
    BdCond = BoundaryConditions[0];
    BdVal = BoundaryValues[0];
    // set homogeneous Neumann boundary conditions
    if ((TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD)
      || (TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM))
    {
      BdCond = BoundaryConditions[1];
      BdVal = BoundaryValues[1];
    }
    // assemble
    Assemble2D(N_fesp, fesp,
      1, SQMATRICES,
      0, NULL,
      1, RHSs, ferhs,
      DiscreteForm,
      &BdCond,
      &BdVal,
      aux);
    if (TDatabase::ParamDB->SC_VERBOSE>1)
      OutPut("Assembling done on mg-level: "<< k <<endl);

    // not for upwind
    if (DiscreteForm != DiscreteForms[2])
    {
      // DG and CIP
      switch (TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS)
      {
        case 1:  Assemble2D_CIP(Coeffs,1, fesp,
          1, SQMATRICES,
          0, NULL,
          1, RHSs, ferhs,
          BoundaryConditions,
          BoundaryValues,
          aux);
        break;
        case 2:  Assemble2D_DG(Coeffs,1, fesp,
          1, SQMATRICES,
          0, NULL,
          1, RHSs, ferhs,
          BoundaryConditions,
          BoundaryValues,
          aux);
        break;
      }

      // reset SUPG parameter
      if (TDatabase::ParamDB->DISCTYPE == SDFEM)
      {
        TDatabase::ParamDB->SDFEM_TYPE = type;
      }
      // LPS ???
      if(TDatabase::ParamDB->DISCTYPE == LOCAL_PROJECTION)
      {
        if(TDatabase::ParamDB->LP_FULL_GRADIENT)
          UltraLocalProjection(MatricesA[k], FALSE, Coeffs);

        if(TDatabase::ParamDB->LP_STREAMLINE)
        {
          flag = 0;
          if (TDatabase::ParamDB->LP_CROSSWIND)
          {
            TDatabase::ParamDB->LP_CROSSWIND = 0;
            flag = 1;
          }
          UltraLocalProjectionStreamlinePLaplacian(MatricesA[k], UArray[k], Coeffs);
          if (flag)
            TDatabase::ParamDB->LP_CROSSWIND = 1;
        }
      }
    }
    // apply upwind
    if (DiscreteForm == DiscreteForms[2])
    {
      UpwindForConvDiff(Coeffs, SQMATRICES[0],RHSs[0], fesp[0],DiscreteForms[2], 
                        //VelocityArray[k]);
                        NULL, NULL, 0);
      cout << "UPWINDING DONE : level " << k << endl;
    }                                             // endif
    delete aux;

    // compute adjoint matrix, only on the fines level
    if ((TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)&&(k==mg_level))
    {
      ComputeAdjointMatrix(MatricesA[k],sqmatrixAadjoint);
    }

    // save rhs for the nonlinear iteration
    if ((TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD)&&(k==mg_level))
    {
      memcpy(rhs_edge,rhs,N_U*SizeOfDouble);
    }

    if ((TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD)
      || (TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM))
    {
      if (!TDatabase::ParamDB->INTERNAL_WRONG_NEUMANN_CHECKED)
      {
        //coll=Domain->GetCollection(It_EQ, k+TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR);
        // coll=conc_space->GetCollection();
        CheckWrongNeumannNodesFct[0](conc_space->GetCollection(), conc_space, N_neum_to_diri, neum_to_diri,
          neum_to_diri_bdry, neum_to_diri_param);
        // for non-multigrid solvers only one matrix there, always the same numbers
        if ((TDatabase::ParamDB->SOLVER_TYPE!=0)&&(k==mg_level))
          TDatabase::ParamDB->INTERNAL_WRONG_NEUMANN_CHECKED = 1;
        //OutPut(k << " check diri " << N_neum_to_diri << endl);
      }
      if (N_neum_to_diri)
      {
        //OutPut("set diri " << N_neum_to_diri << endl);
        SetDirichletNodesFromNeumannNodes(SQMATRICES, RHSs[0], UArray[k]->GetValues(),
          N_neum_to_diri, neum_to_diri,
          neum_to_diri_bdry, neum_to_diri_param,
          BoundaryValues[0]);
      }
    }

    if (DiscreteForm == DiscreteForms[17])
    {
      // update rhs
      LocalProjectionCrossWindCoarseGridQ0(Domain, k, UArray[k],
        pw_const_proj_fefct->GetComponent(1),
        Coeffs, rhs, 1);
      TDatabase::ParamDB->LP_CROSSWIND_COEFF = tmp;
    }
  }
}


/******************************************************************************/
// InitializeParametersForNonlinearMethods
/******************************************************************************/
void InitializeParametersForNonlinearMethods(TFESpace2D **fesp,
TFESpace2D *Uspace, TFESpace2D *sold_space,
TSquareMatrix2D **SQMATRICES,
TSquareMatrix2D *matrixA,
int N_Unknowns, int &m, int &linite,
int &max_it,
int &compute_matrix_D,
int mg_level,
double &lin_red,
double &res_norm_min,
double &nonlin_min_res,
double &omega,
double &omega_max,
double &oldres_stepm1,
double *sol, double *oldsol,
double* &defect, double * rhs_edge)
{
  int sold_parameter_type = TDatabase::ParamDB->SOLD_PARAMETER_TYPE;
  // iteration counter
  m = 0;
  // parameters for solution of linear and nonlinear equation
  linite = 0;
  max_it = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
  TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR_SOLD;
  compute_matrix_D = 1;
  lin_red =  TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR;
  TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR =  TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR_SOLD;
  res_norm_min = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
  TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR_SOLD;
  nonlin_min_res = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR;

  // first damping parameter
  if (TDatabase::ParamDB->P5==1234)
  {
    // fixed damping parameter
    omega = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SCALAR-0.05*mg_level;
    OutPut("constant damping factor changed to " << omega << endl);
    omega_max = omega;
  }
  else
    // from the next coarser level
    omega = omega_max;

  // residual of the previous step
  oldres_stepm1 = 1.0;
  // fe spaces for assembling
  fesp[0] = Uspace;
  if (TDatabase::ParamDB->SOLD_TYPE)
    fesp[1] = sold_space;
  // the matrix
  SQMATRICES[0] = matrixA;
  // save solution
  memcpy(oldsol,sol, N_Unknowns*SizeOfDouble);
  // define and initialize defect
  defect = new double[N_Unknowns];
  memset(defect,0,N_Unknowns*SizeOfDouble);
}


/*********************************************************************/
// AssembleMatrixforBilinearForm_SameCell
// computes an assembled matrix for a bilinear form for a class
// of base functions
// ATTENTION every element is recomputed from skretch!
/*********************************************************************/
// assembly - new Matrix for bilinear form
// coll - set of cells
// ansatz - FE basis for ansatz space
// test - FE basis for test space
// N_Ansatz - number of Ansatz functions
// N_Test - number of test functions
// Coeff - coefficient function
// functional - index of functional to be computed, see below
/*********************************************************************/
void AssembleMatrixforBilinearForm_SameCells(double** Assembly, TCollection *coll, TFEFunction2D **ansatz, TFEFunction2D **test, int N_Ansatz, int N_Test, CoeffFct2D *Coeff, int functional)
{
  const int N_Derivatives = 3;
  int i, j, k, l, m, N_, N_Cells, N_Fespaces, N_LocalUsedElements, N_Points;
  int N_Parameters;
  int N_Edges;
  int *GlobalNumbers_ansatz, *BeginIndex_ansatz, *GlobalNumbers_test, *BeginIndex_test;
  int *N_BaseFunct, * DOF, LocN_BF[2];
  double val_functional, val, value, w;
  double *weights, *xi, *eta, *Values_ansatz, *Values_test, *Orig, *aux1, *aux2, *aux3, *aux4;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], AbsDetjk[MaxN_QuadPoints_2D];
  double FEFunctValues_ansatz[MaxN_BaseFunctions2D],FEFunctValues_test[MaxN_BaseFunctions2D];
  double *Derivatives_ansatz[N_Ansatz][MaxN_QuadPoints_2D], *Derivatives_test[N_Test][MaxN_QuadPoints_2D];
  double *Param[MaxN_QuadPoints_2D], *CoeffArray[MaxN_QuadPoints_2D];
  double **OrigFEValues;
  double scale_coeff[ 5 ];
  double theta1 = TDatabase::TimeDB->THETA1;
  double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  TFESpace2D *space_ansatz, *space_test;
  TBaseCell *cell;
  FE2D LocalUsedElements[2], CurrentElement;
  BaseFunct2D *BaseFuncts, LocBF[2], BaseFunct;
  RefTrans2D RefTrans;
  MultiIndex2D NeededDerivatives[N_Derivatives] = {D00, D10, D01};
  bool SecondDer[2];

  // set zeros
  for( int i = 0; i< N_Test; i++ )
  {
    for( int j = 0; j <N_Ansatz; j++ )
    {
      Assembly[ i ][ j ] = 0.0;
    }
  }

  // get data from the database
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  // get fe spaces
  // assume all functions share same FESpace
  space_ansatz = ansatz[0]->GetFESpace2D();
  space_test   = test[0]->GetFESpace2D();
  N_Fespaces = 2;

  // get arrays for degrees of freedom
  GlobalNumbers_ansatz = space_ansatz->GetGlobalNumbers();
  BeginIndex_ansatz = space_ansatz->GetBeginIndex();
  GlobalNumbers_test = space_test->GetGlobalNumbers();
  BeginIndex_test = space_test->GetBeginIndex();
  // find number of quadrature points
  // it is assumed to be the same on all cells
  cell=coll->GetCell(0);
  // local ansatz space on this cell
  CurrentElement = space_ansatz->GetFE2D(0,cell);
  LocalUsedElements[0] = CurrentElement;
  // local basis function for space_ansatz
  LocN_BF[0] = N_BaseFunct[CurrentElement];
  LocBF[0] = BaseFuncts[CurrentElement];
  SecondDer[0] = true;
  // local test space on this cell
  CurrentElement = space_test->GetFE2D(0,cell);
  LocalUsedElements[1] = CurrentElement;
  // local basis function for space_test
  LocN_BF[1] = N_BaseFunct[CurrentElement];
  LocBF[1] = BaseFuncts[CurrentElement];
  SecondDer[1] = true;
  N_LocalUsedElements = N_Fespaces;

  // calculate values on original element
  // get reference transformation
  RefTrans = TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements,
    coll, cell, SecondDer,
    N_Points, xi, eta, weights, X, Y, AbsDetjk);

  // arrays for derivatives
  aux3 = new double [N_Points*N_Derivatives*N_Ansatz];
  for( i = 0 ; i < N_Ansatz ; i ++ )
  {
    for(j=0;j<N_Points;j++)
    {
      Derivatives_ansatz[ i ][ j ] = aux3 + j*N_Derivatives + i*N_Points*N_Derivatives;
    }
  }
  aux4 = new double [N_Points*N_Derivatives*N_Test];
  for( i = 0 ; i < N_Test ; i ++ )
  {
    for(j=0;j<N_Points;j++)
    {
      Derivatives_test[ i ][ j ] = aux4 + j*N_Derivatives + i*N_Points*N_Derivatives;
    }
  }

  // array for parameters
  N_Parameters = 1;
  aux1 = new double [N_Points*N_Parameters];
  for(j=0;j<N_Points;j++)
    Param[j] = aux1 + j*N_Parameters;

  // array for coefficients, maximal number is 20
  aux2 = new double [N_Points*20];
  for(j=0;j<N_Points;j++)
    CoeffArray[j] = aux2 + j*20;

  N_Cells = coll->GetN_Cells();
  // loop over the mesh cells
  for(i=0;i<N_Cells;i++)
  {
    // get cell
    cell=coll->GetCell(i);

    // local ansatz space on this cell
    CurrentElement = space_ansatz->GetFE2D(i,cell);
    LocalUsedElements[0] = CurrentElement;
    // local basis function for space_ansatz
    LocN_BF[0] = N_BaseFunct[CurrentElement];
    LocBF[0] = BaseFuncts[CurrentElement];
    SecondDer[0] = TRUE;
    // local test space on this cell
    CurrentElement = space_test->GetFE2D(i,cell);
    LocalUsedElements[1] = CurrentElement;
    // local basis function for space_test
    LocN_BF[1] = N_BaseFunct[CurrentElement];
    LocBF[1] = BaseFuncts[CurrentElement];
    SecondDer[1] = TRUE;

    N_LocalUsedElements = N_Fespaces;
    // calculate values on original element
    // get reference transformation
    RefTrans = TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements,
      coll, cell, SecondDer,
      N_Points, xi, eta, weights, X, Y, AbsDetjk);

    // calculate all needed derivatives of this FE function
    // ansatz function
    CurrentElement =  space_ansatz->GetFE2D(i,cell);
    // basis functions
    BaseFunct = BaseFuncts[CurrentElement];
    // no. of basis functions
    N_ = N_BaseFunct[CurrentElement];
    // dof of ansatz function in current mesh cell
    DOF = GlobalNumbers_ansatz + BeginIndex_ansatz[i];

    // build derivatives for all ansatz functions

    for( m = 0; m < N_Ansatz; m ++ )
    {
      Values_ansatz = ansatz[ m ]->GetValues();
      // fe values of ansatz functions
      for(l=0;l<N_;l++)
        FEFunctValues_ansatz[l] = Values_ansatz[DOF[l]];
      // compute values for all derivatives
      // in all quadrature points
      // in original mesh cell
      for(k=0;k<N_Derivatives;k++)                // for all derivatives
      {                                           // get values in original cell
        OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,NeededDerivatives[k]);
        for(j=0;j<N_Points;j++)                   // for all quadrature points
        {
          Orig = OrigFEValues[j];                 // value in original cell
          value = 0;
          for(l=0;l<N_;l++)
          {                                       // for all basis functions
            // accumulate value of derivative in point j
            value += FEFunctValues_ansatz[l] * Orig[l];
          }
          Derivatives_ansatz[m][j][k] = value;    // for k-th derivative
        }                                         // endfor j
      }                                           // endfor k
    }
    // calculate all needed derivatives of this FE function
    // ansatz function
    CurrentElement =  space_test->GetFE2D(i,cell);
    // basis functions
    BaseFunct = BaseFuncts[CurrentElement];
    // no. of basis functions
    N_ = N_BaseFunct[CurrentElement];
    // dof of ansatz function in current mesh cell
    DOF = GlobalNumbers_test + BeginIndex_test[i];

    // build derivatives for all test functions
    for( m = 0; m < N_Test; m ++ )
    {
      Values_test = test[ m ]->GetValues();
      // fe values of ansatz functions
      for(l=0;l<N_;l++)
        FEFunctValues_test[l] = Values_test[DOF[l]];

      // compute values for all derivatives
      // in all quadrature points
      // in original mesh cell
      for(k=0;k<N_Derivatives;k++)                // for all derivatives
      {                                           // get values in original cell
        OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,NeededDerivatives[k]);
        for(j=0;j<N_Points;j++)                   // for all quadrature points
        {
          Orig = OrigFEValues[j];                 // value in original cell
          value = 0;
          for(l=0;l<N_;l++)                       // for all basis functions
          {
            // accumulate value of derivative in point j
            value += FEFunctValues_test[l] * Orig[l];
          }
          Derivatives_test[m][j][k] = value;      // for k-th derivative
        }                                         // endfor j
      }                                           // endfor k
    }
    // get coefficients of pde
    if(Coeff)
    {
      Coeff(N_Points, X, Y, Param, CoeffArray);
    }
    for( int n_a = 0 ; n_a < N_Ansatz; n_a ++ )
    {
      for( int n_t = 0 ; n_t < N_Test ; n_t ++ )
      {
        // loop assembly
        switch (functional)
        {
          // L2 inner product
          case 0:
            // loop over the quadrature points
            for(j=0;j<N_Points;j++)
            {
              val = Derivatives_ansatz[n_a][j][0] * Derivatives_test[n_t][j][0];
              w = AbsDetjk[j]*weights[j];
              Assembly[ n_t ][ n_a ] += w * val;
            }
            break;
            // L2 inner product for gradients times diffusion
          case 1:
            // loop over the quadrature points
            for(j=0;j<N_Points;j++)
            {
              val = CoeffArray[j][0] * (Derivatives_ansatz[n_a][j][1] * Derivatives_test[n_t][j][1]
                + Derivatives_ansatz[n_a][j][2] * Derivatives_test[n_t][j][2]);
              w = AbsDetjk[j]*weights[j];
              Assembly[ n_t ][ n_a ] += w * val;
            }
            break;
            // integral of the function
          case 101:
            // loop over the quadrature points
            for(j=0;j<N_Points;j++)
            {
              val = Derivatives_ansatz[n_a][j][0];
              w = AbsDetjk[j]*weights[j];
              Assembly[ n_t ][ n_a ] += w * val;
            }
            break;
            // bilinear form for convection-diffusion equations
          case 1000:
            // loop over the quadrature points
            for(j=0;j<N_Points;j++)
            {
              // diffusion
              val = CoeffArray[j][0] * (Derivatives_ansatz[n_a][j][1] * Derivatives_test[n_t][j][1]
                + Derivatives_ansatz[n_a][j][2] * Derivatives_test[n_t][j][2]);
              // convection
              val += (CoeffArray[j][1] * Derivatives_ansatz[n_a][j][1] + CoeffArray[j][2] * Derivatives_ansatz[n_a][j][2]) *  Derivatives_test[n_t][j][0];
              // reaction
              val += CoeffArray[j][3] * Derivatives_ansatz[n_a][j][0] *  Derivatives_test[n_t][j][0];
              w = AbsDetjk[j]*weights[j];
              Assembly[ n_t ][ n_a ] += w * val;
            }
            break;
            // right hand side
          case 1001:
            // loop over the quadrature points
            for(j=0;j<N_Points;j++)
            {
              val = CoeffArray[j][4] * Derivatives_test[n_t][j][0];
              w = AbsDetjk[j]*weights[j];
              Assembly[ n_t ][ n_a ] += w * val;
            }
            break;

          default:
            OutPut("functional not defined !!!"<<endl);
            exit(4711);
        }
      }
    }
  }

  // set memory free
  delete[] aux3;
  delete[] aux4;
  delete[] aux1;
  delete[] aux2;
}

/********************************************************************/
//
// ComputeErrorToInterpolantL2Edges for paper with Julia Novo
//
// THIS IS JUST A CLONE OF TCD2DErrorEstimator::GetErrorEstimate !!!
//
/********************************************************************/

double ComputeErrorToInterpolantL2Edges(TFESpace2D *FESpace2D,
TFEFunction2D *FEFunction2D, DoubleFunct2D *Exact,
int N_Derivatives,
MultiIndex2D *NeededDerivatives,
CoeffFct2D *Coeff,
BoundCondFunct2D **BoundaryConds,
BoundValueFunct2D **BoundaryValues,
TAuxParam2D *Aux,
int n_fespaces,
TFESpace2D **fespaces,
double *estimated_global_error)
{
  //const int MaxN_BaseFunctions2D_loc = 16;
  const int N_BaseFuncts2D_loc = 5;
  const int MaxN_QuadPoints_1D_loc = 30;
  int i,j,k,l,n,ij,N_UsedElements, N_LocalUsedElements, found;
  int N_Cells, N_Points, N_Parameters, N_Points1D, N_Edges, N_;
  int Used[N_FEs2D], MaxN_BaseFunctions2D_loc;
  int *N_BaseFunct, BaseFunct_loc;
  BaseFunct2D *BaseFuncts;
  TFESpace2D *fespace;
  FE2D *UsedElements, LocalUsedElements[N_FEs2D], CurrentElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1D;
  BaseFunct2D BaseFunct;
  TBaseFunct2D *bf;
  TCollection *Coll;
  TBaseCell *cell, *neigh;
  BF2DRefElements bf2Drefelements;
  double *weights, *xi, *eta,*weights1D, *zeta;
  double xi1D[N_BaseFuncts2D_loc][4][MaxN_QuadPoints_1D_loc];
  double eta1D[N_BaseFuncts2D_loc][4][MaxN_QuadPoints_1D_loc];
  double *xietaval_ref1D[N_BaseFuncts2D_loc][4][MaxN_QuadPoints_1D_loc];
  double *xideriv_ref1D[N_BaseFuncts2D_loc][4][MaxN_QuadPoints_1D_loc];
  double *etaderiv_ref1D[N_BaseFuncts2D_loc][4][MaxN_QuadPoints_1D_loc];
  //  double xietaval_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D];
  //  double xideriv_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D];
  //  double etaderiv_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D];
  double *xyval_ref1D[4][MaxN_QuadPoints_1D_loc];
  double *xderiv_ref1D[4][MaxN_QuadPoints_1D_loc];
  double *yderiv_ref1D[4][MaxN_QuadPoints_1D_loc];
  double *xyval_1D[4];
  double *xderiv_1D[4];
  double *yderiv_1D[4];
  double *X1D[4], *Y1D[4], val[3], val_exact[4];
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D],*AbsDetjk1D[4];
  RefTrans2D RefTrans;
  double *Param[MaxN_QuadPoints_2D], *aux,*aux1,*aux2;
  double *Derivatives[MaxN_QuadPoints_2D];
  double *AuxArray[MaxN_QuadPoints_2D], *b1, linfb;
  int *DOF, N_DOF;
  double **OrigFEValues, *Orig, value;
  double *FEFunctValues;
  double *Values,max_loc_err;
  int *GlobalNumbers, *BeginIndex;
  double xc, yc, absdetjk1D, hE, weight_neigh, x0, x1, y0, y1, nx, ny, l2_edge_error;
  double estimated_global_errors[7], estimated_local_errors[7], cells_1D, eps=1e-6;
  int LocN_BF[N_BaseFuncts2D];
  BaseFunct2D LocBF[N_BaseFuncts2D];
  bool *SecondDer;
  TJoint *joint;
  TRefDesc *refdesc;
  TVertex *ver0,*ver1;
  const int *TmpEdVer;
  int ee_verbose=1;                               // verbosity

  int memory[3],data_base_memory;
#ifndef __MAC64__ 
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
  data_base_memory=0;
#endif
  //for (i=0;i<FEFunction2D->GetLength();i++)
  // OutPut(i << " fe " << FEFunction2D->GetValues()[i] << endl);

  // ########################################################################
  // store information in local arrays
  // ########################################################################
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  SecondDer = new bool[n_fespaces];

  memset(Used, 0, N_FEs2D*SizeOfInt);
  MaxN_BaseFunctions2D_loc = 0;
  for(i=0;i<n_fespaces;i++)
  {
    fespace = fespaces[i];                        // fe space
    n = fespace->GetN_UsedElements();             // # used finite elements
    UsedElements = fespace->GetUsedElements();    // used finite elements
    for(j=0;j<n;j++)                              // for all finite elements
    {
      CurrentElement = UsedElements[j];
      Used[CurrentElement] = 1;
      k = TFEDatabase2D::GetN_BaseFunctFromFE2D(CurrentElement);
      if (k>MaxN_BaseFunctions2D_loc)
      {
        MaxN_BaseFunctions2D_loc = k;
      }
    }                                             // enfor j
  }                                               // endfor i

  FEFunctValues = new double[MaxN_BaseFunctions2D_loc];
  for (i=0;i<N_BaseFuncts2D_loc;i++)
  {
    for (j=0;j<4;j++)
    {
      for (k=0;k<MaxN_QuadPoints_1D_loc;k++)
      {
        xietaval_ref1D[i][j][k] = new double[MaxN_BaseFunctions2D_loc];
        xideriv_ref1D[i][j][k] = new double[MaxN_BaseFunctions2D_loc];
        etaderiv_ref1D[i][j][k] = new double[MaxN_BaseFunctions2D_loc];
      }
    }
  }

  N_UsedElements = 0;                             // compute number of used elements
  for(i=0;i<N_FEs2D;i++)
    if(Used[i]) N_UsedElements++;

  if (N_UsedElements>N_BaseFuncts2D_loc)
  {
    OutPut("CD2DErrorEstimator: too many finite elements " << N_UsedElements << endl);
    OutPut("Increase N_BaseFuncts2D_loc !!!"<<endl);
    exit(4711);
  }

  UsedElements = new FE2D[N_UsedElements];        // store used finite elements
  j=0;                                            // in array
  for(i=0;i<N_FEs2D;i++)
  {
    if(Used[i])
    {
      UsedElements[j] = (FE2D)i;
      j++;
    }                                             // endif
  }

  if (ee_verbose>1)
  {
    cout << "estimator number of used elements: " << N_UsedElements << endl;
    for(i=0;i<N_UsedElements;i++)
      cout << "UsedElements[" << i << "]: " << UsedElements[i] << endl;
  }

  // ########################################################################
  // calculate values of base functions and derivatives on ref element
  // ########################################################################

  for(i=0;i<N_UsedElements;i++)                   // for used finite elements
  {
    CurrentElement = UsedElements[i];
    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
    LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(22);
    qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1D->GetFormulaData(N_Points1D, weights1D, zeta);
    if (N_Points1D > MaxN_QuadPoints_1D_loc)
    {
      OutPut("CD2DErrorEstimator: too many 1D quadrature points " << N_Points1D << endl);
      OutPut("Increase  MaxN_QuadPoints_1D_loc !!!"<<endl);
      exit(4711);
    }
    BaseFunct = BaseFuncts[CurrentElement];
    bf = TFEDatabase2D::GetBaseFunct2D(BaseFunct);// get base functions
    bf2Drefelements = bf->GetRefElement();
    for (j=0;j<N_UsedElements;j++)
    {
      if ((int) BaseFunct == (int) UsedElements[j])
        break;
    }
    BaseFunct_loc = j;

    switch(bf2Drefelements)                       // compute coordinates of line quadrature
    {                                             // points in reference cell
      // quadrilateral cell
      case BFUnitSquare :                         // edge 0
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunct_loc][0][j] = zeta[j];
          eta1D[BaseFunct_loc][0][j] = -1;
          bf->GetDerivatives(D00, zeta[j], -1, xietaval_ref1D[BaseFunct_loc][0][j]);
          bf->GetDerivatives(D10, zeta[j], -1, xideriv_ref1D[BaseFunct_loc][0][j]);
          bf->GetDerivatives(D01, zeta[j], -1, etaderiv_ref1D[BaseFunct_loc][0][j]);
        }                                         // edge 1
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunct_loc][1][j] = 1;
          eta1D[BaseFunct_loc][1][j] = zeta[j];
          bf->GetDerivatives(D00, 1, zeta[j], xietaval_ref1D[BaseFunct_loc][1][j]);
          bf->GetDerivatives(D10, 1, zeta[j], xideriv_ref1D[BaseFunct_loc][1][j]);
          bf->GetDerivatives(D01, 1, zeta[j], etaderiv_ref1D[BaseFunct_loc][1][j]);
        }                                         // edge 2
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunct_loc][2][j] = -zeta[j];
          eta1D[BaseFunct_loc][2][j] = 1;
          bf->GetDerivatives(D00, -zeta[j], 1, xietaval_ref1D[BaseFunct_loc][2][j]);
          bf->GetDerivatives(D10, -zeta[j], 1, xideriv_ref1D[BaseFunct_loc][2][j]);
          bf->GetDerivatives(D01, -zeta[j], 1, etaderiv_ref1D[BaseFunct_loc][2][j]);
        }                                         // edge 3
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunct_loc][3][j] = -1;
          eta1D[BaseFunct_loc][3][j] = -zeta[j];
          bf->GetDerivatives(D00, -1, -zeta[j], xietaval_ref1D[BaseFunct_loc][3][j]);
          bf->GetDerivatives(D10, -1, -zeta[j], xideriv_ref1D[BaseFunct_loc][3][j]);
          bf->GetDerivatives(D01, -1, -zeta[j], etaderiv_ref1D[BaseFunct_loc][3][j]);
        }
        break;

      case BFUnitTriangle :                       // triangular cell
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunct_loc][0][j] = (zeta[j]+1)/2;
          eta1D[BaseFunct_loc][0][j] = 0;
          bf->GetDerivatives(D00, (zeta[j]+1)/2, 0, xietaval_ref1D[BaseFunct_loc][0][j]);
          bf->GetDerivatives(D10, (zeta[j]+1)/2, 0, xideriv_ref1D[BaseFunct_loc][0][j]);
          bf->GetDerivatives(D01, (zeta[j]+1)/2, 0, etaderiv_ref1D[BaseFunct_loc][0][j]);
        }                                         // edge 1
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunct_loc][1][j] = (-zeta[j]+1)/2;
          eta1D[BaseFunct_loc][1][j] = (zeta[j]+1)/2;
          bf->GetDerivatives(D00, (-zeta[j]+1)/2, (zeta[j]+1)/2, xietaval_ref1D[BaseFunct_loc][1][j]);
          bf->GetDerivatives(D10, (-zeta[j]+1)/2, (zeta[j]+1)/2, xideriv_ref1D[BaseFunct_loc][1][j]);
          bf->GetDerivatives(D01, (-zeta[j]+1)/2, (zeta[j]+1)/2, etaderiv_ref1D[BaseFunct_loc][1][j]);
        }                                         // edge 2
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunct_loc][2][j] = 0;
          eta1D[BaseFunct_loc][2][j] = (-zeta[j] +1)/2;
          bf->GetDerivatives(D00, 0, (-zeta[j]+1)/2, xietaval_ref1D[BaseFunct_loc][2][j]);
          bf->GetDerivatives(D10, 0, (-zeta[j]+1)/2, xideriv_ref1D[BaseFunct_loc][2][j]);
          bf->GetDerivatives(D01, 0, (-zeta[j]+1)/2, etaderiv_ref1D[BaseFunct_loc][2][j]);
        }
        break;
    }
  }                                               // endfor i

  for (i=0;i<4;i++)                               // arrays for coordinates, values and
  {                                               // determinant for 1D quadrature
    X1D[i] = new double[N_Points1D];              // coordinates of edge i
    Y1D[i] = new double[N_Points1D];
                                                  // determinant of affine mapping
    AbsDetjk1D[i] = new double[MaxN_QuadPoints_2D];
    for (j=0;j<N_Points1D;j++)                    // arrays for values in reference cell
    {
      xyval_ref1D[i][j] = new double[MaxN_BaseFunctions2D_loc];
      xderiv_ref1D[i][j] = new double[MaxN_BaseFunctions2D_loc];
      yderiv_ref1D[i][j] = new double[MaxN_BaseFunctions2D_loc];
    }
    xyval_1D[i] = new double[N_Points1D];         // arrays for values in original cell
    xderiv_1D[i] = new double[N_Points1D];
    yderiv_1D[i] = new double[N_Points1D];
  }

  N_Parameters = Aux->GetN_Parameters();          // get number of parameters of equation
  aux = new double [MaxN_QuadPoints_2D*N_Parameters];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Param[j] = aux + j*N_Parameters;

  aux1 = new double [MaxN_QuadPoints_2D*N_Derivatives];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Derivatives[j] = aux1 + j*N_Derivatives;

  // 20 <= number of term
  aux2 = new double [MaxN_QuadPoints_2D*20];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    AuxArray[j] = aux2 + j*20;

  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  BeginIndex = FESpace2D->GetBeginIndex();

  // ########################################################################
  // prepare error estimates
  // ########################################################################

  // all spaces use same Coll
  Coll = FESpace2D->GetCollection();              // collection of mesh cells
  N_Cells = Coll->GetN_Cells();                   // number of mesh cells
  N_DOF = FEFunction2D->GetLength();              // number of global dof
  Values = FEFunction2D->GetValues();             // values of fe function

  for(i=0;i<N_Cells;i++)                          // do for all mesh cells
  {                                               // on the finest level
    cell=Coll->GetCell(i);
    k=cell->GetN_Edges();                         // # edges
    for(j=0;j<k;j++)                              // for all edges
    {
      neigh=cell->GetJoint(j)->GetNeighbour(cell);// neighbour cell
      if(neigh) neigh->SetClipBoard(-1);          // set clipboard to -1
    }
    cell->SetClipBoard(-1);
  }                                               // endfor i
  // non finest neighbours of finest cells have clipboard -1

  for(i=0;i<N_Cells;i++)                          // set clipboard of cells on finest
  {
    cell=Coll->GetCell(i);
    cell->SetClipBoard(i);
  }
  for (i=0;i<7;i++)                               // initialize some quantities
    estimated_global_errors[i]=0.0;
  max_loc_err = 0;
  // ########################################################################
  // loop over all cells
  // ########################################################################

  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
  found = 0;
  for(i=0;i<N_Cells;i++)                          // for all cells on the finest level
  {
    cell = Coll->GetCell(i);                      // next cell

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    for(j=0;j<n_fespaces;j++)
    {
      CurrentElement = fespaces[j]->GetFE2D(i,cell);
      LocalUsedElements[j] = CurrentElement;
      LocN_BF[j] = N_BaseFunct[CurrentElement];   // local basis functions
      LocBF[j] = BaseFuncts[CurrentElement];
      SecondDer[j] = true;                        // with 2nd derivative
    }
    N_LocalUsedElements = n_fespaces;

    // ####################################################################
    // calculate values on original element
    // ####################################################################

    // get reference transformation
    RefTrans = TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements,
      Coll, cell, SecondDer,
      N_Points, xi, eta, weights, X, Y, AbsDetjk);
    if(N_Parameters>0)                            // get parameters of equ.
      Aux->GetParameters(N_Points, Coll, cell, i, xi, eta, X, Y, Param);

    // calculate all needed derivatives of this FE function
    CurrentElement = FESpace2D->GetFE2D(i,cell);  // finite element on cell
    BaseFunct = BaseFuncts[CurrentElement];       // basis functions
    N_ = N_BaseFunct[CurrentElement];             // # basis functions
    DOF = GlobalNumbers + BeginIndex[i];          // dof of current mesh cell

    for(l=0;l<N_;l++)
      FEFunctValues[l] = Values[DOF[l]];          // fe values of dofs

    // compute values for all derivatives
    // in all quadrature points
    // in original mesh cell
    for(k=0;k<N_Derivatives;k++)                  // for all derivatives
    {                                             // get values in original cell
      OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,NeededDerivatives[k]);
      for(j=0;j<N_Points;j++)                     // for all quadrature points
      {
        Orig = OrigFEValues[j];                   // value in original cell
        value = 0;
        for(l=0;l<N_;l++)                         // for all basis functions
          value += FEFunctValues[l] * Orig[l];    // accumulate value of derivative in point j
        Derivatives[j][k] = value;                // for k-th derivative
      }                                           // endfor j
    }                                             // endfor k

    if(Coeff)                                     // get coefficients of pde
    {
      Coeff(N_Points, X, Y, Param, AuxArray);
    }
    // prepare 1D quadrature formula
    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
    LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(22);
    qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1D->GetFormulaData(N_Points1D, weights1D, zeta);
                                                  // update data base
    TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
      ->MakeRefElementData(LineQuadFormula);

    for (j=0;j<N_UsedElements;j++)
    {
      if ((int) BaseFunct == (int) UsedElements[j])
        break;
    }
    BaseFunct_loc = j;
    // get refinement descriptor
    refdesc=cell->GetRefDesc();
    refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer);

    N_Edges=cell->GetN_Edges();                   // # edges
    for(j=0;j<N_Edges;j++)                        // loop over all edges of cell
    {
      // get joint
      joint=cell->GetJoint(j);
      // look for neighbor
      neigh=joint->GetNeighbour(cell);
      if (neigh)
        weight_neigh = 0.5;
      else
        weight_neigh = 1.0;
      // compute length of the edge
      // get refinement descriptor
      refdesc=cell->GetRefDesc();
      refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer);
      // get vertices of face j
      ver0=cell->GetVertex(TmpEdVer[2*j]);
      ver1=cell->GetVertex(TmpEdVer[2*j+1]);
                                                  // coordinates of face j
      x0 = cell->GetVertex(TmpEdVer[2*j])->GetX();
      y0 = cell->GetVertex(TmpEdVer[2*j])->GetY();
      x1 = cell->GetVertex(TmpEdVer[2*j+1])->GetX();
      y1 = cell->GetVertex(TmpEdVer[2*j+1])->GetY();
      // only edges on some coarser level
      if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 120814)
      {
        cells_1D = 6 * pow(2.0,TDatabase::ParamDB->INTERNAL_LEVEL);

        if (((fabs(x0*cells_1D - (int)(x0*cells_1D+1e-8))<eps)&&
          (fabs(x1*cells_1D - (int)(x1*cells_1D+1e-8))<eps) &&
          fabs(x0-x1)<eps) ||
          ((fabs(y0*cells_1D - (int)(y0*cells_1D+1e-8))<eps)&&
          (fabs(y1*cells_1D - (int)(y1*cells_1D+1e-8))<eps)&&
          fabs(y0-y1)<eps))
        {
          //OutPut(found << " " << x0*cells_1D - (int)(x0*cells_1D+1e-8) << " " << x1*cells_1D - (int)(x1*cells_1D+1e-8)<<
          //    " " << y0*cells_1D - (int)(y0*cells_1D+1e-8) << " " << y1*cells_1D - (int)(y1*cells_1D+1e-8)
          //     << "::" << x0 << " " << x1 << " " << y0 << " " << y1  << endl);
          found++;
        }
        else
        {
          continue;
        }
      }
      nx = y1 - y0;                               // compute normal
      ny = x0 - x1;
      hE = sqrt(nx*nx+ny*ny);                     // length of edge
      absdetjk1D = hE/2.0;

      // get original coordinates of edge quad. points
      TFEDatabase2D::GetOrigFromRef(RefTrans,N_Points1D, xi1D[BaseFunct_loc][j],
        eta1D[BaseFunct_loc][j],
        X1D[j], Y1D[j], AbsDetjk1D[j]);

      for(k=0;k<N_Points1D;k++)                   // get values and derivatives in original cell
      {
        TFEDatabase2D::GetOrigValues(RefTrans, xi1D[BaseFunct_loc][j][k],
          eta1D[BaseFunct_loc][j][k],
          TFEDatabase2D::GetBaseFunct2D(BaseFunct),
          Coll, (TGridCell *)cell,
          xietaval_ref1D[BaseFunct_loc][j][k],
          xideriv_ref1D[BaseFunct_loc][j][k],
          etaderiv_ref1D[BaseFunct_loc][j][k],
          xyval_ref1D[j][k],
          xderiv_ref1D[j][k],
          yderiv_ref1D[j][k]);
      }

      if(Coeff)                                   // get coefficients of pde
      {
        Coeff(N_Points1D, X1D[0], Y1D[0], Param, AuxArray);
      }

      for(k=0;k<N_Points1D;k++)                   // for all quadrature points
      {
        val[0]=val[1]=val[2] = 0;
        Exact(X1D[j][k],Y1D[j][k],val_exact);
        for(l=0;l<N_;l++)                         // for all basis functions
        {
                                                  // accumulate value of derivative
          val[0] += FEFunctValues[l] * xyval_ref1D[j][k][l];
                                                  // accumulate value of derivative
          val[1] += FEFunctValues[l] * xderiv_ref1D[j][k][l];
                                                  // accumulate value of derivative
          val[2] += FEFunctValues[l] * yderiv_ref1D[j][k][l];
        }                                         // endfor l
        xyval_1D[j][k]= val[0]-val_exact[0];      // for k-th
        //xyval_1D[j][k]= val[0];                   // for k-th
        xderiv_1D[j][k]= val[1];                  // for k-th
        yderiv_1D[j][k]= val[2];                  // for k-th

        // update the integral
        //if (ee_verbose > 1)
        //{
        // OutPut(k << " " << X1D[j][k] << " " << Y1D[j][k] << " n " <<  weight_neigh << " w " <<  weights1D[k] <<
        //  " det " << absdetjk1D << " fe " << val[0]
        //  << " ex " << val_exact[0] << " " << val[0] - val_exact[0] << endl);
        //}
        // convection
        b1 = AuxArray[k]+1;
        if (fabs(b1[0])> fabs(b1[1]))
          linfb = fabs(b1[0]);
        else
          linfb = fabs(b1[1]);

        estimated_global_errors[0] += weight_neigh * absdetjk1D * weights1D[k] * linfb * xyval_1D[j][k] * xyval_1D[j][k];
      }                                           // endfor k
    }                                             // endfor j
  }                                               // endfor i (loop over cells)

  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
  data_base_memory+= memory[2]-memory[1];

  l2_edge_error = sqrt(estimated_global_errors[0]);
  // set memory free
  delete Param[0];
  delete AuxArray[0];
  delete SecondDer;
  delete UsedElements;
  for (i=0;i<4;i++)
  {
    delete X1D[i];
    delete Y1D[i];
    delete AbsDetjk1D[i];
    for (j=0;j<N_Points1D;j++)
    {
      delete xyval_ref1D[i][j];
      delete xderiv_ref1D[i][j];
      delete yderiv_ref1D[i][j];
    }
    delete xyval_1D[i];
    delete xderiv_1D[i];
    delete yderiv_1D[i];
  }
  delete aux1;

  for (i=0;i<N_BaseFuncts2D_loc;i++)
  {
    for (j=0;j<4;j++)
    {
      for (k=0;k<MaxN_QuadPoints_1D_loc;k++)
      {
        delete xietaval_ref1D[i][j][k];
        delete xideriv_ref1D[i][j][k];
        delete etaderiv_ref1D[i][j][k];
      }
    }
  }
  OutPut("found egdes " << found<<endl);
  delete FEFunctValues;
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;

  if ((memory[1]- memory[0])!=data_base_memory)
    OutPut("WARNING : Error Estimator did not set all memory free !!!" <<  memory[1]- memory[0] << endl);
  return(l2_edge_error);
}                                                 // TCD2DErrorEstimator::GetErrorEstimate



// THIS WORKS ONLY FOR P1 !!!
void JumpTermsForAdjointProblemP1(TFESpace2D *fespace,
        TFEFunction2D *u,
        CoeffFct2D *Coeffs,
        BoundCondFunct2D *BoundaryConditions,
        double *rhs)
{
  int i, j, k, ii, N_Cells, *ColInd, *RowPtr, *GlobalNumbers, *BeginIndex;
  int ActiveBound, *DOF, *DOF_n, N_Edges, boundedge, locdof, found;
  int com00, com01, com10, com11, com20, com21, com_other0, com_other1;
  int loc_vert_n, comp;
  double val[3], val_neigh[3], h, norm_t, x[3], y[3], oldval[3];
  double x_n[3], y_n[3], eps = 1e-6;
  double x0, x1, y0, y1, xs, ys, t1, t2, *coeff, jump, fac0, fac1, fac2;
  double phi0_x, phi0_y, phi1_x, phi1_y, phi2_x, phi2_y, n1, n2, maxjump;
  double phi0_n_x, phi0_n_y, phi1_n_x, phi1_n_y, phi2_n_x, phi2_n_y;
  double phi_n_other_x, phi_n_other_y, p0, p1;
  double sx, sy, tmp, meas, area, rho = 2.0, ansatz, test, area_n, meas_n;
  TBaseCell *cell, *neigh;
  TCollection *coll;
  FE2D CurrentElement;
  TJoint *joint;
  TRefDesc *refdesc;
  TVertex *ver0,*ver1;
  BoundCond BdCond;
  TBoundComp *BdComp;
  TBoundEdge *bound_edge;
  TIsoBoundEdge *isobound_edge;
  const int *TmpEdVer;

  coeff = new double[13];

  // get arrays with the numbering of the dof
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();

  // get start of dirichlet nodes in dof array
  ActiveBound = fespace->GetActiveBound();
  // get collection and number of cells
  coll = fespace->GetCollection();
  N_Cells = coll->GetN_Cells();

  // assign a numbering to the cells
  for(i=0;i<N_Cells;i++)                          // do for all mesh cells
  {                                               // on the finest level
    cell=coll->GetCell(i);
    cell->SetClipBoard(i);
  }                                               // endfor i

  // loop over all cells for computing the jump terms
  for(i=0;i<N_Cells;i++)
  {
    // next cell
    cell = coll->GetCell(i);
    h = cell->GetDiameter();
    meas = cell->GetMeasure();
    // pointer to global indices of dof connected with this cell
    DOF = GlobalNumbers + BeginIndex[i];

    // local dofs are arranged as follows
    // local dof 0 on vertex 0 opposite to edge 1
    // local dof 1 on vertex 1 opposite to edge 2
    // local dof 2 on vertex 2 opposite to edge 0

    CurrentElement = fespace->GetFE2D(i, cell);
    if (CurrentElement!=C_P1_2D_T_A)
    {
        OutPut("JumpTermsForAdjointProblem for element " << CurrentElement <<
         " not implemented !!!"<< endl);
        exit(4711);
    }
    // # of edges
    N_Edges = cell->GetN_Edges();

    sx = sy = 0;
    // compute derivatives for basis functions
    for (j=0;j<N_Edges; j++)
    {
      x[j] = cell->GetVertex(j)->GetX();
      y[j] = cell->GetVertex(j)->GetY();
      sx += x[j];
      sy += y[j];
      //OutPut(x[j] << " " << y[j] << " ");
      u->FindGradientLocal(cell, i, x[j], y[j], val);
      //OutPut("u "<<j << " " << val[1]<<endl);
    }
    sx /= N_Edges;
    sy /= N_Edges;
    //OutPut(endl);
    // compute twice area of triangle
    if (N_Edges==3)
    {
      area = 2 * meas;
      phi0_x = (y[1]-y[2])/area;
      phi0_y = (x[2]-x[1])/area;
      phi1_x = (y[2]-y[0])/area;
      phi1_y = (x[0]-x[2])/area;
      phi2_x = (y[0]-y[1])/area;
      phi2_y = (x[1]-x[0])/area;
    }
    /* OutPut("0 " << phi0_x << " " << phi0_y << endl);
     OutPut("1 " << phi1_x << " " << phi1_y << endl);
     OutPut("2 " << phi2_x << " " << phi2_y << endl);
    */
    // get refinement descriptor
    refdesc=cell->GetRefDesc();
    refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer);

    // compute gradient of current solution (constant)
    u->FindGradientLocal(cell, i, sx, sy, val);

    for(j=0;j<N_Edges;j++)                        // loop over all edges of cell
    {
      joint=cell->GetJoint(j);
      ver0=cell->GetVertex(TmpEdVer[2*j]);        // get vertices of face j
      ver1=cell->GetVertex(TmpEdVer[2*j+1]);
      x0 = ver0->GetX();                          // coordinates of face j
      y0 = ver0->GetY();
      x1 = ver1->GetX();
      y1 = ver1->GetY();
      //OutPut(endl << "ed " << j << " " << x0 << " " << y0 << " ; " << x1 << " " <<y1
      //    << endl);
      // compute tangential
      t1 = x1 - x0;
      t2 = y1 - y0;
      norm_t = sqrt(t1*t1+t2*t2);
      t1 /= norm_t;
      t2 /= norm_t;
      // compute normal
      n1 = t2;
      n2 = -t1;
      //OutPut(t1 << " " << t2 << " " << t1*t1+t2*t2 << endl);
      // compute solution (including derivative) in midpoint of tangential
      // from point of view of this mesh cell
      xs = (x1+x0)/2;
      ys = (y1+y0)/2;
      //OutPut("locd " << x0 << " " << y0 << " : " 
//       << x1 << " " << y1 << " : "  << DOF[(j)%3] << 
      //     " " << DOF[(j+1)%3] << endl);
      //u->FindGradientLocal(cell, i, xs, ys, val);
      //OutPut("grad_i " << val[1] << " " << val[2] << endl);
      // compute solution (including derivative) in midpoint of tangential
      // from point of view of neighbour mesh cell
      // NO ADAPTIVE MESHES ALLOWED
      neigh=joint->GetNeighbour(cell);            // neighbour cell
      if (neigh!=NULL)
      {
        ii =  neigh->GetClipBoard();
        //OutPut("ii " << ii << endl);
  DOF_n = GlobalNumbers + BeginIndex[ii];
        u->FindGradientLocal(neigh, ii, xs, ys, val_neigh);
        boundedge = 0;
  meas_n = neigh->GetMeasure();
  area_n = 2 * meas_n;
  // vertices of neighbour
  for (k=0;k<N_Edges; k++)
  {
      x_n[k] = neigh->GetVertex(k)->GetX();
      y_n[k] = neigh->GetVertex(k)->GetY();
  }
  
  // compute derivatives of basis fcts. in neighbour 
  // cell
  switch(j)
  {
      // edge 0: vertices 0, 1
      // local dof 0, 1
      case 0:
    // find common vertices
    for (k=0;k<N_Edges; k++)
    {
        found = 0;
        if ((fabs(x_n[k] - x[0])<eps) && (fabs(y_n[k] - y[0])<eps))
        {
      com00 = (k+1)%3;
      com01 = (k+2)%3;
      found++;
        }
        if ((fabs(x_n[k] - x[1])<eps) && (fabs(y_n[k] - y[1])<eps))
        {
      com10 = (k+1)%3;
      com11 = (k+2)%3;
      found++;
        }
        if (!found)
        {
      com_other0 = (k+1)%3;
      com_other1 = (k+2)%3;
      loc_vert_n = k;
        }
    }
    phi0_n_x = (y_n[com00]-y_n[com01])/area_n;
    phi0_n_y = (x_n[com01]-x_n[com00])/area_n;
    phi1_n_x = (y_n[com10]-y_n[com11])/area_n;
    phi1_n_y = (x_n[com11]-x_n[com10])/area_n;
    phi2_n_x = 0;
    phi2_n_y = 0;
    phi_n_other_x = (y_n[com_other0] - y_n[com_other1])/area_n;
    phi_n_other_y = (x_n[com_other1] - x_n[com_other0])/area_n;
    //OutPut(phi0_n_x << " " << phi0_n_y << " done " << endl);
    break;
      // edge 1: vertices 1, 2
      // local dof 1, 2
      case 1:
    // find common vertices
    for (k=0;k<N_Edges; k++)
    {
        found = 0;
        if ((fabs(x_n[k] - x[1])<eps) && (fabs(y_n[k] - y[1])<eps))
        {
      com10 = (k+1)%3;
      com11 = (k+2)%3;
      found++;
      //OutPut(x_n[k] << " " << y_n[k] << endl);
        }
        if ((fabs(x_n[k] - x[2])<eps) && (fabs(y_n[k] - y[2])<eps))
        {
      com20 = (k+1)%3;
      com21 = (k+2)%3;
      found++;
        }
        if (!found)
        {
      com_other0 = (k+1)%3;
      com_other1 = (k+2)%3;
      loc_vert_n = k;
        }       
    }
    phi0_n_x = 0;
    phi0_n_y = 0;
    phi1_n_x = (y_n[com10]-y_n[com11])/area_n;
    phi1_n_y = (x_n[com11]-x_n[com10])/area_n;
    phi2_n_x = (y_n[com20]-y_n[com21])/area_n;
    phi2_n_y = (x_n[com21]-x_n[com20])/area_n;
    phi_n_other_x = (y_n[com_other0] - y_n[com_other1])/area_n;
    phi_n_other_y = (x_n[com_other1] - x_n[com_other0])/area_n;
    //OutPut(phi_n_other_x << " " << phi_n_other_y << " other " << endl);
    break;
      // edge 2: vertices 2, 0
      // local dof 2, 0
      case 2:
    // find common vertices
    for (k=0;k<N_Edges; k++)
    {
        found = 0;
        if ((fabs(x_n[k] - x[0])<eps) && (fabs(y_n[k] - y[0])<eps))
        {
      com00 = (k+1)%3;
      com01 = (k+2)%3;
      found++;
      //OutPut(x_n[k] << " " << y_n[k] << endl);
        }
        if ((fabs(x_n[k] - x[2])<eps) && (fabs(y_n[k] - y[2])<eps))
        {
      com20 = (k+1)%3;
      com21 = (k+2)%3;
      found++;
        }
        if (!found)
        {
      com_other0 = (k+1)%3;
      com_other1 = (k+2)%3;
      loc_vert_n = k;
        }       
    }
    phi0_n_x = (y_n[com00]-y_n[com01])/area_n;
    phi0_n_y = (x_n[com01]-x_n[com00])/area_n;
    phi1_n_x = 0;
    phi1_n_y = 0;
    phi2_n_x = (y_n[com20]-y_n[com21])/area_n;
    phi2_n_y = (x_n[com21]-x_n[com20])/area_n;
    phi_n_other_x = (y_n[com_other0] - y_n[com_other1])/area_n;
    phi_n_other_y = (x_n[com_other1] - x_n[com_other0])/area_n;
    //OutPut(phi0_n_x << " " << phi0_n_y << " done " << endl);
    break;
  }
      }
      else
      {
        // boundary edge
  if (cell->GetJoint(j)->GetType() == BoundaryEdge)
  {
      bound_edge = (TBoundEdge *)cell->GetJoint(j);
      BdComp = bound_edge->GetBoundComp();
            bound_edge->GetParameters(p0, p1);      
  }
  if (cell->GetJoint(j)->GetType() == IsoBoundEdge)
  {
            isobound_edge = (TIsoBoundEdge *)joint;
            BdComp = isobound_edge->GetBoundComp();
            isobound_edge->GetParameters(p0, p1);
  }
  // get id of the boundary component
  comp=BdComp->GetID();
  BoundaryConditions(comp, (p0+p1)/2.0, BdCond);
  if (BdCond==NEUMANN)
  {
      boundedge = 2;
      // continuation of solution
      val_neigh[0] = val_neigh[1] = val_neigh[2] = 0;
      // continuation of test functions
      phi0_n_x = phi0_n_y = phi1_n_x = phi1_n_y = phi2_n_x = phi2_n_y;
  }
  if (BdCond==DIRICHLET)
  {
      boundedge = 1;
  }
      }
      //OutPut("grad_ii " << val_neigh[1] << " " << val_neigh[2] << endl);
      // do nothing for Dirichlet edges
      if (boundedge==1)
    continue;
      // compute ansatz factor
      fac0 = val[1] * n1 + val[2] * n2;
      fac1 = val_neigh[1] * n1 + val_neigh[2] * n2;
      //OutPut("jump " << val[1] - val_neigh[1] << " " << val[2] - val_neigh[2] << endl);
      ansatz = fac0 - fac1;
      //OutPut("("<<x0<<"," << y0<<") (" << x1 << "," << y1 << "):  normal jump " << ansatz << endl);
      //OutPut(" a "  << ansatz << " " );
      ansatz *= 2;
      Coeffs(1, &xs, &ys, NULL, &coeff);
      ansatz *= coeff[0]*sqrt(coeff[0]);      
      // length of edge
      ansatz *= norm_t; 
      // weight of error estimator
      fac0 = norm_t/sqrt(coeff[0]);
      if (TDatabase::ParamDB->INTERNAL_COERCIVITY>0)
      {
    if (1.0/sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY)<fac0)
        fac0 = 1.0/sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY); // update weight 
      }
      //OutPut(" fac0 "  << fac0 << " ");
      ansatz *= fac0;
      for (k=0;k<N_Edges; k++)
      {
    // compute jumps of derivative of test function
    switch(k)
    {
        case 0:
      fac0 = phi0_x * n1 + phi0_y * n2;
      fac1 = phi0_n_x * n1 + phi0_n_y * n2;
      break;
        case 1:
      fac0 = phi1_x * n1 + phi1_y * n2;
      fac1 = phi1_n_x * n1 + phi1_n_y * n2;
      break;
        case 2:
      fac0 = phi2_x * n1 + phi2_y * n2;
      fac1 = phi2_n_x * n1 + phi2_n_y * n2;
      break;
    }
    test = fac0 - fac1;
    //OutPut("("<<x0<<"," << y0<<") (" << x1 << "," << y1 << "):  testjump " << test << endl);
    test *= ansatz;
    //OutPut(k << " t " << test << " " );
    
    // update the rhs
    switch(k)
    {
        // dof zero 
        case 0:
      // local dof 0
      locdof = DOF[0];
      break;
      // dof one
        case 1:
      // local dof 1
      locdof = DOF[1];
      break;
      // dof two
        case 2:
      // local dof 2
      locdof = DOF[2];
      break;
    }
    //OutPut(locdof << " integral ("<<x0<<"," << y0<<") (" << x1 << "," << y1 << "):  integral " << test << endl);
    rhs[locdof] += test;
    //OutPut("locd " << locdof << " " <<  rhs[locdof] << endl);
      } // end k (test functions)
      if (boundedge)
    continue;
      // the test functions opposite the edge
      fac1 = phi_n_other_x * n1 + phi_n_other_y * n2;
      test = - fac1;

      test *= ansatz;
      // compute the dof
      locdof = DOF_n[loc_vert_n];
      //OutPut(locdof << " integral ("<<x0<<"," << y0<<") (" << x1 << "," << y1 << "):  integral_bc " 
      //     << test << "::"<< -phi_n_other_x  << " " << -phi_n_other_y <<endl);
      rhs[locdof] += test;
      //OutPut("locd " << locdof << " " <<  rhs[locdof] << endl);
    } // end j (edges)
  }                                               // loop over cells
  /*
    // loop over all cells for computing the edge stabilization
    for(i=0;i<N_Cells;i++)
    {
      // next cell
      cell = coll->GetCell(i);
      // pointer to global indices of dof connected with this cell
      DOF = GlobalNumbers + BeginIndex[i];

      // local dofs are arranged as follows
      // local dof 0 on vertex 0 opposite to edge 1
  // local dof 1 on vertex 1 opposite to edge 2
  // local dof 2 on vertex 2 opposite to edge 0

  // # of edges
  N_Edges = cell->GetN_Edges();

  // compute derivatives for basis functions
  for (j=0;j<N_Edges; j++)
  {
  x[j] = cell->GetVertex(j)->GetX();
  y[j] = cell->GetVertex(j)->GetY();
  OutPut(x[j] << " " << y[j] << " ori " << rhsori[DOF[j]] << " update " << rhs[DOF[j]] <<
  " diff " <<  rhsori[DOF[j]] - rhs[DOF[j]] << endl);
  }
  }
  */
  delete coeff;

}

// =======================================================================
//
// JumpTermsForAdjointProblem
//
// computes jump terms in rhs of adjoint problem for energy norm 
// error estimator
//
// TFESpace2D *fespace                   -- finite element space
// TFEFunction2D *u                      -- finite element function
// CoeffFct2D *Coeff                     -- coefficients of the equation
// BoundCondFunct2D *BoundaryConditions  -- pointer to function for boundary
//                                          conditions
// double *rhs                           -- array for right hand side
//
// =======================================================================

void JumpTermsForAdjointProblem(TFESpace2D *fespace,
        TFEFunction2D *u,
        CoeffFct2D *Coeff,
        BoundCondFunct2D *BoundaryConditions,
        double *rhs)
 {
  const int MaxN_BaseFunctions2D_Ersatz = 100;

  double hK,w,integrand,edge_par,sigma_par,diffusion;
  int out;
  int i,j,k,l,l1,l2,l3,n,n_neigh,m,r,q,dummy,N_UsedElements,N_LocalUsedElements,ii,jj,ll;
  int N_Cells, N_Points, N_Parameters, N_Points1D, N_Edges, N_, N_Hanging;
  int N_Test, N_Ansatz, N_Joints, n_rhs = 1, n_fespaces = 1;
  int Used[N_FEs2D];
  int *N_BaseFunct;
  BaseFunct2D *BaseFuncts;
  TBaseFunct2D *bf;
  FE2D *UsedElements, LocalUsedElements[N_FEs2D], CurrentElement;
  FE2D TestElement, AnsatzElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1D;
  BaseFunct2D BaseFunctCell;
  TCollection *Coll;
  TBaseCell *cell;
  TJoint *joint;
  TBoundEdge *boundedge;
  TIsoBoundEdge *isoboundedge;
  int **RhsGlobalNumbers, **RhsBeginIndex;
  int **TestGlobalNumbers, **TestBeginIndex;
  int **AnsatzGlobalNumbers, **AnsatzBeginIndex;
  TFE2D *ele;
  TFEDesc2D *FEDesc_Obj;
  BF2DRefElements bf2Drefelements;
  double *weights, *xi, *eta, *weights1D, *weights_neigh, *xi_neigh, *eta_neigh, *weights1D_neigh;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], X_neigh[MaxN_QuadPoints_2D], Y_neigh[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D], AbsDetjk_neigh[MaxN_QuadPoints_2D],*AbsDetjk1D[4];
  double *Param[MaxN_QuadPoints_2D];
  double *local_rhs;
  double *righthand;
  double **Matrices, *aux, *aux2, *aux3, *aux4;
  double **Matrix;
  double ***LocMatrices, **LocRhs;
  int LocN_BF[N_BaseFuncts2D];
  BaseFunct2D LocBF[N_BaseFuncts2D];
  double *Coeffs[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last;
  int *TestDOF, *AnsatzDOF;
  double *Entries,*Entries1,*Entries2,*Entries3, *Entries4, *Entries5;
  int *ColInd, *RowPtr;
  int *ColInd1, *RowPtr1,*ColInd2, *RowPtr2, *ColInd3, *RowPtr3;
  int *ColInd4, *RowPtr4, *ColInd5, *RowPtr5;
  double *RHS, *MatrixRow;
  double **HangingEntries, **HangingRhs;
  double *CurrentHangingEntries, *CurrentHangingRhs;
  int *HangingRowPtr, *HangingColInd;
  THangingNode *hn, **HangingNodes;
  HNDesc HNDescr;
  THNDesc *HNDescr_Obj;
  double *Coupling, v;
  TBoundComp *BoundComp;
  double t0, t1, t, s, integral, u_val[4], u_val_neigh[4];
  double *u_jump_x, *u_jump_y;
  int comp, dof_ii,dof_jj, found;
  BoundCond Cond0, Cond1;
  BoundCondFunct2D *BoundaryCondition;
  BoundValueFunct2D *BoundaryValue;
  TNodalFunctional2D *nf;
  int N_EdgePoints;
  double *EdgePoints;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_BaseFunctions2D_Ersatz];
  int *EdgeDOF, N_EdgeDOF;
  int N_LinePoints;
  double *LineWeights, *zeta;
  double x0, x1, y0, y1, hE, nx, ny, tx, ty, x, y, val, eps=1e-12;
  double penetration_penalty, friction_parameter;
  double **JointValues, *JointValue, u1_values[3], u2_values[3];
  double delta;
  bool *SecondDer;

  double *Coefficients1D[MaxN_QuadPoints_2D];
  double *Parameters1D[MaxN_QuadPoints_2D];

  double xi1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D], eta1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D];
  //double xietaval_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];
  //double xideriv_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];
  //double etaderiv_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];
  double**** xietaval_ref1D = new double*** [N_BaseFuncts2D];
  double**** xideriv_ref1D = new double*** [N_BaseFuncts2D];
  double**** etaderiv_ref1D = new double*** [N_BaseFuncts2D];
  double *xyval_ref1D[4][MaxN_QuadPoints_1D];
  double *xderiv_ref1D[4][MaxN_QuadPoints_1D];
  double *yderiv_ref1D[4][MaxN_QuadPoints_1D];
  double *X1D[4], *Y1D[4], *X1D_neigh[4], *Y1D_neigh[4];
  RefTrans2D RefTrans;
  int N_DOF;
  double *Values;
  //double value_basefunct_ref1D[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz],value_basefunct_ori[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz];
  //double xderiv_basefunct_ref1D[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz],xderiv_basefunct_ori[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz];
  //double yderiv_basefunct_ref1D[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz],yderiv_basefunct_ori[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz];
  double*** value_basefunct_ref1D = new double** [N_BaseFuncts2D];
  double*** xderiv_basefunct_ref1D = new double** [N_BaseFuncts2D];
  double*** yderiv_basefunct_ref1D = new double** [N_BaseFuncts2D];
  double *value_basefunct_ori[6];
  double *xderiv_basefunct_ori[6];
  double *yderiv_basefunct_ori[6];
  double x_pos_ref[6];
  double y_pos_ref[6];
  double x_pos[6];
  double y_pos[6];
  double *value_basefunct_ori_neigh[6];
  double *xderiv_basefunct_ori_neigh[6];
  double *yderiv_basefunct_ori_neigh[6];
  double x_pos_neigh[6];
  double y_pos_neigh[6];
  double dummy2[6];

  int neigh_edge, ref_n;
  int neigh_N_,N_Neigh;
  double absdet1D_neigh[MaxN_QuadPoints_2D];
  double xi1DNeigh[N_BaseFuncts2D][MaxN_QuadPoints_1D], eta1DNeigh[N_BaseFuncts2D][MaxN_QuadPoints_1D];
  double *X1DNeigh,*Y1DNeigh;
  TBaseCell *neigh;
  FE2D LocalUsedElements_neigh[N_FEs2D], CurrEleNeigh;
  BaseFunct2D BaseFunctNeigh;
  QuadFormula2D QuadFormulaNeigh;
  TQuadFormula2D *qfNeigh;
  QuadFormula1D LineQuadFormulaNeigh;
  TQuadFormula1D *qf1DNeigh;
  int LocN_BF_neigh[N_BaseFuncts2D];
  BaseFunct2D LocBF_neigh[N_BaseFuncts2D];
  int N_Points1DNeigh,N_PointsNeigh;
  double *weights1DNeigh,*zetaNeigh,*weightsNeigh,*xiNeigh,*etaNeigh;
  TFE2D *eleNeigh;
  RefTrans2D RefTransNeigh;
  BF2DRefElements bf2DrefelementsNeigh;
  int *DOF_neigh;
  double xietaval_refNeigh1D[N_BaseFuncts2D][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];
  double xideriv_refNeigh1D[N_BaseFuncts2D][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];
  double etaderiv_refNeigh1D[N_BaseFuncts2D][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];
  double *xyval_refNeigh1D[MaxN_QuadPoints_1D];
  double *xderiv_refNeigh1D[MaxN_QuadPoints_1D];
  double *yderiv_refNeigh1D[MaxN_QuadPoints_1D];
  double *xderiv_Neigh1D, *yderiv_Neigh1D, *xyval_Neigh1D;

  double jump_xyval[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_xderiv[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_yderiv[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];

  out=1;
  OutPut("JumpTerms"<<endl);

  // ########################################################################
  // store information in local arrays
  // ########################################################################
  // get basis functions from data base
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  // get number of basis functions from data base
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  // for the right hand side (n_rhs == 1)
  if(n_rhs)
  {
      // get information on the global degrees of freedom 
      RhsBeginIndex = new int* [n_rhs];
      RhsGlobalNumbers = new int* [n_rhs];
      for(i=0;i<n_rhs;i++)
      {
    RhsBeginIndex[i] = fespace->GetBeginIndex();
    RhsGlobalNumbers[i] = fespace->GetGlobalNumbers();
      }                                             // endfor
      // allocate arrays for the local computation of the right hand side
      LocRhs = new double* [n_rhs];
      righthand = new double [n_rhs*MaxN_BaseFunctions2D];
      for(i=0;i<n_rhs;i++)
    LocRhs[i] = righthand+i*MaxN_BaseFunctions2D;
  }                                               // endif n_rhs
  
  // no second derivatives necessary
  SecondDer = new bool[1];
  SecondDer[0] = FALSE;

  // for all possible basis functions
  // allocate arrays for values on the actual mesh cells
  for (i=0;i<N_BaseFuncts2D;i++)
  {
    value_basefunct_ref1D[i] = new double* [6];
    xderiv_basefunct_ref1D[i] = new double* [6];
    yderiv_basefunct_ref1D[i] = new double* [6];
    for (j=0;j<6;j++)
    {

      value_basefunct_ref1D[i][j] = new double [MaxN_BaseFunctions2D_Ersatz];
      xderiv_basefunct_ref1D[i][j] = new double [MaxN_BaseFunctions2D_Ersatz];
      yderiv_basefunct_ref1D[i][j] = new double [MaxN_BaseFunctions2D_Ersatz];

      memset( value_basefunct_ref1D[i][j] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
      memset( xderiv_basefunct_ref1D[i][j] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
      memset( yderiv_basefunct_ref1D[i][j] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );

    }
  }

  // for all possible basis functions
  // allocate arrays for values on the reference mesh cell
  for (i=0;i<N_BaseFuncts2D;i++)
  {
    xietaval_ref1D[i] = new double** [4];
    xideriv_ref1D[i] = new double** [4];
    etaderiv_ref1D[i] = new double** [4];
    for (j=0;j<4;j++)
    {
      xietaval_ref1D[i][j] = new double* [MaxN_QuadPoints_1D];
      xideriv_ref1D[i][j] = new double* [MaxN_QuadPoints_1D];
      etaderiv_ref1D[i][j] = new double* [MaxN_QuadPoints_1D];
      for (n=0;n<MaxN_QuadPoints_1D;n++)
      {
        xietaval_ref1D[i][j][n] = new double [MaxN_BaseFunctions2D_Ersatz];
        xideriv_ref1D[i][j][n] = new double [MaxN_BaseFunctions2D_Ersatz];
        etaderiv_ref1D[i][j][n] = new double [MaxN_BaseFunctions2D_Ersatz];

        memset( xietaval_ref1D[i][j][n] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
        memset( xideriv_ref1D[i][j][n] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
        memset( etaderiv_ref1D[i][j][n] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
      }
    }
  }

  memset(Used, 0, N_FEs2D*SizeOfInt);

  // for the finite element space
  for(i=0;i<n_fespaces;i++)
  {
    n = fespace->GetN_UsedElements();             /* # used finite elements */
    UsedElements = fespace->GetUsedElements();    /* used finite elements */
    for(j=0;j<n;j++)                              /* for all finite elements */
    {
      CurrentElement = UsedElements[j];
      Used[CurrentElement] = 1;
    }                                             // enfor j
  }                                               // endfor i

  N_UsedElements = 0;                             /* compute number of used elements */
  for(i=0;i<N_FEs2D;i++)
    if(Used[i]) N_UsedElements++;

  UsedElements = new FE2D[N_UsedElements];        /* store used finite elements */
  j=0;                                            /* in array */
  for(i=0;i<N_FEs2D;i++)
    if(Used[i])
  {
    UsedElements[j] = (FE2D)i;
    if (out==2)
  OutPut("element " << UsedElements[j] << endl);
    j++;
  }                                               // endif

  // ########################################################################
  // calculate values of base functions and derivatives on ref element
  // ########################################################################
  if (out==2)
  {
      OutPut("N_UsedElements: " << N_UsedElements << " " << endl); fflush(0);
      OutPut("N_BaseFuncts2D: " << N_BaseFuncts2D << " " << endl);
      OutPut("MaxN_QuadPoints_1D: " << MaxN_QuadPoints_1D << " " << endl);
      OutPut("MaxN_BaseFunctions2D_Ersatz: " << MaxN_BaseFunctions2D_Ersatz << " " << endl);
  }

  for(n=0;n<N_UsedElements;n++)                   // for used finite elements
  {
    CurrentElement = UsedElements[n];
    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
    LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1D->GetFormulaData(N_Points1D, weights1D, zeta);
    BaseFunctCell = BaseFuncts[CurrentElement];
                                                  // get base functions
    bf = TFEDatabase2D::GetBaseFunct2D(BaseFunctCell);
    bf2Drefelements = bf->GetRefElement();
    switch(bf2Drefelements)                       // compute coordinates of line quadrature
    {                                             // points in reference cell
      // quadrilateral cell
      case BFUnitSquare :                         // edge 0

        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][0][j] = zeta[j];
          eta1D[BaseFunctCell][0][j] = -1;
          bf->GetDerivatives(D00, zeta[j], -1, xietaval_ref1D[BaseFunctCell][0][j]);
          bf->GetDerivatives(D10, zeta[j], -1, xideriv_ref1D[BaseFunctCell][0][j]);
          bf->GetDerivatives(D01, zeta[j], -1, etaderiv_ref1D[BaseFunctCell][0][j]);
        }                                         // edge 1
        x_pos_ref[0] = -1;
        y_pos_ref[0] = 1;
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][1][j] = 1;
          eta1D[BaseFunctCell][1][j] = zeta[j];
          bf->GetDerivatives(D00, 1, zeta[j], xietaval_ref1D[BaseFunctCell][1][j]);
          bf->GetDerivatives(D10, 1, zeta[j], xideriv_ref1D[BaseFunctCell][1][j]);
          bf->GetDerivatives(D01, 1, zeta[j], etaderiv_ref1D[BaseFunctCell][1][j]);
        }                                         // edge 2
        x_pos_ref[1] = 1;
        y_pos_ref[1] =-1;
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][2][j] = -zeta[j];
          eta1D[BaseFunctCell][2][j] = 1;
          bf->GetDerivatives(D00, -zeta[j], 1, xietaval_ref1D[BaseFunctCell][2][j]);
          bf->GetDerivatives(D10, -zeta[j], 1, xideriv_ref1D[BaseFunctCell][2][j]);
          bf->GetDerivatives(D01, -zeta[j], 1, etaderiv_ref1D[BaseFunctCell][2][j]);
        }                                         // edge 3
        x_pos_ref[2] = 1;
        y_pos_ref[2] = 1;
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][3][j] = -1;
          eta1D[BaseFunctCell][3][j] = -zeta[j];
          bf->GetDerivatives(D00, -1, -zeta[j], xietaval_ref1D[BaseFunctCell][3][j]);
          bf->GetDerivatives(D10, -1, -zeta[j], xideriv_ref1D[BaseFunctCell][3][j]);
          bf->GetDerivatives(D01, -1, -zeta[j], etaderiv_ref1D[BaseFunctCell][3][j]);
        }
        x_pos_ref[3] = -1;
        y_pos_ref[3] = -1;
        
        ref_n=4;
        break;

      case BFUnitTriangle :                       // triangular cell

        bf->GetDerivatives(D00, 0, 0, value_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D10, 0, 0, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 0, 0, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[0] = 0;
        y_pos_ref[0] = 0;
        bf->GetDerivatives(D00, 1, 0, value_basefunct_ref1D[BaseFunctCell][1]);
        bf->GetDerivatives(D10, 1, 0, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 1, 0, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[1] = 1;
        y_pos_ref[1] = 0;
        bf->GetDerivatives(D00, 0, 1, value_basefunct_ref1D[BaseFunctCell][2]);
        bf->GetDerivatives(D10, 0, 1, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 0, 1, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[2] = 0;
        y_pos_ref[2] = 1;

        bf->GetDerivatives(D00, 0.5, 0, value_basefunct_ref1D[BaseFunctCell][3]);
        bf->GetDerivatives(D10, 0.5, 0, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 0.5, 0, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[3] = 0.5;
        y_pos_ref[3] = 0;
        bf->GetDerivatives(D00, 0.5, 0.5, value_basefunct_ref1D[BaseFunctCell][4]);
        bf->GetDerivatives(D10, 0.5, 0.5, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 0.5, 0.5, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[4] = 0.5;
        y_pos_ref[4] = 0.5;
        bf->GetDerivatives(D00, 0, 0.5, value_basefunct_ref1D[BaseFunctCell][5]);
        bf->GetDerivatives(D10, 0, 0.5, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 0, 0.5, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[5] = 0;
        y_pos_ref[5] = 0.5;
        
        ref_n=6;
        for (j=0;j<N_Points1D;j++)                // for all quadrature poin
        {
          xi1D[BaseFunctCell][0][j] = (zeta[j]+1)/2;
          eta1D[BaseFunctCell][0][j] = 0;
          bf->GetDerivatives(D00, (zeta[j]+1)/2, 0, xietaval_ref1D[BaseFunctCell][0][j]);
          bf->GetDerivatives(D10, (zeta[j]+1)/2, 0, xideriv_ref1D[BaseFunctCell][0][j]);
          bf->GetDerivatives(D01, (zeta[j]+1)/2, 0, etaderiv_ref1D[BaseFunctCell][0][j]);
        }                                         // edge 1
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][1][j] = (-zeta[j]+1)/2;
          eta1D[BaseFunctCell][1][j] = (zeta[j]+1)/2;
          bf->GetDerivatives(D00, (-zeta[j]+1)/2, (zeta[j]+1)/2, xietaval_ref1D[BaseFunctCell][1][j]);
          bf->GetDerivatives(D10, (-zeta[j]+1)/2, (zeta[j]+1)/2, xideriv_ref1D[BaseFunctCell][1][j]);
          bf->GetDerivatives(D01, (-zeta[j]+1)/2, (zeta[j]+1)/2, etaderiv_ref1D[BaseFunctCell][1][j]);
        }                                         // edge 2
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][2][j] = 0;
          eta1D[BaseFunctCell][2][j] = (-zeta[j] +1)/2;
          bf->GetDerivatives(D00, 0, (-zeta[j]+1)/2, xietaval_ref1D[BaseFunctCell][2][j]);
          bf->GetDerivatives(D10, 0, (-zeta[j]+1)/2, xideriv_ref1D[BaseFunctCell][2][j]);
          bf->GetDerivatives(D01, 0, (-zeta[j]+1)/2, etaderiv_ref1D[BaseFunctCell][2][j]);
        }
        break;
    }
  }                                               // endfor n
  if (out==2)
      OutPut("basefunct" << endl);

  for(l=0;l<ref_n;l++)
  {
    value_basefunct_ori[l] = new double[MaxN_BaseFunctions2D_Ersatz];
    xderiv_basefunct_ori[l]  = new double[MaxN_BaseFunctions2D_Ersatz];
    yderiv_basefunct_ori[l]  = new double[MaxN_BaseFunctions2D_Ersatz];
    value_basefunct_ori_neigh[l] = new double[MaxN_BaseFunctions2D_Ersatz];
    xderiv_basefunct_ori_neigh[l]  = new double[MaxN_BaseFunctions2D_Ersatz];
    yderiv_basefunct_ori_neigh[l]  = new double[MaxN_BaseFunctions2D_Ersatz];
  }

  for(m=0;m<4;m++)                                // arrays for coordinates, values and
  {                                               // determinant for 1D quadrature
    X1D[m] = new double[N_Points1D];              // coordinates of edge i
    Y1D[m] = new double[N_Points1D];
                                                  // determinant of affine mapping
    AbsDetjk1D[m] = new double[MaxN_QuadPoints_2D];
    for (j=0;j<N_Points1D;j++)                    // arrays for values in reference cell
    {
      xyval_ref1D[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];
      xderiv_ref1D[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];
      yderiv_ref1D[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];
    }
  }                                               // endfor m
  // arrays for jump of the gradient
  u_jump_x = new double[N_Points1D]; 
  u_jump_y = new double[N_Points1D]; 

  for (j=0;j<N_Points1D;j++)                      // arrays for values in reference cell
  {
    xyval_refNeigh1D[j] = new double[MaxN_BaseFunctions2D_Ersatz];
    xderiv_refNeigh1D[j] = new double[MaxN_BaseFunctions2D_Ersatz];
    yderiv_refNeigh1D[j] = new double[MaxN_BaseFunctions2D_Ersatz];
  }

  // ########################################################################
  // Arrays for Parameters
  // ########################################################################

  if (out==2)
      OutPut("coeff" << endl);
  // 20 <= number of term
  aux2 = new double [MaxN_QuadPoints_2D*20];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Coeffs[j] = aux2 + j*20;

  aux4 = new double [MaxN_QuadPoints_2D*20];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Coefficients1D[j] = aux4 + j*20;

  // ########################################################################
  // prepare loop over cells
  // ########################################################################

  // all spaces use same Coll
  Coll = fespace->GetCollection();            // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)                          // set clipboard of cells on finest
  {
    cell=Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  // ########################################################################
  // loop over all cells
  // ########################################################################
  for(i=0;i<N_Cells;i++)                          // for all cells on the finest level
  {
    cell = Coll->GetCell(i);                      // next cell

    if (out==2)
  OutPut("cell " << i << endl);

  // calculate all needed derivatives of this FE function
  CurrentElement = fespace->GetFE2D(i,cell);  // finite element on cell

  BaseFunctCell = BaseFuncts[CurrentElement]; // basis functions
  N_ = N_BaseFunct[CurrentElement];           // # basis functions
  DOF = RhsGlobalNumbers[0] + RhsBeginIndex[0][i];  // dof of current mesh cell
  if (out==2)
  {
      OutPut("N_ basefct " << N_ << " dof: ");
      for (j=0;j<N_;j++)
    OutPut(DOF[j] << " ");
      OutPut(endl);
  }
  LocalUsedElements[0] = CurrentElement;
  LocN_BF[0] = N_BaseFunct[CurrentElement];   // local basis functions
  LocBF[0] = BaseFuncts[CurrentElement];
  SecondDer[0] = FALSE;
  RefTrans = TFEDatabase2D::GetOrig(1, LocalUsedElements,
            Coll, cell, SecondDer,
            N_Points, xi, eta, weights, X, Y, AbsDetjk);

  // get coefficients of pde
  if(Coeff) Coeff(N_Points, X, Y, Param, Coeffs);
  // prepare 1D quadrature formula
  l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
  if(out==2)
      OutPut("Polynomial degree on cell: " << l << endl);
  LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
  qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
  qf1D->GetFormulaData(N_Points1D, weights1D, zeta);

  if(out==2)
  {
      for(j=0;j<N_Points1D; j++)
      {
          OutPut("weights1D["<<j<<"]:" <<  weights1D[j] << endl);
      }
      OutPut(endl);
  }

  // update data base
  TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
      ->MakeRefElementData(LineQuadFormula);
  N_Edges=cell->GetN_Edges();                 // # edges
  if (out==2)
  {
      for(r=0;r<N_Edges;r++)
      {
    cell->GetVertex(r)->GetCoords(x0, y0);
    cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1);
    OutPut("Local edge r: " << r << " vertex A " << x0 << " " << y0 << " vertex B " << 
     x1 << " " << y1 << endl);
      }
  }
 

  for(r=0;r<N_Edges;r++)                      // loop over all edges of cell
  {                                           // get original coordinates of edge quad. points
      TFEDatabase2D::GetOrigFromRef(RefTrans,N_Points1D, xi1D[BaseFunctCell][r],
            eta1D[BaseFunctCell][r],
            X1D[r], Y1D[r], AbsDetjk1D[r]);
      
      for(j=0;j<N_Points1D;j++)                 // get values and derivatives in original cell
      {
          TFEDatabase2D::GetOrigValues(RefTrans, xi1D[BaseFunctCell][r][j],
               eta1D[BaseFunctCell][r][j],
               TFEDatabase2D::GetBaseFunct2D(BaseFunctCell),
               Coll, (TGridCell *)cell,
               xietaval_ref1D[BaseFunctCell][r][j],
               xideriv_ref1D[BaseFunctCell][r][j],
               etaderiv_ref1D[BaseFunctCell][r][j],
               xyval_ref1D[r][j],
               xderiv_ref1D[r][j],
               yderiv_ref1D[r][j]);
      }
      
  }                                           // endfor r
  
  TFEDatabase2D::GetOrigFromRef(RefTrans,ref_n,x_pos_ref,y_pos_ref,x_pos,y_pos,dummy2);
  for(l=0;l<ref_n;l++)
  {
      TFEDatabase2D::GetOrigValues(RefTrans, x_pos_ref[l],
           y_pos_ref[l],
           TFEDatabase2D::GetBaseFunct2D(BaseFunctCell),
           Coll, (TGridCell *)cell,
           value_basefunct_ref1D[BaseFunctCell][l],
           xderiv_basefunct_ref1D[BaseFunctCell][l],
           yderiv_basefunct_ref1D[BaseFunctCell][l],
           value_basefunct_ori[l],
           xderiv_basefunct_ori[l],
           yderiv_basefunct_ori[l]);
      // OutPut("Hallo: x_pos_ref[l]: " << x_pos_ref[l] << "value_basefunct_ref1D[BaseFunctCell][l]: " << value_basefunct_ref1D[BaseFunctCell][l] << endl);
  }

  for(r=0;r<N_Edges;r++)
  {                                           
      // for each edge, get the corresponding neighbour cell.
      neigh=cell->GetJoint(r)->GetNeighbour(cell);
      //#######################################################################//
      // get coefficients on edges
      //only implemented for coeffs that do not depend on the params
      //#######################################################################//
      
      //  if(N_Parameters>0)                // get parameters of equ.
      // Parameters->GetParameters(N_Points1D, Coll, cell, i, xi1D[BaseFunctCell][r], eta1D[BaseFunctCell][r], X1D[r], Y1D[r], Param1D);
      
      if(Coeff) Coeff(N_Points1D, X1D[r], Y1D[r], Param, Coefficients1D);
      //#######################################################################//
      // If there is a neighbour to the edge, do...
      if(neigh)
      {                                         
          // get the number of this neigbbour cell from the clipboard
          q = neigh->GetClipBoard();
    if (out==2)
        OutPut("neighbor " << q << endl);
    // only neighbors with larger number
    // to ensure that each inner edge is treated only once
          if(1)//i<q)
          {
        // calculate all needed derivatives of this FE function
        // finite element on neighbour
        CurrEleNeigh = fespace->GetFE2D(q,neigh);
        BaseFunctNeigh = BaseFuncts[CurrEleNeigh];
        eleNeigh =  TFEDatabase2D::GetFE2D(CurrEleNeigh);
        //BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();    // basis functions on neighbour
        N_Neigh = eleNeigh->GetN_DOF();       // number of basis functions on neighbour
        // dof of current mesh cell on neighbour cell
        DOF_neigh = RhsGlobalNumbers[0] + RhsBeginIndex[0][q];
        
        if (out==2)
        {
      OutPut("neigh N_ basefct " << N_Neigh << " dof: ");
      for (j=0;j<N_;j++)
          OutPut(DOF_neigh[j] << " ");
      OutPut(endl);
        }
        LocalUsedElements_neigh[0] = CurrEleNeigh;
        // local basis functions
        LocN_BF_neigh[0] = N_BaseFunct[CurrEleNeigh];
        LocBF_neigh[0] = BaseFuncts[CurrEleNeigh];
        
        RefTransNeigh = TFEDatabase2D::GetOrig(1, LocalUsedElements_neigh,
                 Coll, neigh, SecondDer,
                 N_Points, xi_neigh, eta_neigh, 
                 weights_neigh, X_neigh, Y_neigh, AbsDetjk_neigh);
        
        // get edge of the neighbour cell which is the edge r
        neigh_edge=0;
        while(neigh->GetJoint(neigh_edge)->GetNeighbour(neigh)!=cell) neigh_edge ++;
        if (out==2)
      OutPut("neigh_edge " << neigh_edge << endl);

        // arrays for coordinates on neighbour cell
        for(m=0;m<4;m++)                      
        {
      X1D_neigh[m] = new double[N_Points1D];
      Y1D_neigh[m] = new double[N_Points1D];
        }
        
        // get original coordinates of edge quad. points of neighbour cell
        TFEDatabase2D::GetOrigFromRef(RefTransNeigh,N_Points1D, xi1D[BaseFunctNeigh][neigh_edge],
              eta1D[BaseFunctNeigh][neigh_edge],
              X1D_neigh[neigh_edge], Y1D_neigh[neigh_edge], 
              AbsDetjk1D[neigh_edge]);
        
        // get values and derivatives on original neighbour cell on edge neigh_edge
        for (j=0;j<N_Points1D;j++)
        {                                     
      if(out==2)
      {  
          OutPut("X1D[r][j]: " << X1D[r][j] << " Y1D[r][j]: "  <<  Y1D[r][j] <<  
           " X1D[neigh_edge][j] " << X1D_neigh[neigh_edge][j] << 
           " Y1D[neigh_edge][j]: " <<  Y1D_neigh[neigh_edge][j] <<  endl);
      }
        }
        
        if(X1D_neigh[neigh_edge][0] == X1D[r][0] && Y1D_neigh[neigh_edge][0] == Y1D[r][0] )
        {
      if(out==2)
      {
          OutPut("Quadrature points on neighbour edge in the correct order." << endl);
      }
      for (j=0;j<N_Points1D;j++)
      {
          TFEDatabase2D::GetOrigValues(RefTransNeigh, xi1D[BaseFunctNeigh][neigh_edge][j],
               eta1D[BaseFunctNeigh][neigh_edge][j],
               TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
               Coll, (TGridCell *)neigh,
               xietaval_ref1D[BaseFunctNeigh][neigh_edge][j],
               xideriv_ref1D[BaseFunctNeigh][neigh_edge][j],
               etaderiv_ref1D[BaseFunctNeigh][neigh_edge][j],
               xyval_refNeigh1D[j],
               xderiv_refNeigh1D[j],
               yderiv_refNeigh1D[j]);
      }                                   //endfor j
        }                                     //endif
        else
        {
      if(out==2)
      {
          OutPut("Inverse the order of the quadrature points on neighbour edge !" << endl);
      }
      for (j=0;j<N_Points1D;j++)
      {
          TFEDatabase2D::GetOrigValues(RefTransNeigh, xi1D[BaseFunctNeigh][neigh_edge][j],
               eta1D[BaseFunctNeigh][neigh_edge][j],
               TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
               Coll, (TGridCell *)neigh,
               xietaval_ref1D[BaseFunctNeigh][neigh_edge][N_Points1D-j-1],
               xideriv_ref1D[BaseFunctNeigh][neigh_edge][N_Points1D-j-1],
               etaderiv_ref1D[BaseFunctNeigh][neigh_edge][N_Points1D-j-1],
               xyval_refNeigh1D[j],
               xderiv_refNeigh1D[j],
               yderiv_refNeigh1D[j]);
      }                                   //endfor j
        }                                     //endelse
        
        TFEDatabase2D::GetOrigFromRef(RefTransNeigh,ref_n,x_pos_ref,y_pos_ref,
              x_pos_neigh,y_pos_neigh,dummy2);

        for(l=0;l<ref_n;l++)
        {
      TFEDatabase2D::GetOrigValues(RefTrans, x_pos_ref[l],
                 y_pos_ref[l],
                 TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
                 Coll, (TGridCell *)neigh,
                 value_basefunct_ref1D[BaseFunctNeigh][l],
                 xderiv_basefunct_ref1D[BaseFunctNeigh][l],
                 yderiv_basefunct_ref1D[BaseFunctNeigh][l],
                 value_basefunct_ori_neigh[l],
                 xderiv_basefunct_ori_neigh[l],
                 yderiv_basefunct_ori_neigh[l]);
        }

        if (out>2)
        {
      for(k=0;k<N_; k++)
      {
          for(l=0;l<ref_n;l++)
          {
        OutPut("basis fct: "<< DOF[k] << 
         " (x,y)-coordinate: " << x_pos[l] << " " << y_pos[l] << 
         " value: " <<  value_basefunct_ori[l][k] << endl);
          }
      }
        }
        
        for(k=0;k<N_Neigh; k++)
        {
      if (out>2)
      {
          for(l=0;l<ref_n;l++)
          {
        if(out>2)
        {
            OutPut("Basisfkt neigh: "<< DOF_neigh[k] <<"  (x,y)-coordinate: " << 
             x_pos_neigh[l] << " " << y_pos_neigh[l] << " value: " <<  
             value_basefunct_ori_neigh[l][k] << endl);
        }
          }
      }
        }                                     //endfor k
      
        // compute the jumps of the basis functions
        // and of their derivatives in the quadrature points on edge r
        // first for the basis functions of cell i
        
        for(k=0;k<N_; k++)
        {
      dummy = 0;
      l=0;
      // Check if basis function k of cell i is in the FE-Space of neighbour cell q
      while(l<N_Neigh && dummy == 0)
      {
          if(DOF[k] == DOF_neigh[l])dummy=1;
          l++;
      }
      l = l-1;
      // if basis function k of cell i is in the local FE-Space of neighbour cell q do
      if(dummy ==1 )
      {   
                      // assumption: N_Points1D cell =  N_Points1D neighbour !!!!
          for(j=0;j<N_Points1D;j++)
          {
        jump_xyval[j][k] = xyval_ref1D[r][j][k]  -  xyval_refNeigh1D[j][l];
        jump_xderiv[j][k] = xderiv_ref1D[r][j][k] - xderiv_refNeigh1D[j][l];
        jump_yderiv[j][k] = yderiv_ref1D[r][j][k] - yderiv_refNeigh1D[j][l];
        
        if(out>2)
        {
            OutPut("xyval_cell: "<< xyval_ref1D[r][j][k] <<" xyval_Neigh: " << xyval_refNeigh1D[j][l] << " jump= " <<  jump_xyval[j][k] << " of basefunction: "<< DOF[k] << " in cell: " << i << " on edge (local): " << r << " in quadrature point: " << j << endl);
            OutPut("xderiv_cell: "<< xderiv_ref1D[r][j][k] <<  " xderiv_Neigh: " << xderiv_refNeigh1D[j][l] << " jump= " << jump_xderiv[j][k] << " of basefunction: "<< DOF[k] << " in cell: " << i <<  " on edge: " << r << " in quadrature point: " << j << endl);
            OutPut("yderiv_cell: "<< yderiv_ref1D[r][j][k] <<  " yderiv_Neigh: " << yderiv_refNeigh1D[j][l] << " jump= " << jump_yderiv[j][k] << " of basefunction: "<< DOF[k]  << " in cell: " << i <<  " on edge: " << r << " in quadrature point: " << j << endl);
            OutPut(endl);
        }
          }
      }
      //endif

      // if basis function k of cell i is NOT in the local FE-Space of neighbour cell q 
      // extend them by zero to cell q
      if (dummy == 0)
      {
          for(j=0;j<N_Points1D;j++)
          {
        jump_xyval[j][k]  = xyval_ref1D[r][j][k] ;
        jump_xderiv[j][k] = xderiv_ref1D[r][j][k];
        jump_yderiv[j][k] = yderiv_ref1D[r][j][k];
        
        if(out>2)
        {
            OutPut("No Neighbour: xyval_cell: "<< xyval_ref1D[r][j][k] << " jump= " <<  jump_xyval[j][k] <<  " of basefunction: "<< DOF[k] << " in cell: " << i << " on edge (local): " << r << " in quadrature point: " << j << endl);
            OutPut("No Neighbour: x-deriv-jump= " <<  jump_xderiv[j][k] << " of basefunction: "<< DOF[k] << " in cell: " << i <<  " on edge: " << r << " in quadrature point: " << j << endl);
            OutPut("No Neighbour: y-deriv-jump= " <<  jump_yderiv[j][k] << " of basefunction: "<< DOF[k] << " in cell: " << i <<  " on edge: " << r << " in quadrature point: " << j << "\n" << endl);
            OutPut(endl);
        }
          }                                 //endfor j
      }                                   //endif
        }                                     //endfor k
        
        // now for the basis functions of neighbour cell q
        // which are not in the local FE-Space of cell i
        for(l=0;l<N_Neigh; l++)
        {
      dummy = 0;
      k=0;
      while(k<N_ && dummy == 0 )
      {
          if(DOF_neigh[l] == DOF[k]) dummy=1 ;
          k++;
      }
      k=k-1;

      // if basis function l of neighbour cell q is NOT in the local FE-Space of cell i do
      // extend them by zero to cell i
      if( dummy == 0)
      {         
          for(j=0;j<N_Points1D;j++)
          {
        jump_xyval[j][l+N_] = -xyval_refNeigh1D[j][l] ;
        jump_xderiv[j][l+N_]= -xderiv_refNeigh1D[j][l];
        jump_yderiv[j][l+N_]= -yderiv_refNeigh1D[j][l];
        if(out>2)
        {
            OutPut("Neighbour!!" << "xyval_Neigh: " << xyval_refNeigh1D[j][l] << " jump= " <<  jump_xyval[j][l+N_]<<  " of basefunction: "<< DOF_neigh[l] << " in cell: " << q << " on edge (local): " << r << " in quadrature point: " << j << endl);
            OutPut("Neighbour!! " << "x-deriv-jump: " << jump_xderiv[j][l+N_] << " of basefunction: "<< DOF_neigh[l] << " in cell: " << q <<  " on edge: " << r << " in quadrature point: " << j << endl);
            OutPut("Neighbour!! y-deriv-jump= " <<  jump_yderiv[j][l+N_]<< " of basefunction: "<< DOF_neigh[l] << " in cell: " << q <<  " on edge: " << r << " in quadrature point: " << j << "\n" << endl);
            OutPut(endl);
        }
          }                                 //endfor j
      }                                   //endif
        }                                     //endfor l
        
        // compute jumps of gradient of current solution in the
        // quadrature points on the edge
        for(j=0;j<N_Points1D;j++)
        {
      // quadrature point s
      if (out==2)
          OutPut("quad point " << X1D_neigh[neigh_edge][j] << 
           " " << Y1D_neigh[neigh_edge][j] << endl);
      // compute values on cell i
      u->FindGradientLocal(cell,i,X1D_neigh[neigh_edge][j],Y1D_neigh[neigh_edge][j],
               u_val);
      // compute values on cell q
      u->FindGradientLocal(neigh,q,X1D_neigh[neigh_edge][j],Y1D_neigh[neigh_edge][j],
               u_val_neigh);
      // compute jumps (value of cell i - value of q)
      // same direction as for basis functions
      u_jump_x[j] = u_val[1] - u_val_neigh[1];
      u_jump_y[j] = u_val[2] - u_val_neigh[2];
      
      if (out==2)
          OutPut("jump of gradient " << u_jump_x[j] << " " << u_jump_y[j] << endl);
        }
        // #################################################################################
        // Compute the edge integrals with the jumps of the basis functions and their derivatives
        // #################################################################################
        
        // get vertices of the edge
        // the edge is counter-clockwise wrt cell i
        cell->GetVertex(r)->GetCoords(x0, y0);
        cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1);
        
        if(out==2)
      OutPut("vertex A " << x0 << " " << y0 << " vertex B " << x1 << " " << y1 << endl);
        // compute length of the edge
        hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
        // compute normal vector to this edge (normalized)
        // this normal is pointing outward of cell i (an inward of cell q)
        nx = (y1-y0)/hE;
        ny = (x0-x1)/hE;
        // tangential normal vector to this edge (normalized)
        // directed from x0 to x1
        tx = (x1-x0)/hE;
        ty = (y1-y0)/hE;
        
        // compute weight of the jump term
        sigma_par = TDatabase::ParamDB->INTERNAL_COERCIVITY;
        diffusion = Coeffs[0][0];
        edge_par = hE/sqrt(diffusion);
        if (sigma_par > 0)
        {
      sigma_par = 1/sqrt(sigma_par);
      if (sigma_par < edge_par)
          edge_par = sigma_par;
        }
        edge_par = 2 * edge_par/sqrt(diffusion);
        if (out==2){ OutPut("weight of jump term: " << edge_par << " " << Coeffs[0][0] << endl)};
        
        ActiveBound = fespace->GetActiveBound();
        if(out>2)
        {
      OutPut("ActiveBound "<< ActiveBound  << endl);
      for(j=0;j<N_Points1D; j++)
      {
          OutPut("Det["<<j<<"]:" <<  AbsDetjk1D[r][j] <<" b1: "<<Coefficients1D[j][1] <<" b2: " << Coefficients1D[j][2] << endl);
      }
        }
    
        // edge integral: test function from cell
        if(out==2)
        {  
      OutPut("testfct. cell" << endl);
        }
        for (ii=0;ii<N_;ii++)                 //ii - test function of cell
        {
      // look for 'ii'-th row in all matrices
      dof_ii = DOF[ii];
      // Dirichlet node
      if (dof_ii>=ActiveBound)
          continue;
      
      // initialize edge integral
      integral=0;
      // compute edge integral
      for (j=0;j<N_Points1D;j++)       
      {
          if (out==2)
          {
        OutPut("(" << x0 << "," << y0 << ") (" << x1 << "," << y1 <<"): normal jump " << u_jump_x[j]*nx + u_jump_y[j]*ny << endl);
          OutPut("(" << x0 << "," << y0 << ") (" << x1 << "," << y1 <<"): testjump " << jump_xderiv[j][ii]*nx + jump_yderiv[j][ii]*ny << " : " << jump_xderiv[j][ii]*nx << " " <<  jump_yderiv[j][ii]*ny << endl);
          }
          integrand = edge_par * diffusion * diffusion 
        * (u_jump_x[j]*nx + u_jump_y[j]*ny) *
        (jump_xderiv[j][ii]*nx + jump_yderiv[j][ii]*ny);
          // hE/2 : determinant of reference transformation
          w = weights1D[j]*hE/2;
          integral += w*integrand;        // integral on the edge
      }
      if (out==2)
          OutPut(dof_ii << " integral " << "(" << x0 << "," << y0 << ") (" << x1 << "," << y1 <<"): "
           << integral << endl);
      // update rhs
      rhs[dof_ii] += integral;
        }                                     // end outer loop over dof (ii)
        
        // edge integral: test function from neigh cell q
        // test functions that belong also to cell i are already treated
        if(out==2)
        {  
      OutPut("testfct. neigh cell" << endl);
        }
        for (ii=0;ii<N_;ii++)                 //ii - test function of cell
        {
      // look for 'ii'-th row in all matrices
      dof_ii = DOF_neigh[ii];
      // Dirichlet node
      if (dof_ii>=ActiveBound)
          continue;
      // check if test function belongs also to mesh cell i
      dummy = 0;
      k=0;
      while(k<N_ && dummy == 0)
      {
          if(dof_ii == DOF[k])
        dummy=1;
          k++;
      }
      if (dummy) 
          continue;
      // initialize edge integral
      integral=0;
      // compute edge integral
      for (j=0;j<N_Points1D;j++)       
      {
          if (out==2)
          {
        OutPut("(" << x0 << "," << y0 << ") (" << x1 << "," << y1 <<"): normal jump " << u_jump_x[j]*nx + u_jump_y[j]*ny << endl);
          OutPut("(" << x0 << "," << y0 << ") (" << x1 << "," << y1 <<"): testjump " << jump_xderiv[j][N_+ii]*nx + jump_yderiv[j][N_+ii]*ny << " : " << jump_xderiv[j][N_+ii]*nx << " " <<  jump_yderiv[j][N_+ii]*ny << endl);
          }
          integrand = edge_par * diffusion * diffusion 
        *(u_jump_x[j]*nx + u_jump_y[j]*ny) *
        (jump_xderiv[j][N_+ii]*nx + jump_yderiv[j][N_ + ii]*ny);
          // hE/2 : determinant of reference transformation
          w = weights1D[j]*hE/2;
          integral += w*integrand;        // integral on the edge
      }
      if (out==2)
          OutPut(dof_ii <<" integral b" << "(" << x0 << "," << y0 << ") (" << x1 << "," << y1 <<"): "
           << integral << " :: " << jump_xderiv[0][N_ + ii]  << " " << jump_yderiv[0][N_ + ii] << endl);
      // update rhs
      rhs[dof_ii] += integral;
        }                                     // end outer loop over dof (ii)
        for (m=0;m<4;m++)
        {
      delete X1D_neigh[m];
      delete Y1D_neigh[m];
        }                                     //endfor m
    } // end if (i<q)
      } // end if(neigh)
  } // end loop over the edges (r)
  } // end loop over the cells (i)

  if (out==2)
      OutPut("free memory " << endl);

  delete SecondDer;

  delete u_jump_x;
  delete u_jump_y;

  delete UsedElements;
  for (i=0;i<4;i++)
  {
    delete X1D[i];
    delete Y1D[i];
    delete AbsDetjk1D[i];
    for (j=0;j<N_Points1D;j++)
    {
      delete xyval_ref1D[i][j];
      delete xderiv_ref1D[i][j];
      delete yderiv_ref1D[i][j];
    }
  }
  for(l=0;l<ref_n;l++)
  {
    delete value_basefunct_ori[l];
    delete xderiv_basefunct_ori[l];
    delete yderiv_basefunct_ori[l];
    delete value_basefunct_ori_neigh[l];
    delete xderiv_basefunct_ori_neigh[l];
    delete yderiv_basefunct_ori_neigh[l];
  }

  for (i=0;i<N_BaseFuncts2D;i++)
  {
     for (j=0;j<ref_n;j++)
    {
      delete value_basefunct_ref1D[i][j];
      delete xderiv_basefunct_ref1D[i][j];
      delete yderiv_basefunct_ref1D[i][j];
    }
     delete value_basefunct_ref1D[i];
     delete xderiv_basefunct_ref1D[i];
     delete yderiv_basefunct_ref1D[i];
  }


  for (i=0;i<N_Points1D;i++)
  {
    delete xyval_refNeigh1D[i];
    delete  xderiv_refNeigh1D[i];
    delete   yderiv_refNeigh1D[i];
  }

  delete aux2;
  delete aux4;

  if(n_rhs)
  {
    delete righthand;
    delete LocRhs;
    delete RhsBeginIndex;
    delete RhsGlobalNumbers;
  }
  for(i=0; i < N_BaseFuncts2D; i++)
  {
    for(j=0; j < 4; j++)
      {
  for(m=0; m < MaxN_QuadPoints_1D; m++)
          {
      //for(n=0; n < MaxN_BaseFunctions2D_Ersatz; n++)
      //{
      delete xietaval_ref1D[i][j][m];
      delete xideriv_ref1D[i][j][m];
      delete etaderiv_ref1D[i][j][m];
      // }
    }
  delete xietaval_ref1D[i][j];
  delete xideriv_ref1D[i][j];
  delete etaderiv_ref1D[i][j];
      }
    delete xietaval_ref1D[i];
    delete xideriv_ref1D[i];
    delete  etaderiv_ref1D[i];
  }

  delete value_basefunct_ref1D;
  delete xderiv_basefunct_ref1D;
  delete yderiv_basefunct_ref1D;

  delete xietaval_ref1D;
  delete xideriv_ref1D;
  delete etaderiv_ref1D;

  int N_Rows;
  // ####################################################################
  // print jump term
  // ####################################################################
  if(out==2)
  {
    N_Rows = fespace->GetN_DegreesOfFreedom();
    for(i=0;i<N_Rows;i++)
      OutPut(setw(5) << i << setw(20) << rhs[i] << endl);
  }

}                                                 // end of JumpTermsForAdjointProblem


// ========================================================================
// parameters:  partial derivatives
// ========================================================================
void Params_All_Deriv(double *in, double *out)
{
  out[0] = in[2];                                 // u
  out[1] = in[3];                                 // u_x
  out[2] = in[4];                                 // u_y
  out[3] = in[5];                                 // u_xx
  out[4] = in[6];                                 // u_yy
}

// ========================================================================
// parameters:  partial derivatives and velocity field
// ========================================================================
void Params_All_Deriv_And_Velo(double *in, double *out)
{
  out[0] = in[2];                                 // u
  out[1] = in[3];                                 // u_x
  out[2] = in[4];                                 // u_y
  out[3] = in[5];                                 // u_xx
  out[4] = in[6];                                 // u_yy
  out[5] = in[7];                                 // v_1
  out[6] = in[8];                                 // v_2
}

// ========================================================================
// parameters: two fe values + pw constant projection
// ========================================================================
void Params_Sol_And_Pw_Const_Proj(double *in, double *out)
{
  out[0] = in[2];                                 // u_x
  out[1] = in[3];                                 // u_y
  out[2] = in[4];                                 // first projection
  out[3] = in[5];                                 // second projectin
}


// ========================================================================
// parameters: two fe values + pw constant projection + velocity field
// ========================================================================
void Params_Sol_And_Pw_Const_Proj_And_Velo(double *in, double *out)
{
  out[0] = in[2];                                 // u_x
  out[1] = in[3];                                 // u_y
  out[2] = in[4];                                 // first projection
  out[3] = in[5];                                 // second projectin
}

