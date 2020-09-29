// =======================================================================
// @(#)MainRoutines.h        
//
// Purpose: contains routines which are called from the main program 
//
// Author: Volker John 
//
// History: start of implementation 22.09.2009
//
// =======================================================================

#ifndef __MAINROUTINES2D__
#define __MAINROUTINES2D__
#include <list>
#include <ConvDiff2D.h>

void SetParametersCDAdapt2D();

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
		    double* &neum_to_diri_param);

void AllocateArraysCD_2D(int LEVELS, TFEFunction2D **&UArray, 
                         TFEVectFunct2D **&VelocityArray, TFESpace2D **&USpaces,
                         TSquareMatrix2D **&MatricesA, TMultiGrid2D *&MG,
                         double **&RhsArray,
                         double *&l2, double *&h1, double *&energy,
                         double *&sd, double *&l_inf, double *&lp,
                         double *Parameters, int *&N_Uarray,
                         int &mg_type, int &mg_level);

void RefineGridCD_2D(TDomain *Domain, TCollection *coll,
         int LEVELS, int curr_level, 
         double *&eta_K, double &maximal_local_error,
         double *estimated_global_error,
         int current_estimator);
         
void GetCollectionCD_2D(TDomain *Domain, TCollection *&coll, 
      int i, int level, int mg_level, int &N_Cells);

void PrepareDataOutputCD_2D(
    TOutput2D *Output, 
    TDomain *Domain, 
    TFEFunction2D *u, 
    TFEVectFunct2D *pw_const_proj_fefct);

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
    double **RhsArray);
      
void AllocateMatricesCD_2D(TFESpace2D *fe_space,TSquareStructure2D *&sqstructureA,
    TSquareMatrix2D *&sqmatrixA, TSquareMatrix2D **MatricesA,
    int k, double *&matrix_D_Entries);

void PrepareSolverCD_2D(TSquareMatrix2D *sqmatrixA, TMultiGrid2D *MG, TMGLevel2D *&MGLevel,
    int level, int mg_type, int mg_level, int &low,
    double *sol, double *rhs);

void PrepareOutPutCD_2D(TDomain *Domain, TFEFunction2D *conc, TFEVectFunct2D *velocity,
    TFEFunction2D *conc_adjoint,
    TFEVectFunct2D *pw_const_proj_fefct, TFEVectFunct2D *pw_const_param_fe,
    TFEVectFunct2D *pw_const_param_deriv_fe, TOutput2D *&Output,
    int level, int mg_level);

void GetInformationFromCoarseGridCD_2D(TFESpace2D **USpaces, TFESpace2D *old_u_space,
    TFEFunction2D **UArray, TFEFunction2D *old_u,
    int level, int mg_type, int mg_level, int FirstSolve,
    double *sol, double *oldsol);

void MeasureErrorsCD_2D(TFESpace2D *fe_space, TFEFunction2D *conc,
    DoubleFunct2D *Exact,
    BoundCondFunct2D **BoundaryConditions,
    BoundValueFunct2D **BoundaryValues,
    CoeffFct2D *BilinearCoeffs, int level, int mg_level, int FirstSolve,
    int *N_Uarray, double *sol, double *rhs,
    double *l2, double *h1, double *sd, double *energy, double *l_inf, double *lp);

/******************************************************************************/
// measures the error to the solution on the currently finest grid
/****************F**************************************************************/

void MeasureErrorsToFinestGridSolutionCD_2D(TFESpace2D **fe_spaces, TFEFunction2D **concs,
    DoubleFunct2D *Exact,
    BoundCondFunct2D **BoundaryConditions,
    BoundValueFunct2D **BoundaryValues,
    CoeffFct2D *BilinearCoeffs, int level, int mg_level, int FirstSolve,
    int *N_Uarray, double *sol, double *rhs,
    double *l2, double *h1, double *sd, double *energy, double *l_inf, double *lp,
    double **&sol_array, int &first_comp);

/******************************************************************************/
// restricts the solution on the currently finest grid and saves it
/****************F**************************************************************/

void SaveRestrictedFinestGridSolutionCD_2D(TFEFunction2D **concs, int mg_level,
    int FirstSolve);

/******************************************************************************/
// reads the restricted solution from the finest grid that was saved
/****************F**************************************************************/

void ReadRestrictedFinestGridSolutionCD_2D(TFEFunction2D *concs);


void Assemble_CD_2D(TDomain *domain,
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
        int problem_id, TFEFunction2D **param_funct);


void Solver(TSquareMatrix **sqmatrices, TMatrix **matrices,
	    double *rhs, double *sol, 
	    MatVecProc *MatVect, DefectProc *Defect,
	    TMultiGrid2D *MG, 
	    int N_Unknowns, int ns_type);

void OutputData2D(std::ostringstream& os, TOutput2D *Output, int counter);

void ComputeErrorEstimate(TCollection *coll, TFEFunction2D *u,
			  CoeffFct2D *Coeffs, BoundCondFunct2D **BoundaryConditions,
			  BoundValueFunct2D **BoundaryValues,
			  double* &eta_K,
			  double *maximal_local_error,
			  double *estimated_global_error,
			  double l2, double h1, double cd, int N_Unknowns);

void CheckWrongNeumannNodesUnitSquareDiri(TFESpace2D *fespace,
            int &N_neum_to_diri, int* &neum_to_diri,
            int* &neum_to_diri_bdry, 
            double* &neum_to_diri_param);
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
               double* &defect, double * rhs_edge);


void ComputeAdjointMatrix(TSquareMatrix2D *A,TSquareMatrix2D *AT);


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
			     int N_U, int N_Cells);

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
			   double *neum_to_diri_param);

void ComputeDerivativeOfEstimator2D(TCollection *coll,
				    CoeffFct2D *Coeffs,
				    TFEFunction2D *velo,
				    TFEFunction2D *u_adjoint,
				    double *pw_const_param_deriv);

void  RestrictParametersToAdmissibleValues(TCollection *coll,
					   TFESpace2D *pw_const_param_space,
                                           CoeffFct2D *Coeffs,
					   double* pw_const_param);

void CheckWrongNeumannNodesUnitSquareDiri(TCollection *Coll, TFESpace2D *fespace,
					  int &N_neum_to_diri, int* &neum_to_diri,
					  int* &neum_to_diri_bdry, 
					  double* &neum_to_diri_param);

void CheckWrongNeumannNodesUnitSquareDiri(TFESpace2D *fespace,
            int &N_neum_to_diri, int* &neum_to_diri,
            int* &neum_to_diri_bdry, 
            double* &neum_to_diri_param);
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
               double* &defect, double * rhs_edge);

void AssembleMatrixforBilinearForm_SameCells(double** Assembly, 
                                             TCollection *coll, 
                                             TFEFunction2D **ansatz,
                                             TFEFunction2D **test,
                                             int N_Ansatz, int N_Test,
                                             CoeffFct2D *Coeff, int functional);

void AssembleMatricesTCD_2D(double*** Matrices, int N_Matrices, double* rhs,
                            double* l2_test_bnd, double* blf_test_bnd,
                            TCollection *coll, TFEFunction2D **ansatz,
                            TFEFunction2D **test,
                            TFEFunction2D* boundary_ansatz, int N_Ansatz,
                            int N_Test, CoeffFct2D *Coeff, int functional);

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
double *estimated_global_error);

void JumpTermsForAdjointProblemP1(TFESpace2D *fespace, TFEFunction2D *u,
                                  CoeffFct2D *Coeffs, 
                                  BoundCondFunct2D *BoundaryConditions,
                                  double *rhs);

void JumpTermsForAdjointProblem(TFESpace2D *fespace, TFEFunction2D *u,
                                CoeffFct2D *Coeffs,
                                BoundCondFunct2D *BoundaryConditions,
                                double *rhs);


/** @todo we really need to change the way the discrete forms are created, this 
 *        is very poor programming style */


					  /** all derivatives of current N_FESpaces_All_Deriv as params */
void Params_All_Deriv(double *in, double *out);

int N_FESpaces_All_Deriv = 1;
int N_Fct_All_Deriv = 1;
int N_ParamFct_All_Deriv = 1;
int N_FEValues_All_Deriv = 5;
int N_Params_All_Deriv = 5;
int FEFctIndex_All_Deriv[5] = { 0, 0, 0, 0, 0 };
ParamFct *Fct_All_Deriv[1] = { Params_All_Deriv };
int BeginParam_All_Deriv[1] = { 0 };

/** all derivatives of current solution as params */
void Params_Sol(double *in, double *out);

int N_FESpaces_Sol = 1;
int N_Fct_Sol = 1;
int N_ParamFct_Sol = 1;
int N_FEValues_Sol = 1;
int N_Params_Sol = 1;
int FEFctIndex_Sol[1] = { 0 };
ParamFct *Fct_Sol[1] = { Params_Sol };
int BeginParam_Sol[1] = { 0 };

// parameters:  partial derivatives
void Params_All_Deriv(double *in, double *out);

// parameters:  partial derivatives and velocity field
void Params_All_Deriv_And_Velo(double *in, double *out);

void Params_Sol_And_Pw_Const_Proj(double *in, double *out);

// parameters: two fe values + pw constant projection + velocity field
void Params_Sol_And_Pw_Const_Proj_And_Velo(double *in, double *out);

#endif // __MAINROUTINES2D__
