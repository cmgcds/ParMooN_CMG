// =======================================================================
// @(#)ConvDiff2D_Adjoint.h
//
// Purpose: contains routines which are called for parameter optimization
//          with adjoint problem
//
// Author: Volker John
//
// History: start of implementation 15.08.2011
//
// =======================================================================

#ifndef __CONVDIFF2D_ADJOINT__
#define __CONVDIFF2D_ADJOINT__

/******************************************************************************/
//
// ROUTINES FOR A POSTERIORI COMPUTATION OF STABILIZATION PARAMETERS
//
/******************************************************************************/

void ComputeAdjointMatrix(TSquareMatrix2D *A,TSquareMatrix2D *AT);

void PrepareAdjointProblem2D(TCollection *coll,
                             TFEFunction2D* &u_adjoint,
                             TSquareMatrix2D* &sqmatrixAadjoint,
                             TFESpace2D* &pw_const_param_space,
                             TFEVectFunct2D* &pw_const_param_fe,
                             TFEVectFunct2D* &pw_const_param_deriv_fe,
                             TFESpace2D *velocity_space,
                             TSquareStructure2D *sqstructureA,
                             BoundCondFunct2D *BoundCondition,
                             CoeffFct2D *Coeffs,
                             double* &sol_adjoint,
                             double* &rhs_edge,
                             double* &pw_const_param,
                             double* &pw_const_param_deriv,
                             double* &pw_const_param_deriv_old,
                             double* &pw_const_param_old,
                             double* &pw_const_param_zero,
                             double* &search_direction, 
                             double* &search_direction_old,
                             double* &functional_adjoint, 
                             double* &weight_local_adjoint, 
                             double** &bfgs_s,
                             double** &bfgs_y,
                             int* &indicator_type_adjoint,
                             int N_U, int N_Cells);

/******************************************************************************/
// SolveAdjointProblem2D
//  - assembles right hand side for adjoint problem
//  - solves adjoint problem
//  - computes derivatives of error estimator wrt stabilization parameters
/******************************************************************************/
void SolveAdjointProblem2D(TCollection *coll,
                           TDiscreteForm2D *DiscreteForm,
                           TFESpace2D *conc_space,
                           TFEFunction2D *conc,
                           TFEVectFunct2D *velocity_array,
                           TFEFunction2D *conc_adjoint,
                           TSquareMatrix2D* &sqmatrixAadjoint,
                           CoeffFct2D *Coeffs,
                           BoundCondFunct2D **BoundaryConditions,
                           BoundValueFunct2D **BoundaryValues,
                           double *rhs, double *rhs_edge,
                           double *sol_adjoint,
                           double *pw_const_param,
                           double *pw_const_param_zero,
                           double *pw_const_param_deriv,
                           int N_U, int N_Active, int N_neum_to_diri,
                           int *neum_to_diri, int *neum_to_diri_bdry,
                           double *neum_to_diri_param, int N_Cells);

void ComputeDerivativeOfResidual_wrt_Parameter(TCollection *coll,
                                               CoeffFct2D *Coefficients,
                                               TFEFunction2D *u, 
                                               TFEVectFunct2D *velocity_array,
                                               TFEFunction2D *u_adjoint,
                                               double *pw_const_param_deriv);

/******************************************************************************/
// restrict parameters to admitted values
/******************************************************************************/

void  RestrictParametersToAdmissibleValues(TCollection *coll,
                                           TFESpace2D *pw_const_param_space,
                                           CoeffFct2D *Coeffs,
                                           double* pw_const_param);

/******************************************************************************/
// computes the value of the target functional
/******************************************************************************/
double ComputeValueOfTargetFunctional(TCollection *coll, TFEFunction2D *u,
                                      TFEVectFunct2D *velocity_array,
                                      CoeffFct2D *Coefficients, 
                                      DoubleFunct2D *Exact,
                                      BoundCondFunct2D **BoundaryConditions,
                                      BoundValueFunct2D **BoundaryValues,
                                      int N_Unknowns, int N_Cells,
                                      int *indicator_type_adjoint,
                                      double *pw_const_param,
                                      double *pw_const_param_zero,
                                      double *functional_local_adjoint,
                                      double *weight_local_adjoint,
                                      int loop_count, int inner_loop_count);

/******************************************************************************/
// checks if a stopping criterion is fulfilled
/******************************************************************************/
int CheckStopOfIteration(int &loop_count, int &inner_loop_count, int &solve_adj,
                         int N_adjoint_funct_values, int N_Cells,
                         double observed_functional, 
                         double *adjoint_funct_values, double *pw_const_param,
                         double *pw_const_param_zero, 
                         double *pw_const_param_deriv,
                         double *functional_local_adjoint, 
                         double *weight_local_adjoint, 
                         int *indicator_type_adjoint,
                         int &stage_adjoint, int &bfgs_count);

/******************************************************************************/
// prepares the second stage of the optimization
/******************************************************************************/
void PrepareSecondStageOfOptimization(int N_Cells, double observed_functional, 
                                      double *functional_local_adjoint,
                                      double *pw_const_param,
                                      int *indicator_type_adjoint,
                                      int &bfgs_count);

//******************************************************************************/
// remove mesh cells from optimization 
// routine for a posteriori computation of stabilization parameters
/******************************************************************************/
void RemoveMeshCellsFromOptimization(TCollection *Coll, int N_Cells,
                                     double observed_functional, 
                                     double *functional_local_adjoint,
                                     double *pw_const_param,
                                     double *pw_const_param_zero,
                                     double *weight_local_adjoint,
                                     int *indicator_type_adjoint,
                                     int &bfgs_count);

/******************************************************************************/
// change the weights of the local contributions to the target functional
/******************************************************************************/
void ChangeWeights(int N_Cells, double observed_functional, 
                   double *functional_local_adjoint,
                   double *weight_local_adjoint, int &bfgs_count);

/******************************************************************************/
// computes the value of the penalty term
/******************************************************************************/
double ComputeValueOfPenalty(int N_, double *pw_const_param, 
                             double *pw_const_param_zero);

/******************************************************************************/
// computes the value of the penalty term
/******************************************************************************/
void ComputeDerivativeOfPenalty(int N_, double *pw_const_param, 
                                double *pw_const_param_zero,
                                double *pw_const_param_deriv);
/******************************************************************************/
// computes the total variation and similar expressions
/******************************************************************************/
double ComputeTotalVariation(TCollection *coll,
                             TFEFunction2D *u, TFEVectFunct2D *velocity_array,
                             int type, int *indicator_type_adjoint,
                             double *functional_local_adjoint,
                             double *weight_local_adjoint,
                             CoeffFct2D *Coefficients);

/******************************************************************************/
// prepares iteration for optimizing the parameter
/******************************************************************************/
void  PrepareIterationForOptimization(TCollection *coll,
                                      double *pw_const_param_old, 
                                      double *pw_const_param,
                                      double *pw_const_param_deriv_old,
                                      double *pw_const_param_deriv,
                                      double *search_direction,
                                      int *indicator_type_adjoint,
                                      double observed_functional,
                                      int loop_count,
                                      double &observed_functional_old,
                                      int &solve_adj,
                                      int &succ_step,
                                      double &step_length_adjoint,
                                      BoundCondFunct2D *BoundaryCond,
                                      CoeffFct2D *Coefficients);

int CheckFirstWolfeCondition(int N_Cells, double initial_functional, 
                             double final_functional, double alpha,
                             double *pw_const_param_deriv,
                             double *descent_direction);

int CheckSecondWolfeCondition(int N_Cells, double *pw_const_param_deriv, 
                              double *pw_const_param_deriv_old,
                              double *descent_direction);

int CheckGoldsteinConditions(int N_Cells, double initial_functional, 
                             double final_functional, double alpha,
                             double *pw_const_param_deriv_old,
                             double *search_direction);

void ComputeSearchDirection(double *pw_const_param,
                            double *pw_const_param_deriv_old,
                            double *pw_const_param_deriv,
                            double *search_direction,
                            double *search_direction_old,
                            double **bfgs_s, double **bfgs_y,
                            int N_Cells, int loop_count, int &bfgs_count, 
                            double step_length_adjoint);

int ComputeStepLengthParameter(double *pw_const_param_old, 
                               double *pw_const_param,
                               double *pw_const_param_deriv_old,
                               double *pw_const_param_deriv,
                               double *adjoint_funct_values,
                               double *search_direction_old,
                               double *search_direction,
                               double *oldsol,
                               double *sol,
                               double &observed_functional,
                               double observed_functional_initial,
                               int N_Cells, int N_Unknowns, int loop_count,
                               int bfgs_count,
                               int &inner_loop_count,
                               int N_adjoint_funct_values, 
                               double &step_length_adjoint,
                               double &step_length_adjoint_initial, 
                               double &observed_functional_old,
                               int &succ_step,
                               double *observed_functionals, 
                               double  *step_lengths);

void ComputeNextIterate(TCollection *coll,
                        TFESpace2D *pw_const_param_space,
                        CoeffFct2D *Coefficients,
                        int N_Cells,
                        double &step_length_adjoint,
                        double *pw_const_param_old,
                        double *pw_const_param,
                        double *pw_const_param_zero,
                        double *pw_const_param_deriv,
                        double *search_direction,
                        double *adjoint_funct_values,
                        double **bfgs_s,double **bfgs_y,
                        double observed_functional_old,
                        double initial_observed_functional,
                        double *step_lengths,
                        int &solve_adj,
                        int &loop_count, 
                        int &inner_loop_count, 
                        int N_adjoint_funct_values, 
                        int status_of_iteration);

int L_BFGS_SearchDirection(int N_, int loop_count, 
                           double *deriv_of_functional,
                           double *search_direction, 
                           double **bfgs_s, double **bfgs_y);

void L_BFGS(TCollection *coll, TFESpace2D *pw_const_param_space,
            CoeffFct2D *Coefficients, int N_, int &loop_count,
            double &step_length_adjoint, double *current_iterate,
            double *deriv_of_functional, double *search_direction,
            double *deriv_of_functional_old, 
            double **bfgs_s, double **bfgs_y);

void L_BFGS_Two_Loop_Recursion(int N_, int m, 
                               double gamma, double *deriv_of_functional,
                               double* &search_direction, double **bfgs_s, 
                               double **bfgs_y);
/******************************************************************************/
// DiscreteFormForAdjointProblem
// assigns the DiscreteForm for the adjoint problem
// the DiscreteForms are saved in the array DiscreteForms
/******************************************************************************/

void  DiscreteFormForAdjointProblem(TDiscreteForm2D *&DiscreteForm, 
                                    TDiscreteForm2D **DiscreteForms);



void ClearAdjointProblem2D(TFEFunction2D* &u_adjoint,
                           TSquareMatrix2D* &sqmatrixAadjoint,
                           TFESpace2D* &pw_const_param_space,
                           TFEVectFunct2D* &pw_const_param_fe,
                           TFEVectFunct2D* &pw_const_param_deriv_fe,
                           double* &sol_adjoint,
                           double* &rhs_edge,
                           double* &pw_const_param,
                           double* &pw_const_param_deriv,
                           double* &pw_const_param_deriv_old,
                           double* &pw_const_param_old,
                           double* &pw_const_param_zero,
                           double* &search_direction, 
                           double* &search_direction_old, 
                           double* &functional_adjoint, 
                           double* &weight_local_adjoint, 
                           double** &bfgs_s,
                           double** &bfgs_y,
                           int* &indicator_type_adjoint);

/******************************************************************************/
// DiscreteFormForAdjointProblem
// assigns the DiscreteForm for the adjoint problem
// the DiscreteForms are saved in the array DiscreteForms
/******************************************************************************/

void  DiscreteFormForAdjointProblem(TDiscreteForm2D *&DiscreteForm, 
                                    TDiscreteForm2D **DiscreteForms);

#endif
