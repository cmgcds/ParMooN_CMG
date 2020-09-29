// ======================================================================
// @(#)ConvDiff.h        12/06/26
//
// common declaration for all convection diffusion problems
// ======================================================================

#ifndef __CONVDIFF__
#define __CONVDIFF__

/******************************************************************************/
// SetParametersCD
// sets parameters of the data base for the main programs
void SetParametersCD(int &nonlinear_method);

double Mesh_size_in_convection_direction(double hK, double b1, double b2);

double Mesh_size_in_convection_direction_without_storing(double hK, double b1, 
                                                         double b2);


double Mesh_size_in_convection_direction(double hK, double b1, double b2,
                                         double b3);
double Mesh_size_in_convection_direction_without_storing(double hK, double b1,
                                                         double b2, double b3);



double Compute_SDFEM_delta(double hK, double eps, double b1, double b2,
#ifdef __3D__
                           double b3,
#endif
                           double react, double linfb);



double Compute_SOLD_sigma(double hK, double eps, double b1, double b2, 
#ifdef __3D__
                          double b3, 
#endif
                          double c, double f, double linfb, double deltaK,
                          double *param, double residual, int residual_computed,
                          int time_dependent_problem);

/** coercivity constant of problem 
 * used eg for residual based estimator of Verf"uhrt 2005 */
double EstimateCoercivityConstant(TCollection *Coll, 
#ifdef __2D__
                                  CoeffFct2D *Coeff
#else // 3D
                                  CoeffFct3D *Coeff
#endif
                                 );

/** it should be i>=0 and i<= 31, otherwise error*/
void SetSoldParameters(int i);


#ifdef __2D__
void EdgeStabilization(TFESpace2D *fespace,  TFEFunction2D *u, 
                       CoeffFct2D *Coeffs, double *rhs, int time_dependent,
                       double *time_step, TFEFunction2D *old_u);
#endif

double ComputeAlpha(double hK);

#endif // __CONVDIFF__
