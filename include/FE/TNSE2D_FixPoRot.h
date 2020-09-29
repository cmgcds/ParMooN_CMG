// ======================================================================
// TNSE2D_FixPoRot.h     06/03/20
//
// common declaration for all time dependent Navier-Stokes problems
// rotation form of nonlinear term
// ======================================================================
// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType3SmagorinskyRot(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
			       double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType3SmagorinskyRotDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
				 double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType4SmagorinskyRot(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
			       double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType4SmagorinskyRotDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
				 double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLSmagorinskyRot(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
				   double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLSmagorinskyRotDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
				     double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Galerkin (grad u, grad v)
// ======================================================================
void TimeNSType3GalerkinRot(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Galerkin, D(u):D(v)
// ======================================================================
void TimeNSType3GalerkinRotDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Galerkin, (grad u, grad v)
// ======================================================================
void TimeNSType4GalerkinRot(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 4, Galerkin, D(u):D(v)
// ======================================================================
void TimeNSType4GalerkinRotDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Galerkin, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinRot(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 3, Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Galerkin, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinRotDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, Standard Galerkin + Div Term
// ======================================================================
void TimeNSType1GalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 2, Standard Galerkin + div term
// ======================================================================
void TimeNSType2GalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

// ======================================================================
// Type 1, Standard Galerkin + div term, only nonlinear part
// Type 2, Standard Galerkin + div term, only nonlinear part
// ======================================================================
void TimeNSType1_2NLGalerkinDiv(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

void TStokes_PSPG_GRADDIV(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

 // ======================================================================
// Type 14, PSPG for transient Stokes, nonlinear loop 
// ======================================================================
void TStokes_PSPG_GRADDIV_NL(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);

void TStokes_PSPG_GRADDIV_Rhs(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs);
