#ifdef __2D__
void ApproximateRFBSolutionQuadNSE2D(TCollection *Coll, TFEFunction2D *u1, 
				     TFEFunction2D *u2, CoeffFct2D *coeff,
				     double *rhs);

void ApproximateRFBSolutionQuad_Q2_NSE2D(TCollection *Coll, TFEFunction2D *u1,
					 TFEFunction2D *u2, CoeffFct2D *Coeffs,
					 double *rhs);
#endif

#ifdef __3D__
void Compute_Q2_Value(double *coeff, double x, double y, double z, double *val);
void Compute_Q2_Value_Gradient(double *coeff, double x, double y, double z, double *val);

void ApproximateRFBSolutionQuadNSE3D(TCollection *Coll, TFEFunction3D *u1,
TFEFunction3D *u2, TFEFunction3D *u3, CoeffFct3D *Coeffs,
				     double *rhs);



void ApproximateTimeRFBSolutionQuadNSE3D(TCollection *Coll, TFEFunction3D *u1,
TFEFunction3D *u2, TFEFunction3D *u3, TFEFunction3D *p, CoeffFct3D *Coeffs,
double *rhs);

void ApproximateTimeRFBSolutionQuad_Q2_NSE3D(TCollection *Coll, TFEFunction3D *u1,
TFEFunction3D *u2, TFEFunction3D *u3, TFEFunction3D *p, CoeffFct3D *Coeffs, double *rhs);

void ApproximateTimeRFB_coupled_SolutionQuad_Q2_NSE3D(TCollection *Coll, TFEFunction3D *u1,
     TFEFunction3D *u2, TFEFunction3D *u3, TFEFunction3D *p, CoeffFct3D *Coeffs,
     double *rhs);

void ApproximateTimeRFB_coupled_cn_SolutionQuad_Q2_NSE3D(TCollection *Coll, TFEFunction3D *u1,
     TFEFunction3D *u2, TFEFunction3D *u3, TFEFunction3D *p, CoeffFct3D *Coeffs,
     double *old_small_scales, double *rhs);

void ApproximateTimeRFBSolutionQuad_cn_NSE3D(TCollection *Coll, TFEFunction3D *u1,
TFEFunction3D *u2, TFEFunction3D *u3, TFEFunction3D *p, CoeffFct3D *Coeffs,
					     double *old_small_scales, double *rhs);
#endif

