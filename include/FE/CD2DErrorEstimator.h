// =======================================================================
// @(#)CD2DErrorEstimator.h        1.3 11/15/99
//
// Class:       TCD2DErrorEstimator
//
// Purpose:     define error estimator
//
// Author:      Volker John
//
// History:     29.01.98 start implementation
//
// =======================================================================

#ifndef __CD2D_ERROR_ESTIMATOR__
#define __CD2D_ERROR_ESTIMATOR__

#ifdef __2D__

#include <FEFunction2D.h>

#define cd_gradient_indicator                      0
#define cd_residual_estimator_h1                   1
#define cd_residual_estimator_l2                   2
#define cd_residual_estimator_energy_quasi_robust  3
#define cd_gradient_recovery                       4
#define cd_implicit_estimator_neumann              5

class TCD2DErrorEstimator
{
  protected:
    /** cell collection which the spaces are based on */
    TCollection *Collection;

    /** ansatz space */
    TFESpace2D *FESpace2D;

    /** FE function */
    TFEFunction2D *FEFunction2D;

    /** error estimator routine */
    int FELocalEstimator;

    /** do or not do error control */
    int ErrorControl;

public:
    /** initialize error estimator */
    TCD2DErrorEstimator(int fe_local_estimator,
                       TFEFunction2D *fe_function, int error_control);

    /** return Collection */
    TCollection *GetCollection()
      { return Collection; };

    /** return ansatz space */
    TFESpace *GetAnsatzSpace()
      { return FESpace2D; };

    /** return */
    TFEFunction2D *GetFEFunction2D()
      { return FEFunction2D;};

    /** return element error estimator */
    int GetFELocalEstimator()
      { return FELocalEstimator;};

    /** return error control */
    int GetErrorControl()
      { return ErrorControl;}


    /** the estimator */
    void GetErrorEstimate(int N_Derivatives,
                          MultiIndex2D *NeededDerivatives,
                          CoeffFct2D *Coeff,
                          BoundCondFunct2D **BoundaryConds,
                          BoundValueFunct2D **BoundaryValues,
                          TAuxParam2D *Aux,
                          int n_fespaces,
                          TFESpace2D **fespaces,
                          double *eta_K,
                          double *eta_max,
                          double *estimated_global_error);

    /** elementwise estimator */
    void  EstimateCellError_new(TFESpace2D *fespace,
                            TBaseCell *cell,
                            int N_Points,
                            double *X,
                            double *Y,
                            double *AbsDetjk,
                            double *weights,
                            double **Derivatives,
                            double **AuxArray,
                            BoundCondFunct2D **BoundaryConds,
                            BoundValueFunct2D **BoundaryValues,
                            int N_Points1D,
                            double *zeta,
                            double *X1D[4],
                            double *Y1D[4],
                            double *weights1D,
                            double *xyval_ref1D[4],
                            double *xderiv_ref1D[4],
                            double *yderiv_ref1D[4],
			    int MaxN_BaseFunctions2D_loc,
                            int *GlobalNumbers,
                            int *BeginIndex,
                            int *DOF,
                            double *Values,
                            double *local_error);

    void  EstimateCellError(TFESpace2D *fespace,
                            TBaseCell *cell,
                            int N_Points,
                            double *X,
                            double *Y,
                            double *AbsDetjk,
                            double *weights,
                            double **Derivatives,
                            double **AuxArray,
                            BoundCondFunct2D **BoundaryConds,
                            BoundValueFunct2D **BoundaryValues,
                            int N_Points1D,
                            double *zeta,
                            double *X1D[4],
                            double *Y1D[4],
                            double *weights1D,
                            double *xyval_ref1D[4],
                            double *xderiv_ref1D[4],
                            double *yderiv_ref1D[4],
			    int *GlobalNumbers,
                            int *BeginIndex,
                            int *DOF,
                            double *Values,
                            double *local_error);

};
#endif

#endif // #ifdef __2D__