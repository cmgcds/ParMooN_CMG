/****************************************************************************************
 *                                                                                      *
 *                         FDM_test.h                                                   *
 *                        -------------                                                 *
 *                                                                                      *
 ***************************************************************************************/

#ifndef __FDM_TEST__
#define __FDM_TEST__

void FWE_FDM_Upwind_2D(TCollection *coll,
		       TFEFunction2D *fesol,
		       CoeffFct2D *Coeffs, BoundValueFunct2D *bound_val,
		       int *dof_conversion, int first);

#endif


