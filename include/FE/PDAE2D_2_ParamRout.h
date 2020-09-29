// ======================================================================
// @(#)TNSE2D_ParamRout.h        1.2 05/05/00
//
// common declaration for time dependent semiconductor device equation
// ======================================================================

#ifndef __PDAE2D_2_PARAMROUT__
#define __PDAE2D_2_PARAMROUT__

// ======================================================================
// setting for error calculation for all types
// ======================================================================
MultiIndex2D TimeNSAllDerivatives[3] = { D00, D10, D01 };

// ========================================================================
// parameter routines
// ========================================================================

// ========================================================================
// parameters: u1old, u2old
// ========================================================================
void PDAE2D_2_Params2(double *in, double *out);

int PDAE_N_FESpaces2 = 1;
int PDAE_N_Fct2 = 2;
int PDAE_N_ParamFct2 = 1; // richtig
int PDAE_N_FEValues2 = 6;
int PDAE_N_Params2 = 6;
int PDAE_FEFctIndex2[6] = { 0, 0, 0, 1, 1, 1 };
MultiIndex2D PDAE_FEMultiIndex2[6] = { D10, D01, D00, D10, D01, D00 };
ParamFct *PDAE_Fct2[1] = { PDAE2D_2_Params2 };
int PDAE_BeginParam2[1] = { 0 };

#endif
