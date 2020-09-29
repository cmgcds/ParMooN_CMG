// ======================================================================
// NSE2D_EquOrd_FixPo.h
//
// common declaration for all Navier-Stokes problems
// ======================================================================

#ifndef __NSE2D_EQUORD_FIXPO__
#define __NSE2D_EQUORD_FIXPO__

#include <Enumerations.h>

// ======================================================================
// Type 4, SDFEM, (grad u, grad v)
// ======================================================================
void NSType4SDFEMEquOrd(double Mult, double *coeff,
                        double *param, double hK,
                        double **OrigValues, int *N_BaseFuncts,
                        double ***LocMatrices, double **LocRhs);

void NSType4NLSDFEMEquOrd(double Mult, double *coeff,
                          double *param, double hK,
                          double **OrigValues, int *N_BaseFuncts,
                          double ***LocMatrices, double **LocRhs);

#endif // __NSE2D_EQUORD_FIXPO__
