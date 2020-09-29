
#include <stdlib.h>
#include <math.h>
#include <Database.h>


#include <Convolution.h>
#include <MooNMD_Io.h>
#include <ConvDiff.h>

#include <TNSE3D_Routines.h>
#include <MainUtilities.h>

void SolveGridEquation(double **Entries, double *sol, double *rhs,
                       int *KCol, int *RowPtr, int N_DOF);
// ======================================================================
// declaration for InitializeDiscreteFormGrid()
// ======================================================================

static int GridN_Terms_3D = 3;
static MultiIndex3D GridDerivatives_3D[3] = { D100, D010, D001 };
static int GridSpaceNumbers_3D[3] = { 0, 0 , 0 };
static int GridN_Matrices_3D = 9;
static int GridRowSpace_3D[9] = { 0, 0 , 0,  0, 0 , 0 ,  0, 0 , 0  };
static int GridColumnSpace_3D[9] = { 0, 0 , 0,  0, 0 , 0 ,  0, 0 , 0  };
static int GridN_Rhs_3D = 0;
static int *GridRhsSpace_3D = NULL;

void GridAssemble4(double Mult, double *coeff,
                  double *param, double hK,
                  double **OrigValues, int *N_BaseFuncts,
                  double ***LocMatrices, double **LocRhs);