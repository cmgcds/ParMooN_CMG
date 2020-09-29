// ======================================================================
// TCD3D.h
//
// common declaration for time dependent convection diffusion problems
// ====================================================================== 

#include <Constants.h>
#include <Enumerations.h>

namespace TCD3D_Surf
{
  static int N_Terms = 4;
  static MultiIndex3D Derivatives[4] = { D100, D010, D001, D000 };
  static int FESpaceNumbers[4] = { 0, 0, 0, 0 };
  static int N_Matrices = 2;
  static int N_Rhs = 0;
  static int RowSpaces[2] = { 0, 0 };
  static int ColumnSpaces[2] = { 0, 0 };
  static int *RhsSpaces = NULL;
  
  void MatricesAssemble (double Mult, double *coeff, double *param, double hK,
			 double **OrigValues, int *N_BaseFuncts,
			 double ***LocMatrices, double **LocRhs);
			 
  // parameters
  static int N_Parameters = 4;
  static int N_FEFct = 3;
  static int N_FEVectFct = 1;
  static int FctIndex[4] = { 0, 1, 2, 3};
  static MultiIndex3D ParamDerivatives[4] = { D100, D010, D001, D000 };
  static int IsGlobal[4] = { 1, 1, 1, 1 };
}