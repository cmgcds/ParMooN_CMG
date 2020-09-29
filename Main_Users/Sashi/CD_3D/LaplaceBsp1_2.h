// ======================================================================
// Laplace problem 3D
// f"ur PDE-Vorlesung, Bsp. 1.2
// ======================================================================
#include <ConvDiff3D.h>

void ExampleFile()
{
  OutPut("Example: LaplaceBsp1_2.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
    values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
    if (fabs(y)<1e-6)
	value = 1;
    else
	value = 0;
}

void BilinearCoeffs(int n_points, double *x, double *y, double *z,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = 1;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
    coeff[4] = 0;
    coeff[5] = 0;
  }
}

