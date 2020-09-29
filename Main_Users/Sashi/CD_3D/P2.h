// ======================================================================
// Sine problem 3D
// ======================================================================
#include <ConvDiff3D.h>

void ExampleFile()
{
  OutPut("Example: P2.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  values[0] = x*x+y*z;
  values[1] = 2*x;
  values[2] = z;
  values[3] = y;
  values[4] = 2;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
  value = x*x+y*z;
}

void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;
  double x, y, z;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];
    x = X[i];
    y = Y[i];
    z = Z[i];

    coeff[0] = eps;
    coeff[1] = 4;
    coeff[2] = 3;
    coeff[3] = 2;
    coeff[4] = 1;
    coeff[5] = -eps * 2 + 4 * 2*x + 3 * z + 2 * y + 1 * (x*x+y*z);
  }
}

