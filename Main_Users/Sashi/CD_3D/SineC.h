// ======================================================================
// Sine problem 3D
// ======================================================================
#include <ConvDiff3D.h>

void ExampleFile()
{
  OutPut("Example: SineC.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  values[0] = sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[1] = Pi*cos(Pi*x)*sin(Pi*y)*sin(Pi*z);
  values[2] = Pi*sin(Pi*x)*cos(Pi*y)*sin(Pi*z);
  values[3] = Pi*sin(Pi*x)*sin(Pi*y)*cos(Pi*z);
  values[4] = -3*Pi*Pi*sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
  value = sin(Pi*x)*sin(Pi*y)*sin(Pi*z);
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
    coeff[5] = (3*Pi*Pi*eps+1)*sin(Pi*x)*sin(Pi*y)*sin(Pi*z)
               +4*Pi*cos(Pi*x)*sin(Pi*y)*sin(Pi*z)
               +3*Pi*sin(Pi*x)*cos(Pi*y)*sin(Pi*z)
               +2*Pi*sin(Pi*x)*sin(Pi*y)*cos(Pi*z);
  }
}

