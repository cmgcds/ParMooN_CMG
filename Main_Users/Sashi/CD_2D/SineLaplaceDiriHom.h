// ======================================================================
// Sine problem
// ======================================================================
#include <ConvDiff2D.h>

void ExampleFile()
{
  OutPut("Example: SineLaplaceDiriHom.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = sin(Pi*x)*sin(Pi*y);
  values[1] = Pi*cos(Pi*x)*sin(Pi*y);
  values[2] = Pi*sin(Pi*x)*cos(Pi*y);
  values[3] = -2*Pi*Pi*sin(Pi*x)*sin(Pi*y);
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;

    value = 0;
}

void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;

    coeff[4] = 2*Pi*Pi*eps*sin(Pi*x[i])*sin(Pi*y[i]);
  }
}

