// ======================================================================
// problem in unit circle
// ======================================================================
#include <ConvDiff2D.h>

void ExampleFile()
{
  OutPut("Example: Circle.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 1-x*x-y*y;
  values[1] = -2*x;
  values[2] = -2*y;
  values[3] = -4;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  if(BdComp == 0)
    value = 0;
  else
    value = -3;
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  double a=0, b=0, c=1;
  int i;
  double *coeff, *param;
  double x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = a;
    coeff[2] = b;
    coeff[3] = c;

    coeff[4] = 4*eps-2*a*x-2*b*y+c*(1-x*x-y*y);
  }
}

