// ======================================================================
// instationary problem
// ======================================================================
#include <TimeConvDiff2D.h>

void ExampleFile()
{
  OutPut("Example: Time1.h" << endl);
}

// exact solution
void Exact(double x, double y, double *values)
{
  double r2, t;

  t = TDatabase::TimeDB->CURRENTTIME;
  r2 = (x-0.5-0.2*t)*(x-0.5-0.2*t) + (y-0.5-0.1*t)*(y-0.5-0.1*t);

  if( r2 < 0.0625 )
    values[0] = 0.25 - sqrt(r2);
  else
    values[0] = 0;

  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  double a=0.2, b=0.1, c=0;
  int i;
  double *coeff, *param;
  double x, y;

  eps = 0;
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

    coeff[4] = 0;
  }
}

// exact solution
void Initial(double x, double y, double *values)
{
  double r2;

  r2 = (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5);

  if( r2 < 0.0625 )
    values[0] = 0.25 - sqrt(r2);
  else
    values[0] = 0;
}

