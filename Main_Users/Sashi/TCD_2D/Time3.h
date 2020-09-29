// ======================================================================
// instationary problem
// ======================================================================
#include <TimeConvDiff2D.h>

void ExampleFile()
{
  OutPut("Exzmple: Time3.h" << endl);
}

// exact solution
void Exact(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 1+2*x*t*t+3*t*y;
  values[1] = 2*t*t;
  values[2] = 3*t;
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double Param, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0:
      value = 1+2*Param*t*t;
    break;
    case 1:
      value = 1+2*t*t+3*Param*t;
    break;
    case 2:
      value = 1+3*t+2*(1-Param)*t*t;
    break;
    case 3:
      value = 1+3*t*(1-Param);
    break;
  } // endswitch
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  double a=0, b=0, c=0;
  int i;
  double *coeff, *param;
  double x, y;
  double t=TDatabase::TimeDB->CURRENTTIME;

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

    coeff[4] = c*(1+2*x*t*t+3*t*y)
              +2*a*t*t
              +3*t*b
              +3*y+4*x*t;
  }
}

// exact solution
void Initial(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 1+2*x*t*t+3*t*y;
  values[1] = 2;
  values[2] = 3*t;
  values[3] = 0;
}

