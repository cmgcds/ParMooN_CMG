// ======================================================================
// Example from P.W. Hemker
// ======================================================================
#include <ConvDiff2D.h>

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void BoundCondition(int BdComp, double t, BoundCond &cond)
{
  switch(BdComp)
  {
    case 1:
      cond = NEUMANN;
      break;
    default:
      cond = DIRICHLET;
  }
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 1:
      value = 0;
      break;
    case 4:
      value = 1;
      break;
    default:
      value = 0;
  }
}

void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  double angle = 0, v1, v2;
  int i;
  double *coeff, *param;

  v1 = cos(angle);
  v2 = sin(angle);

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = v1;
    coeff[2] = v2;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}

