// ======================================================================
// parabolic layers
// ======================================================================
#define __PARABOLIC_LAYERS__

#include <ConvDiff3D.h>
#include <Joint.h>
#include <BoundFace.h>
#include <BoundComp.h>
#include <FE3D.h>
#include <FEDesc3D.h>

void ExampleFile()
{
  OutPut("Example: ParabolicLayers3D.h" << endl) ;
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
    value = 0;
}

void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 1;
    coeff[2] = 0;
    coeff[3] = 0;
    coeff[4] = 0;
    coeff[5] = 1;
  }
}

