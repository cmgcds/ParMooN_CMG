// ======================================================================
// Laplace-Gleichung zur Illustration der Mehrgitter-Vorlesung
// f = 1, homogene Dirichlet R.B.
// ======================================================================
#include <ConvDiff2D.h>
#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>
#include <FE2D.h>
#include <FEDesc2D.h>

void ExampleFile()
{
  OutPut("Example: MG_Vorlesung_rhs1.h" << endl) ;
}
// exact solution (this is the solution for eps = 0)
void Exact(double x, double y, double *values)
{

    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}

// kind of boundary condition
void BoundCondition(int i, double t, BoundCond &cond)
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
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = 1.0;
    coeff[1] = 0.0;
    coeff[2] = 0.0;
    coeff[3] = 0;
    coeff[4] = 1.0;
  }
}

