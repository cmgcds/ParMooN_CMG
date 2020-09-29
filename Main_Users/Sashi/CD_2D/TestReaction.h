// ======================================================================
// test example for reaction-dominated problem
// ======================================================================
#include <ConvDiff2D.h>
#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>
#include <FE2D.h>
#include <FEDesc2D.h>

void ExampleFile()
{
  OutPut("Example: TestReaction.h" << endl) ;
}
// exact solution (this is the solution for eps = 0)
void Exact(double x, double y, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3]= 0;
}

// kind of boundary condition
void BoundCondition(int i, double t, BoundCond &cond)
{
    if (i==0)
	cond = DIRICHLET;
    else
	cond = NEUMANN;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
    if ((BdComp==0) && (Param>0.47) && (Param<0.53))
	value = 1;
    else
	value = 0;
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
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
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 1;
    coeff[4] = 0;
  }
}
