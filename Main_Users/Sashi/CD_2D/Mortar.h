// ======================================================================
// Plane problem
// ======================================================================
#include <ConvDiff2D.h>

void ExampleFile()
{
  OutPut("Example: Mortar.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition_Mortar(int i, double t, BoundCond &cond)
{
   switch(i)
    {
	case 0:
	case 1:
	case 3:
	    cond = NEUMANN;
	 break;
	default:
	    cond = DIRICHLET;
    }
}

void BoundCondition_NonMortar(int i, double t, BoundCond &cond)
{
   switch(i)
    {
	case 1:
	case 2:
	case 3:
	    cond = NEUMANN;
	 break;
	default:
	    cond = DIRICHLET;
    }
}


// value of boundary condition
void BoundValue_Mortar(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
	case 0:
	case 1:
	case 3:
	    value = 0;
	 break;
	default:
	    value = 1.- Param;
  }
}

// value of boundary condition
void BoundValue_NonMortar(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
	case 1:
	case 2:
	case 3:
	    value = 0;
	 break;
	default:
	    value = 1.- Param;
  }
}

void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  static double eps=1./TDatabase::ParamDB->RE_NR;
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

    coeff[4] = 0;
  }
}

