// ======================================================================
// convection skew to the domain
// Hughes, Mallet, Mizukami 1986
// ======================================================================

#include <ConvDiff2D.h>

void ExampleFile()
{
  OutPut("Example: SinSin.h" << endl) ;
}
// exact solution (this is the solution for eps = 0)
void Exact(double x, double y, double *values)
{
    values[0] = sin(Pi*x)*sin(Pi*y);
    values[1] = Pi*cos(Pi*x)*sin(Pi*y);
    values[2] = Pi*sin(Pi*x)*cos(Pi*y);
    values[3]= -Pi*Pi*sin(Pi*x)*sin(Pi*y)-Pi*Pi*sin(Pi*x)*sin(Pi*y);
}

// kind of boundary condition
void BoundCondition(int i, double t, BoundCond &cond)
{
    if (i==1)
      cond = NEUMANN;
    else 
      cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;

    switch(BdComp)
    {
	case 0:
	case 3:
	    value = 0;
	    break;
	case 1:
	    value = eps*Pi*cos(Pi)*sin(Pi*Param);
	    break;
	case 2:
	    value = sin(Pi*(1-Param))*sin(Pi);
	    break;
    }
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;
  double x, y, arg;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];
    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = 2;
    coeff[2] = 3;
    coeff[3] = 1;
    arg =  coeff[3] * sin(Pi*x)*sin(Pi*y);
    arg += coeff[1] * Pi*cos(Pi*x)*sin(Pi*y);
    arg += coeff[2] * Pi*sin(Pi*x)*cos(Pi*y);
    arg -= coeff[0] * (-2 * Pi*Pi* sin(Pi*x)*sin(Pi*y));
    coeff[4] = arg;
  }
}

void CheckWrongNeumannNodes(TCollection *Coll, TFESpace2D *fespace,
			    int &N_neum_to_diri, int* &neum_to_diri,
			    int* &neum_to_diri_bdry, 
			    double* &neum_to_diri_param)
{
    exit(1);
}
