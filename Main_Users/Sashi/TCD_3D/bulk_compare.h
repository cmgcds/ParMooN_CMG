// ======================================================================
// Sine problem 3D
// ======================================================================

void ExampleFile()
{
  OutPut("Example: bulk_compare.h" << endl) ;
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

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{
    values[0] = 0;
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
    if ((fabs(x-0.5)<1e-6)&&(fabs(y-0.5)<1e-6)&&(fabs(z)<1e-6))
	value = 1e5;
}

void BilinearCoeffs(int n_points, double *x, double *y, double *z,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR, X, Y;
  int i;
  double *coeff, *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];
    X = x[i];
    Y = y[i];


    coeff[0] = eps;
    coeff[1] = Y;
    coeff[2] = X;
    coeff[3] = 2e-8;
    coeff[4] = 0;
    coeff[5] = 0;
  }
}

