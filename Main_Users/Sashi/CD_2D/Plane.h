// ======================================================================
// Plane problem
// ======================================================================
// #include <ConvDiff2D.h>

void ExampleFile()
{
  OutPut("Example: Plane.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 1+2*x+3*y;
  values[1] = 2;
  values[2] = 3;
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=1+2*Param;
            break;
    case 1: value=3+3*Param;
            break;
    case 2: value=6-2*Param;
            break;
    case 3: value=4-3*Param;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps;
    coeff[1] = 1;
    coeff[2] = 2;
    coeff[3] = 0;

    coeff[4] = 1*2+2*3;
  }
}

