// Navier-Stokes problem, solution in ansatz space
// velocity pw quadratic, pressure linear
//
//
void ExampleFile()
{
  OutPut("Example: NoFlow.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 3*x-sin(y+4*z)-1.5-sin(1.0)*cos(1.0)
    -4*cos(1.0)*cos(1.0)*cos(1.0)*cos(1.0)*sin(1.0)
    +sin(1.0)*cos(1.0)*cos(1.0)-2*sin(1.0)*sin(1.0)*sin(1.0)
    +2*sin(1.0)*cos(1.0)*cos(1.0)*cos(1.0)+2*sin(1.0);
  values[1] = 3;
  values[2] = -cos(y+4*z);
  values[3] = -4*cos(y+4*z);
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y, z, u1, u2, u3, ux, uy, uz;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    x = X[i];
    y = Y[i];
    z = Z[i];
    coeff[0] = eps;
    coeff[1] = 3; // f1
    coeff[2] = - cos(y+4*z); // f2
    coeff[3] = - 4*cos(y+4*z);   // f3
  }
}
