// time dependent Navier-Stokes problem 3D, ansatz
//
//u1(t,x,y,z) = t^4*(y*(1-y)+z*(1-z))
//u2(t,x,y,z) = t^4*(x*(1-x)-z*(1-z))
//u3(t,x,y,z) = t^4*(2*x*(1-x)+y*(1-y))
//p(t,x,y,z) = 0

void ExampleFile()
{
  OutPut("Example: Bsp2.h" << endl);
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = t*t*t*t*(y*(1-y)+z*(1-z));
values[1] = 0;
values[2] = t*t*t*t*(1-2*y);
values[3] = t*t*t*t*(1-2*z);
values[4] = -4*t*t*t*t;
}

void InitialU2(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = t*t*t*t*(x*(1-x)-z*(1-z));
values[1] = t*t*t*t*(1-2*x);
values[2] = 0;
values[3] = t*t*t*t*(-1+2*z);
values[4] = 0;
}

void InitialU3(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = t*t*t*t*(2*x*(1-x)+y*(1-y));
values[1] = t*t*t*t*(2-4*x);
values[2] = t*t*t*t*(1-2*y);
values[3] = 0;
values[4] = -6*t*t*t*t;
}

void InitialP(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = 0;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = t*t*t*t*(y*(1-y)+z*(1-z));
values[1] = 0;
values[2] = t*t*t*t*(1-2*y);
values[3] = t*t*t*t*(1-2*z);
values[4] = -4*t*t*t*t;
}

void ExactU2(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = t*t*t*t*(x*(1-x)-z*(1-z));
values[1] = t*t*t*t*(1-2*x);
values[2] = 0;
values[3] = t*t*t*t*(-1+2*z);
values[4] = 0;
}

void ExactU3(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = t*t*t*t*(2*x*(1-x)+y*(1-y));
values[1] = t*t*t*t*(2-4*x);
values[2] = t*t*t*t*(1-2*y);
values[3] = 0;
values[4] = -6*t*t*t*t;
}

void ExactP(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = 0;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = t*t*t*t*(y*(1-y)+z*(1-z));
}

void U2BoundValue(double x, double y, double z, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = t*t*t*t*(x*(1-x)-z*(1-z));
}

void U3BoundValue(double x, double y, double z, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = t*t*t*t*(2*x*(1-x)+y*(1-y));
}

void U1BoundValue_diff(double x, double y, double z, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = 4*t*t*t*(y*(1-y)+z*(1-z));
}

void U2BoundValue_diff(double x, double y, double z, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = 4*t*t*t*(x*(1-x)-z*(1-z));
}

void U3BoundValue_diff(double x, double y, double z, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = 4*t*t*t*(2*x*(1-x)+y*(1-y));
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y, z;
  double t=TDatabase::TimeDB->CURRENTTIME;
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];
    z = Z[i];
    coeff[0] = eps;
    coeff[1] = 4*t*t*t*(y*(1-y)+z*(1-z))+4*eps*t*t*t*t+t*t*t*t*t*t*t*t*(x*(1-x)-z*(1-z))*(1-2*y)+t*t*t*t*t*t*t*t*(2*x*(1-x)+y*(1-y))*(1-2*z);
    coeff[2] = 4*t*t*t*(x*(1-x)-z*(1-z))+t*t*t*t*t*t*t*t*(y*(1-y)+z*(1-z))*(1-2*x)+t*t*t*t*t*t*t*t*(2*x*(1-x)+y*(1-y))*(-1+2*z);
    coeff[3] = 4*t*t*t*(2*x*(1-x)+y*(1-y))+6*eps*t*t*t*t+t*t*t*t*t*t*t*t*(y*(1-y)+z*(1-z))*(2-4*x)+t*t*t*t*t*t*t*t*(x*(1-x)-z*(1-z))*(1-2*y);
    coeff[4] = 12*t*t*(y*(1-y)+z*(1-z))+16*eps*t*t*t+8*t*t*t*t*t*t*t*(x*(1-x)-z*(1-z))*(1-2*y)+8*t*t*t*t*t*t*t*(2*x*(1-x)+y*(1-y))*(1-2*z);
    coeff[5] = 12*t*t*(x*(1-x)-z*(1-z))+8*t*t*t*t*t*t*t*(y*(1-y)+z*(1-z))*(1-2*x)+8*t*t*t*t*t*t*t*(2*x*(1-x)+y*(1-y))*(-1+2*z);
    coeff[6] = 12*t*t*(2*x*(1-x)+y*(1-y))+24*eps*t*t*t+8*t*t*t*t*t*t*t*(y*(1-y)+z*(1-z))*(2-4*x)+8*t*t*t*t*t*t*t*(x*(1-x)-z*(1-z))*(1-2*y);
  }
}
