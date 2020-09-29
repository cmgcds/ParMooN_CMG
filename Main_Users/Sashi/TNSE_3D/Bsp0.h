// time dependent Navier-Stokes problem 3D, ansatz
//
//u1(t,x,y,z) = 2*t
//u2(t,x,y,z) = 3
//u3(t,x,y,z) = 4
//p(t,x,y,z) = 0

void ExampleFile()
{
  OutPut("Example: Bsp0.h" << endl);
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = 2*t;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
}

void InitialU2(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = 3;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
}

void InitialU3(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = 4;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
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

values[0] = 2*t;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
}

void ExactU2(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = 3;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
}

void ExactU3(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

values[0] = 4;
values[1] = 0;
values[2] = 0;
values[3] = 0;
values[4] = 0;
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
value = 2*t;
}

void U2BoundValue(double x, double y, double z, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = 3;
}

void U3BoundValue(double x, double y, double z, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = 4;
}

void U1BoundValue_diff(double x, double y, double z, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = 2;
}

void U2BoundValue_diff(double x, double y, double z, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = 0;
}

void U3BoundValue_diff(double x, double y, double z, double &value)
{
double t=TDatabase::TimeDB->CURRENTTIME;
value = 0;
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
    coeff[1] = 2;
    coeff[2] = 0;
    coeff[3] = 0;
    coeff[4] = 0;
    coeff[5] = 0;
    coeff[6] = 0;
  }
}
