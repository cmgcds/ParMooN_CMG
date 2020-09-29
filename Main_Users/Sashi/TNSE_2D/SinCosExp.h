// time dependent Navier-Stokes problem 
//
// Tafti, Computers & Fluids (25), p. 647 - 665, 1996
// 
// u_1(x,y) = - cos(n pi x) sin(n pi y) exp(-2 (n pi)^2 nu t)
// u_2(x,y) = sin(n pi x) cos(n pi y) exp(-2 (n pi)^2 nu t)
// p(x,y) = -1/4 (cos(2n Pi x) + cos(2n Pi y))  exp(-4 (n pi)^2 nu t)
//
// n - number of oscillations

#define N_OSCILLATIONS    4

void ExampleFile()
{
  OutPut("Example: SinCosExp.h Number of oscillations: ");
  OutPut(N_OSCILLATIONS << endl) ;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
  double nu=1/TDatabase::ParamDB->RE_NR;
 
  fac = exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
  values[0] = -cos(N_OSCILLATIONS*pi*x)*sin(N_OSCILLATIONS*pi*y)*fac;
 
  //values[0] = -cos(N_OSCILLATIONS*pi*x)*sin(N_OSCILLATIONS*pi*y);
}

void InitialU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
  double nu=1/TDatabase::ParamDB->RE_NR;
 
  fac = exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
  values[0] = sin(N_OSCILLATIONS*pi*x)*cos(N_OSCILLATIONS*pi*y)*fac;
}

void InitialP(double x, double y, double *values)
{
  // values[0] = -(cos(2*N_OSCILLATIONS*pi*x)+cos(2*N_OSCILLATIONS*pi*y))/4.0;
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
  double nu=1/TDatabase::ParamDB->RE_NR;
 
  fac = exp(-4*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
  values[0] = -(cos(2*N_OSCILLATIONS*pi*x)+cos(2*N_OSCILLATIONS*pi*y))*fac/4.0;
}

f:=
// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
  double nu=1/TDatabase::ParamDB->RE_NR;
 
  fac = exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
  values[0] = -cos(N_OSCILLATIONS*pi*x)*sin(N_OSCILLATIONS*pi*y)*fac;
  values[1] = N_OSCILLATIONS*pi*sin(N_OSCILLATIONS*pi*x)*sin(N_OSCILLATIONS*pi*y)*fac;
  values[2] = -N_OSCILLATIONS*pi*cos(N_OSCILLATIONS*pi*x)*cos(N_OSCILLATIONS*pi*y)*fac;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
  double nu=1/TDatabase::ParamDB->RE_NR;
 
  fac = exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
  values[0] = sin(N_OSCILLATIONS*pi*x)*cos(N_OSCILLATIONS*pi*y)*fac;
  values[1] = N_OSCILLATIONS*pi*cos(N_OSCILLATIONS*pi*x)*cos(N_OSCILLATIONS*pi*y)*fac;
  values[2] = -N_OSCILLATIONS*pi* sin(N_OSCILLATIONS*pi*x)*sin(N_OSCILLATIONS*pi*y)*fac;;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
  double nu=1/TDatabase::ParamDB->RE_NR;
 
  fac = exp(-4*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
  values[0] = -(cos(2*N_OSCILLATIONS*pi*x)+cos(2*N_OSCILLATIONS*pi*y))*fac/4.0;
  values[1] = 2*N_OSCILLATIONS*pi*sin(2*N_OSCILLATIONS*pi*x)*fac/4.0;
  values[2] = 2*N_OSCILLATIONS*pi*sin(2*N_OSCILLATIONS*pi*y)*fac/4.0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  double nu=1/TDatabase::ParamDB->RE_NR;
  double fac;
  double pi = 3.14159265358979;

  switch(BdComp)
  {
  case 0: 
  case 2:
    value=0;
    break;
    case 1: 
      fac = exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
      value = -cos(N_OSCILLATIONS*pi)*sin(N_OSCILLATIONS*pi*Param)*fac;
      break;
    case 3: 
      fac = exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
      value = -sin(N_OSCILLATIONS*pi*(1-Param))*fac;
      break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U1BoundValue_diff(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  double nu=1/TDatabase::ParamDB->RE_NR;
  double fac;
  double pi = 3.14159265358979;

  switch(BdComp)
  {
  case 0: 
  case 2:
    value=0;
    break;
    case 1: 
      fac = -2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
      value = -cos(N_OSCILLATIONS*pi)*sin(N_OSCILLATIONS*pi*Param)*fac;
      break;
    case 3: 
      fac = -2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
      value = -sin(N_OSCILLATIONS*pi*(1-Param))*fac;
      break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  double nu=1/TDatabase::ParamDB->RE_NR;
  double fac;
  double pi = 3.14159265358979;
  value = 0;
  switch(BdComp)
  {
  case 1: 
  case 3:
    value=0;
    break;
  case 0: 
    fac = exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
    value = sin(N_OSCILLATIONS*pi*Param)*fac;
    break;
  case 2: 
    fac = exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
    value = sin(N_OSCILLATIONS*pi*(1-Param))*cos(N_OSCILLATIONS*pi)*fac;
    break;
  default: cout << "wrong boundary part number" << endl;
    break;
  }
}
void U2BoundValue_diff(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  double nu=1/TDatabase::ParamDB->RE_NR;
  double fac;
  double pi = 3.14159265358979;
  value = 0;
  switch(BdComp)
  {
  case 1: 
  case 3:
    value=0;
    break;
  case 0: 
    fac = -2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
    value = sin(N_OSCILLATIONS*pi*Param)*fac;
    break;
  case 2: 
    fac = -2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*nu*t);
    value = sin(N_OSCILLATIONS*pi*(1-Param))*cos(N_OSCILLATIONS*pi)*fac;
    break;
  default: cout << "wrong boundary part number" << endl;
    break;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff, x, y; 
  double nu=1/TDatabase::ParamDB->RE_NR;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = nu;

    coeff[1] = 0;
    coeff[2] = 0;
  }
}


