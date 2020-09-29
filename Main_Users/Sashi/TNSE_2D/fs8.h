// Navier-Stokes problem
// 
// u(x,y) = (t^2 y, tx)
// p(x,y) = 0

void ExampleFile()
{
  OutPut("Example: Bsp6.h" << endl) ;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*t*y*y;
}

void InitialU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 2*t*t*x*x;
}

void InitialP(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*x+y-0.5*(1+t);
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*t*y*y;
  values[1] = 0;
  values[2] = t*t*2*y;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*t*2*x*x;
  values[1] = t*t*4*x;
  values[2] = 0;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = t*x+y-0.5*(1+t);
  values[1] = t;
  values[2] = 1;
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
  double t=TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0:
      value = 0;
    break;
    case 1:
      value = t*t*Param*Param;
    break;
    case 2:
      value = t*t*1;
    break;
    case 3:
      value = t*t*(1-Param)*(1-Param);
    break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0:
      value = t*t*2*Param*Param;
    break;
    case 1:
      value = t*t*2;
    break;
    case 2:
      value = t*t*2*(1-Param)*(1-Param);
    break;
    case 3:
      value = 0;
    break;
  }
}

void U1BoundValue_diff(int BdComp, double Param, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0:
      value = 0;
    break;
    case 1:
      value = 2*t*Param;
    break;
    case 2:
      value = 2*t;
    break;
    case 3:
      value = 2*t*(1-Param);
    break;
  }
  value = 0;
}

void U2BoundValue_diff(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0:
      value = Param;
    break;
    case 1:
      value = 1;
    break;
    case 2:
      value = (1-Param);
    break;
    case 3:
      value = 0;
    break;
  }
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  static double nu = 1/TDatabase::ParamDB->RE_NR;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double t1, t2, t3, t5, t6, t8, t9, t10, t14, t16;
  double t17, t22, t23, t35, t38, t40;
  int i;
  double *coeff, x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = nu;

    coeff[1] = 2*t*y*y -2*nu*t*t + 4*x*x*y*t*t*t*t + t;
    coeff[2] = 4*t*x*x -4*nu*t*t + 4*y*y*x*t*t*t*t + 1;
    coeff[3] = 1;
    coeff[4] = 0;
    
  }
}


