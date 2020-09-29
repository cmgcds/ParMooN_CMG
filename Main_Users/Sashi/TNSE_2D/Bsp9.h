// Navier-Stokes problem
// 
// u(x,y) = a(t)*(4*y*(1-y), 0)
// p(x,y) = 0
// a(t) = sin(Pi/2*t), wenn t<=1
//      = 1, sonst

void ExampleFile()
{
  OutPut("Example: Bsp9.h" << endl) ;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  double a;

  if(t<=1)
    a = sin(Pi/2*t);
  else
    a = 1;

  values[0] = a*4*y*(1-y);
  values[1] = 0;
  values[2] = a*(4-8*y);
  values[3] = -a*8;
}

void ExactU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
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
  double a;

  if(t<=1)
    a = sin(Pi/2*t);
  else
    a = 1;

  switch(BdComp)
  {
    case 0:
      value = 0;
    break;
    case 1:
      value = 4*a*Param*(1-Param);
    break;
    case 2:
      value = 0;
    break;
    case 3:
      value = 4*a*Param*(1-Param);
    break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

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

    if(t<=1)
    {
      coeff[1] = Pi/2*cos(Pi/2*t)*4*y*(1-y)+nu*8*sin(Pi/2*t);
      coeff[2] = 0;
    }
    else
    {
      coeff[1] = nu*8;
      coeff[2] = 0; 
    }
  }
}


