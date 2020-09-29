// Stokes problem, exact solution is
// 
// u(x,y) = (sin(x)sin(y),cos(x)cos(y))
// p(x,y) = 2 cos(x)sin(y) - 2 sin(1)(1-cos(1))
//
// from: D.Braess and R.Sarazin, An efficient smoother for the Stokes
//          problem, Appl. Num. Math., 23: 3-19, 1997
//

void ExampleFile()
{
  OutPut("Example: PressureGradient.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================

void ExactU1(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = y*y*y;
  values[1] = 0;
  values[2] = 3*y*y;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  switch(i)
  {
  case 0: 
  case 2: 
    cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
    TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
    //   cond = DIRICHLET;
    break;
  case 1: 
  case 3: 
    cond = DIRICHLET;
    break;
  }
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  value=0;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value=0;
}


// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff,x,y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = 0;               // f1
    coeff[2] = 3*y*y; // f2
  }
}

