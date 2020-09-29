// Navier-Stokes problem, Ferro Magnet
// 
// u(x,y) = unknown
// p(x,y) = unknown

void ExampleFile()
{
  OutPut("Example: FerroMagnet.h" << endl) ;
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
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  switch(i)
  {
  case 2:
  case 3:
  case 4:
  case 5:
  case 8:
    cond = DIRICHLET;
    break;
  case 0:
  case 1:
  case 6:
  case 7:
    cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
    TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
    TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
    cond = DIRICHLET;
    break;
  default:
    OutPut("Wrong Boundary Number !!!" << endl);
    exit(1);
  }    
}

void U1BoundValue(int BdComp, double Param, double &value)
{
     value = 0;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 1; // f1
    coeff[2] = 0; // f2
  }
}

