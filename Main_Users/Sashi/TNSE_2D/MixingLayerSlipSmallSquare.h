// Navier-Stokes problem, Driven Cavity
// 
// u(x,y) = ?
// p(x,y) = ?

#define U_INFTY 1   

void ExampleFile()
{
  OutPut("Example: MixingLayerSlipSmallSquare.h" << endl) ;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double z;

  z = (2*y-1)*28;
  if (z>=0)
    values[0] = U_INFTY * (1-exp(-2*z))/(1+exp(-2*z));
  else
    values[0] = U_INFTY * (exp(2*z)-1)/(exp(2*z)+1);    
  z = z/2;
  values[0] -= 0.1* U_INFTY *exp(-z*z)*2*(y-0.5)*cos(8*Pi*x)*28*28;
  values[0] -= 0.1* U_INFTY *exp(-z*z)*2*(y-0.5)*cos(20*Pi*x)*28*28;
}

void InitialU2(double x, double y, double *values)
{
  double z;

  z = (y-0.5)*28;
  values[0] = 0.1*U_INFTY*exp(-z*z)*sin(8*Pi*x)*8*Pi;
  values[0] += 0.1*U_INFTY*exp(-z*z)*sin(20*Pi*x)*20*Pi;
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
  cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
  TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff, x, y; 
  static double nu=1/TDatabase::ParamDB->RE_NR;

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


