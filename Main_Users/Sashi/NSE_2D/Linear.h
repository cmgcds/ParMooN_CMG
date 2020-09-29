// Navier-Stokes problem, exact solution is
// 
// u(x,y) = (y,x)
// p(x,y) = 0

void ExampleFile()
{
  OutPut("Example: Linear.h" << endl) ;
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = OSEEN_PROBLEM;
  OutPut("Example: Linear.h with INTERNAL_PROBLEM_IDENTITY " << 
	 TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY <<  endl) ;
}

// ========================================================================
// multi indices used for various things
// ========================================================================

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = y;
  values[1] = 0;
  values[2] = 1;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = x;
  values[1] = 1;
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
  cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=Param;
            break;
    case 2: value=1;
            break;
    case 3: value=1-Param;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=Param;
            break;
    case 1: value=1;
            break;
    case 2: value=1-Param;
            break;
    case 3: value=0;
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
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = x; // f1
    coeff[2] = y; // f2
    // for Oseen equations: solution is convection
    coeff[3] = y;
    coeff[4] = x;    
    if (TDatabase::ParamDB->STOKES_PROBLEM)
    {
	coeff[1] = 0; // f1
	coeff[2] = 0; // f2
	coeff[3] = 0; // b1
	coeff[4] = 0; // b2
    }
  }
}


