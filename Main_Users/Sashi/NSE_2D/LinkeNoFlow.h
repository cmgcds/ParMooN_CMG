// Navier-Stokes problem from Friedhelms Habil
// 

void ExampleFile()
{
  OutPut("Example: LinkeNoFlow.h, ") ;
  TDatabase::ParamDB->STOKES_PROBLEM = 1;
  OutPut(" set TDatabase::ParamDB->STOKES_PROBLEM to 1, "); 
  OutPut(" amplification factor (P9)" << TDatabase::ParamDB->P9 << endl);   
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
  values[0] = TDatabase::ParamDB->P9*(x*x*x+y*y*y+x-1);
  values[1] = TDatabase::ParamDB->P9*(3*x*x+1);
  values[2] = TDatabase::ParamDB->P9*3*y*y;
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

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  double nu=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    x = X[i];
    y = Y[i];

    coeff[0] = nu;
    coeff[1] = TDatabase::ParamDB->P9*(3*x*x+1);
    coeff[2] = TDatabase::ParamDB->P9*3*y*y;
  }
}

