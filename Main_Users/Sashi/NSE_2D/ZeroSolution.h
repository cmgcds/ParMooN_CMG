// Stokes problem, exact solution is
// 
// u(x,y) = (0,0)
// p(x,y) = 0
//

void ExampleFile()
{
  OutPut("Example: ZeroSolution.h" << endl) ;
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
  cond = DIRICHLET;
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

    coeff[0] = 1;
    coeff[1] = 0;               // f1
    coeff[2] = 0;               // f2
  }
}

