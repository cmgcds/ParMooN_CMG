// Navier-Stokes problem, Poiseuille-Problem
// 
// u(x,y) = (4*y*(1-y), 0)
// p(x,y) = x-1/2

// ========================================================================
// multi indices used for various things
// ========================================================================
MultiIndex2D StokesATestDerivatives[] = { D10, D01 };
MultiIndex2D StokesAAnsatzDerivatives[] = { D10, D01 };

MultiIndex2D NavierStokesATestDerivatives[] = { D10, D01, D00, D00 };
MultiIndex2D NavierStokesAAnsatzDerivatives[] = { D10, D01, D10, D01 };

MultiIndex2D L2H1ErrorDerivatives[] = { D00, D10, D01 };

MultiIndex2D NavierStokesB1TestDerivatives[] = { D00 };
MultiIndex2D NavierStokesB1AnsatzDerivatives[] = { D10 };
MultiIndex2D NavierStokesB2TestDerivatives[] = { D00 };
MultiIndex2D NavierStokesB2AnsatzDerivatives[] = { D01 };

MultiIndex2D NavierStokesRhsDerivatives[] = { D00 };

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 4*y*(1-y);
  values[1] = 0;
  values[2] = 4-8*y;
  values[3] = -8;
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
  values[0] = x-0.5;
  values[1] = 1;
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
    case 1: value=4*Param*(1-Param);
            break;
    case 2: value=0;
            break;
    case 3: value=4*Param*(1-Param);
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
// routine for evaluating errors
// ========================================================================
void L2H1Errors(int N_Points, double *X, double *Y, double *AbsDetjk, 
                double *Weights, double **Der, double **Exact,
                double *LocError)
{
  int i;
  double *deriv, *exactval, w, t;

  LocError[0] = 0.0;
  LocError[1] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    t = deriv[0]-exactval[0];
    LocError[0] += w*t*t;

    t = deriv[1]-exactval[1];
    LocError[1] += w*t*t;
    t = deriv[2]-exactval[2];
    LocError[1] += w*t*t;
  } // endfor i

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

// ========================================================================
// bilinear forms for matrix blocks
// ========================================================================
void StokesACoeffs(int n_points, double *x, double *y,
                   double **parameters, double **coeffs)
{
  static double eps = 1;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = eps;
  }
}

void NavierStokesACoeffs(int n_points, double *x, double *y,
                   double **parameters, double **coeffs)
{
  static double eps = 1;
  int i;
  double *coeff, *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps;
    coeff[1] = eps;
    coeff[2] = param[0];
    coeff[3] = param[1];
    // coeff[2] = y[i];
    // coeff[3] = x[i];
    
/*
    cout << "param[0]: " << param[0] << endl;
    cout << "y[i]: " << y[i] << endl;
    cout << "param[1]: " << param[1] << endl;
    cout << "x[i] " << x[i] << endl;
    cout << endl;
*/
  }
}

void NavierStokesB1Coeffs(int n_points, double *x, double *y,
            double **parameters, double **coeffs)
{
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = -1;
  }
}

void NavierStokesB2Coeffs(int n_points, double *x, double *y,
            double **parameters, double **coeffs)
{
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = -1;
  }
}

// ========================================================================
// linear forms for right-hand sides
// ========================================================================
void NavierStokesRhs1Coeffs(int n_points, double *X, double *Y,
                double **parameters, double **coeffs)
{
  static double eps=1;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    coeff[0] = 1+8*eps;
  }
}

void NavierStokesRhs2Coeffs(int n_points, double *X, double *Y,
                double **parameters, double **coeffs)
{
  static double eps = 1;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    coeff[0] = 0;
  }
}

// ========================================================================
// parameter routine
// ========================================================================
void NavierStokesParams(double *in, double *out)
{
  out[0] = in[2];
  out[1] = in[3];
}

