// Navier-Stokes problem with sine and cosine functions
// 

void ExampleFile()
{
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = OSEEN_PROBLEM;
  OutPut("Example: ExpLayer.h with INTERNAL_PROBLEM_IDENTITY " << 
	 TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY <<  endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
    double nu=1/TDatabase::ParamDB->RE_NR;
    double u, v, w;
    
    v = exp(1.0/nu);
    u = 1-v;
    u = 1/u;
    w = exp(y/nu);
    
    values[0] = u*(w-v);
    values[1] = 0;
    values[2] = u*w/nu;
    values[3] = u*w/(nu*nu);
    return;
}

void ExactU2(double x, double y, double *values)
{
    double nu=1/TDatabase::ParamDB->RE_NR;
    double u, v, w;
    
    v = exp(1.0/nu);
    u = 1-v;
    u = 1/u;
    w = exp(x/nu);
    
    values[0] = u*(w-v);
    values[1] = u*w/nu;
    values[2] = 0;
    values[3] = u*w/(nu*nu);
}

void ExactP(double x, double y, double *values)
{
    values[0] = x-y;
    values[1] = 1;
    values[2] = -1;
    values[3] = 0;
  return;
}

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
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
    double nu=1/TDatabase::ParamDB->RE_NR;
    double y;
  switch(BdComp)
  {
      case 0: y = 0;
            break;
      case 1: y = Param;
            break;
      case 2: y = 1;
            break;
      case 3: y = 1-Param;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  value = (exp(y/nu)-exp(1/nu))/(1-exp(1/nu));
  return;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
    double nu=1/TDatabase::ParamDB->RE_NR;
    double x;
  switch(BdComp)
  {
      case 0: x = Param;
            break;
      case 1: x = 1;
            break;
      case 2: x = 1-Param;
            break;
      case 3: x = 0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  value = (exp(x/nu)-exp(1/nu))/(1-exp(1/nu));
  return;
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
  double u1, u1x, u1y, u1lap, u2, u2x, u2y, u2lap, px, py;
  double u, v, w;

  v = exp(1.0/nu);
  u = 1-v;
  u = 1/u;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    x = X[i];
    y = Y[i];

    w = exp(y/nu);
    coeff[0] = nu;
    // prescribed solution
    u1 =   u*(w-v);
    u1x = 0;
    u1y =  u*w/nu;
    u1lap =  u*w/(nu*nu);

    w = exp(x/nu);    
    u2 =   u*(w-v);
    u2x =  u*w/nu;
    u2y = 0;
    u2lap = u*w/(nu*nu);
    px =  1;
    py = -1;

    coeff[1] = -nu*u1lap+u1*u1x+u2*u1y+px;
    coeff[2] = -nu*u2lap+u1*u2x+u2*u2y+py;

/*    coeff[3] = u1;
      coeff[4] = u2;*/
    coeff[3] = 1;
    coeff[4] = 1;
 
  }
}
