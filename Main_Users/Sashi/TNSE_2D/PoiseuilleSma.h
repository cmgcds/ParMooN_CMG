// Navier-Stokes problem, Poiseuille-Problem
// 
// u(x,y) = (4*y*(1-y), 0)
// p(x,y) = nu*(4-8*x)

void ExampleFile()
{
  OutPut("Example: PoiseuilleSma.h" << endl) ;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  values[0] = 4*y*(1-y);
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y, double *values)
{
  double nu=1/TDatabase::ParamDB->RE_NR;

  values[0] = nu*(4-8*x);
}


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
  static double nu=1/TDatabase::ParamDB->RE_NR;
  values[0] = nu*(4-8*x);
  values[1] = -8*nu;
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
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;
  double c,tmp,delta,d12,d21,normD,dnormD_dy,dd12_dy,x,y;

  tmp = 1;
  delta =  CharacteristicFilterWidth(tmp);

  c = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT*delta*delta;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;

    x = X[i];
    y = Y[i];

    d12 = 2 - 4*y;
    d21 = d12;
    normD = sqrt(d12*d12+d21*d21);
    if (y<0.5)      
      dnormD_dy = -4*sqrt(2.0);
    else
    {
      if (y>0.5)
        dnormD_dy = 4*sqrt(2.0);
      else
        dnormD_dy = 0;
    }
    dd12_dy = -4;

    coeff[1] = c*(-normD*dd12_dy-dnormD_dy*d12); // f1
    coeff[2] = 0; // f2
  }
}

