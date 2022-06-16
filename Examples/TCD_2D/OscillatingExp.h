// ======================================================================
// instationary problem
// ======================================================================

#define __SINCOS1__

/// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  OutPut("Example: Temperature.h" << endl); 
}
// exact solution
void Exact(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{

    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;
  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value = 0;
            break;
    case 2: 
            value = 0;

            break;
    case 3: value = 0;
            break;
  }
}

// initial conditon
void InitialCondition(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0;
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1./TDatabase::ParamDB->PE_NR;
  
  int i;
  double *coeff, *param;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;


  double b1=0., b2=0., c=0.;

  // previous discrete time
  tau = t-tau;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = b1;
    coeff[2] = b2;
    coeff[3] = c;
    coeff[4] = sin(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*t)*0.5;
  }
}
