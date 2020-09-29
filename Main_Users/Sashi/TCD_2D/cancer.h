// ======================================================================
// instationary problem
// ======================================================================

/// ========================================================================
// example file
// ========================================================================

#define __CANCER__

void ExampleFile()
{
  OutPut("Example: cancer.h" << endl); 
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] =  0;
}

// =========================================================================
// kind of boundary condition (for FE space needed)
// =========================================================================
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
 cond = NEUMANN;
}

void V_BoundCondition(int BdComp, double t, BoundCond &cond)
{
 cond = NEUMANN;
}

void W_BoundCondition(int BdComp, double t, BoundCond &cond)
{
 cond = NEUMANN;
}

// =========================================================================
// value of boundary condition
// =========================================================================
void BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

void V_BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

void W_BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

// =========================================================================
// initial conditon
// =========================================================================
void InitialCondition(double x, double y, double *values)
{
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = exp(t)*(sin(2*Pi*x)*sin(2*Pi*y));
  
}

void V_InitialCondition(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = exp(t)*(sin(2*Pi*x)*sin(2*Pi*y));
  
}
void W_InitialCondition(double x, double y, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  values[0] = exp(t)*(sin(2*Pi*x)*sin(2*Pi*y));
  
}
// =========================================================================
// linear coeffs
// =========================================================================
void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->PE_NR;
  double b1=0, b2=0, c=0;
  int i;
  double *coeff, *param;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

//     x = X[i];
//     y = Y[i];

    coeff[0] = eps;
    coeff[1] = b1;
    coeff[2] = b2;
    coeff[3] = c;
    coeff[4] = 0.;
  }
}

void V_BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->PE_NR;
  double b1=0, b2=0, c=0;
  int i;
  double *coeff, *param;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

//     x = X[i];
//     y = Y[i];

    coeff[0] = eps;
    coeff[1] = b1;
    coeff[2] = b2;
    coeff[3] = c;
    coeff[4] = 0.;
  }
}

void W_BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->PE_NR;
  double b1=0, b2=0, c=0;
  int i;
  double *coeff, *param;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

//     x = X[i];
//     y = Y[i];

    coeff[0] = eps;
    coeff[1] = b1;
    coeff[2] = b2;
    coeff[3] = c;
    coeff[4] = 0.;
  }
}


// =========================================================================

