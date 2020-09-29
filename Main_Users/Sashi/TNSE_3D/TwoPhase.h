void ExampleFile()
{
  OutPut("Example: TwoPhase.h" << endl);
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = x;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = y;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = z;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  
  values[0] = x;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;

}

/// ========================================================================
/// initial solution
/// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0;
}

void InitialU2(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0;
}

void InitialU3(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0;

}

void InitialP(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0;
}

/// ========================================================================
/// boundary conditions
/// ========================================================================
/// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
  value = 0;    
}

/// ========================================================================
/// coefficients for Stokes form: A, B1, B2, f1, f2
/// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double t=TDatabase::TimeDB->CURRENTTIME;
  double *coeff, x, y, z, ratio_mu, ratio_rho;
  int phase = (int) coeffs[0][0];

//   cout << phase << endl;
  
  if ( phase == 1 )
  {
    ratio_mu = TDatabase::ParamDB->P5;
    ratio_rho = TDatabase::ParamDB->P6;
  }
  else 
  {
    ratio_mu = 1.;
    ratio_rho = 1.;
  }
  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    x = X[i];
    y = Y[i];
    z = Z[i];
    coeff[0] = ratio_mu*eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
    coeff[3] = ratio_rho*(-10); // f3
  }  
}