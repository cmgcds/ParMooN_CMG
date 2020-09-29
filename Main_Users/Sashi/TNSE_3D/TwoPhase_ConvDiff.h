// ========================================================================
// exact solution
// ========================================================================
void ExactC(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

/// ========================================================================
/// initial solution
/// ========================================================================
void InitialC(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = 0;
}


/// ========================================================================
/// boundary conditions
/// ========================================================================
/// kind of boundary condition (for FE space needed)
void ConcenBoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = NEUMANN;
}

// value of boundary condition
void CBoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

void ConcenCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double t=TDatabase::TimeDB->CURRENTTIME;
  double *coeff, x, y, z;
  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    x = X[i];
    y = Y[i];
    z = Z[i];
    coeff[0] = eps;
  }  
} 
