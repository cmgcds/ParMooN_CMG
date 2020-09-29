// ==========================================================================
// instationary problem
// ==========================================================================

//===========================================================================
// example file
// =========================================================================

void ExampleFile()
{
  OutPut("Example: TEST_RFB.h" << endl);
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
    values[4] = 0;
}

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{
    values[0] = 0;
}

// kind of boundary condition
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
    cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
    value = 0.0;
}

void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
        double **parameters, double **coeffs)
{
  double eps = 1.0/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;                                  // *param;
  double x, y, z, c, a[3], b[3], s[3], h;
  // double t = TDatabase::TimeDB->CURRENTTIME;
  
  // norm of convection vector
  c = 1.0/16 + 1.0/64 + 1;
  c = sqrt(c);

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    // param = parameters[i];

    x = X[i];
    y = Y[i];
    z = Z[i];

    // diffusion
    coeff[0] = 1;
    // convection in x direction
    coeff[1] = 0;
    // convection in y direction
    coeff[2] = 0;
    // convection in z direction
    coeff[3] = 0;
    // reaction
    coeff[4] = 0;
    // rhs
    coeff[5] = 1;
    // rhs from previous time step
    coeff[6] = 0;
  }
}

