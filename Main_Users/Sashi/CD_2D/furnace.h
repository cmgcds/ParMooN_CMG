// ======================================================================
// Sine problem
// ======================================================================


void ExampleFile()
{
  OutPut("Example: furnace.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
  if(BdComp==3)
    cond = DIRICHLET;
  else
    cond = NEUMANN;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  static double eps=1/TDatabase::ParamDB->PE_NR;

  if(BdComp==3)
   {
    value = 90;
   }
  else if (BdComp==2)
   {
    value = 750*(1./TDatabase::ParamDB->PE_NR);
   }
   else
   { value = 0; }
}

void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  static double eps=1./TDatabase::ParamDB->PE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    //double *param = parameters[i];

    coeff[0] = eps;
    coeff[1] = 1.;
    coeff[2] = 0;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}

