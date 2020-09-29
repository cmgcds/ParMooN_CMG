// Generalized Stokes problem, exact solution is
// 
// u(x,y) = (sin(x)sin(y),cos(x)cos(y))
// p(x,y) = 2 cos(x)sin(y) - 2 sin(1)(1-cos(1))
//
// from: D.Braess and R.Sarazin, An efficient smoother for the Stokes
//          problem, Appl. Num. Math., 23: 3-19, 1997
//

void ExampleFile()
{
  OutPut("Example: BraessSarazinGeneralized.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================

void ExactU1(double x, double y, double *values)
{
  values[0] = sin(x)*sin(y);
  values[1] = cos(x)*sin(y);
  values[2] = sin(x)*cos(y);
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = cos(x)*cos(y);
  values[1] = -sin(x)*cos(y);
  values[2] = -cos(x)*sin(y);
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 2*(cos(x)*sin(y)-(double)sin(1.0)+(double)sin(1.0)*(double)cos(1.0));
  values[1] = -2*sin(x)*sin(y);
  values[2] = 2*cos(x)*cos(y);
  values[3] = 0;
}

void ExactNull(double x, double y, double *values)
{
  values[0] =0;
  values[1] =0;
  values[2] =0;
  values[3] =0;
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
  case 0: 
    value=0;
    break;
  case 1: 
    value=(double)sin(1.0)*sin(Param);
    break;
  case 2: 
    value=sin(1-Param)*(double)sin(1.0);
    break;
  case 3: 
    value=0;
    break;
  default: cout << "wrong boundary part number" << endl;
    break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
  case 0: 
    value=cos(Param);
    break;
  case 1: 
    value=(double)cos(1.0)*cos(Param);
    break;
  case 2: 
    value=cos(1.0-Param)*(double)cos(1.0);
    break;
  case 3: 
    value=cos(1.0-Param);
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
  double sigma = TDatabase::ParamDB->P8;
  int i;
  double *coeff,x,y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];

    coeff[0] = 1;
    coeff[1] = sigma *sin(x)*sin(y) ;               // f1
    coeff[2] = 4*cos(x)*cos(y) + sigma *cos(x)*cos(y) ; // f2
  }
}

