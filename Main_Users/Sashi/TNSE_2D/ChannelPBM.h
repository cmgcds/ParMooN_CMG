// Navier-Stokes problem
// 
// u(x,y) = unknown
// p(x,y) = unknown

void ExampleFile()
{
  OutPut("Example: ChanelPBM.h" << endl) ;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=0;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
//   if(y>=1./3. && y<=2./3.)
//    values[0] = 6.*(2.-3.*y)*(3.*y-1.);
//   else 
   values[0] = 0.0;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0.0;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0.;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0.;
  values[1] = 0.;
  values[2] = 0.;
  values[3] = 0.;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
   switch(i)
   {
     case 0:   cond = DIRICHLET;
            break;
    case 1:   cond = NEUMANN;
            break;
    case 2:   cond = DIRICHLET;
            break;
    case 3:   cond = DIRICHLET;
            break;
    default: cout << "wrong boundary part number: " << i << endl;
   }
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0: value = 0.;
            break;
    case 1: value = 0.;
            break;
    case 2: value = 0.;
            break;
    case 3:
	if((Param>=1./3.)&&(Param<=2./3.))
	 value = 6.*(2.-3.*Param)*(3.*Param-1.);
	else 
	 value = 0.0;
	    break;
        default: cout << "wrong boundary part number: " << BdComp << endl;
  }

}

void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value = 0.0;
            break;
    case 1: value = 0.0;
            break;
    case 2: value = 0.0;
            break;
    case 3: value = 0.0;
	   break;
        default: cout << "wrong boundary part number: " << BdComp << endl;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  static double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0.; // f1
    coeff[2] = 0.; // f2
  }
}

