
void ExampleFile()
{
  OutPut("Example: SimPaTurS_UnitSquare.h" << endl) ;
  TDatabase::ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD=1;
}

// ========================================================================
// exact solution
// ========================================================================

void ExactU1(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
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
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void InitialU1(double x, double y, double *values)
{
  values[0] = 0.125*y*(1.-y);
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0.0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition_NSE(int i, double t, BoundCond &cond)
{

  switch(i)
  {
   case 0: 
   case 2: 
   case 3: 
     cond = DIRICHLET;
   break;
   case 1: 
     cond = NEUMANN;
        TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=0;
   break;
   default: 
      cout << "wrong boundary part number" << endl;
      exit(0);
   break;
  }

}

void U1BoundValue(int BdComp, double Param, double &value)
{

  double velo= TDatabase::ParamDB->P1/0.25;

  switch(BdComp)
  {
  case 0: 
  case 1: 
  case 2: 
    value=0.;
  break;
  case 3: 
    value= velo*Param*(1.-Param);
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
  case 1: 
  case 2: 
  case 3:  
    value=0.;
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
  double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff,x,y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0;               // f1
    coeff[2] = 0; // f2
  }
}


// ========================================================================
// for scalar equations
// ========================================================================
void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  double a=0.2, b=0., c=0.;
  int i;
  double *coeff, *param;
  double x, y;

  eps = 0;
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = param[0];  // u1
    coeff[2] = param[1];  // u2
    coeff[3] = c;

    coeff[4] = 0.; // f
  }
}


void TInitial(double x, double y, double *values)
{
  values[0] = 0.5;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}


// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{

  switch(i)
  {
   case 0: 
   case 2: 
     cond = NEUMANN;
   break;
   case 3: 
     cond = DIRICHLET;
   break;
   case 1: 
     cond = NEUMANN;
   break;
   default: 
      cout << "wrong boundary part number" << endl;
      exit(0);
   break;
  }

}

void TBoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
  case 0: 
  case 2: 
    value=0.5;
    value=0.;
  break;
  case 1: 
    value=0;
  break;
  case 3: 
    value=1.;
   break;
  default: cout << "wrong boundary part number" << endl;
    break;
  }
}


