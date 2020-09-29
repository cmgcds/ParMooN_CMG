
void ExampleFile()
{
  OutPut("Example: CirculateNSE.h" << endl) ;
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

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  switch(i)
  {

   case 0:
     cond = NEUMANN;
   break;

   case 1:
   case 2:
   case 3:
     cond = DIRICHLET;
   break;

   default: 
      cout << "wrong boundary part number" << endl;
      exit(0);
   break;
  }

}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double val;
  double U_In = 1.0;

  switch(BdComp)
  {
  case 0:
     value=0.;
  break;

  case 1:
  case 2:
    value=0.;
  break;
  case 3:
      if ((Param>1.0/3.0)&& (Param<2.0/3.0))
       {
        val = 3.*(Param-1./3.);
        value = U_In*val*( 1 - val );
       }
      else
       {
        value = 0;
       }

  break;

  default: cout << "wrong boundary part number" << endl;
    break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
    value=0;
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

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = 0.;  // f1
    coeff[2] = 0.; // f2


  }
}

