// Stokes problem, exact solution is
// 
// u(x,y) = (sin(2 pi x)cos(2 pi y), - cos(2 pi x)sin(2 pi y))
// p(x,y) = x^2 + y^2 + C
//
// from: Cai, Tang, Vassilevski, Wang: A mixed finite element 
//       method for Stokes Equations Based on 
//       Pseudostress-Velocity formulation

void ExampleFile()
{
  OutPut("Example: CTVW05.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================

void ExactU1(double x, double y, double *values)
{
  values[0] = sin(2*Pi*x)*cos(2*Pi*y);
  values[1] = 2*Pi*cos(2*Pi*x)*cos(2*Pi*y);
  values[2] = -2*Pi*sin(2*Pi*x)*sin(2*Pi*y);
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = -cos(2*Pi*x)*sin(2*Pi*y);
  values[1] = 2*Pi*sin(2*Pi*x)*sin(2*Pi*y);;
  values[2] = -2*Pi*cos(2*Pi*x)*cos(2*Pi*y);;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = x*x+y*y-2.0/3.0;
  values[1] = 2*x;
  values[2] = 2*y;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
    cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double t, double &value)
{
    double x,y;
    
    switch(BdComp)
    {
	case 0:
	    x = t; y = 0;
	    break;
	case 1:
	    x = 1; y = t;
	    break;
	case 2:
	    x = 1 - t; y = 1;
	    break;
	case 3:
	    x = 0; y = 1 - t;
	    break;
    }
    value = sin(2*Pi*x)*cos(2*Pi*y);
}
void U2BoundValue(int BdComp, double t, double &value)
{

    double x,y;
    
    switch(BdComp)
    {
	case 0:
	    x = t; y = 0;
	    break;
	case 1:
	    x = 1; y = t;
	    break;
	case 2:
	    x = 1 - t; y = 1;
	    break;
	case 3:
	    x = 0; y = 1 - t;
	    break;
    }
    value = -cos(2*Pi*x)*sin(2*Pi*y);
}




// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff,x,y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = 2*(4*Pi*Pi*sin(2*Pi*x)*cos(2*Pi*y)+x);               // f1
    coeff[2] = 2*(-4*Pi*Pi*cos(2*Pi*x)*sin(2*Pi*y)+y); // f2
  }
}

