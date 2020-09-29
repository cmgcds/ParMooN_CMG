// Navier-Stokes problem 
// for example with analytical solution
// comparison of strong no-slip and weak no-slip
// strong no-slip case


#include <math.h>
#include <Constants.h>

void ExampleFile()
{
  OutPut("Example: atife_slip_sincos.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
// first component of the velocity u1
// values[0] = u1
// values[1] = \partial u1/ \partial x
// values[2] = \partial u1/ \partial y
// values[3] = Laplacian of u1
// if the functions are not known or if they are not needed, set them to 0

void ExactU1(double x, double y, double *values)
{

  values[0] = 2*Pi*sin(Pi*y)*cos(Pi*y)*sin(Pi*x)*sin(Pi*x)*sin(Pi*x);
  values[1] = 6*Pi*Pi*sin(Pi*y)*cos(Pi*y)*sin(Pi*x)*sin(Pi*x)*cos(Pi*x);
  values[2] = 2*Pi*Pi*cos(Pi*y)*cos(Pi*y)*sin(Pi*x)*sin(Pi*x)*sin(Pi*x)
   -2*Pi*Pi*sin(Pi*y)*sin(Pi*y)*sin(Pi*x)*sin(Pi*x)*sin(Pi*x);
  values[3] = 0;
}

// second component of the velocity u2
// values[0] = u2
// values[1] = \partial u2/ \partial x
// values[2] = \partial u2/ \partial y
// values[3] = Laplacian of u2
// if the functions are not known or if they are not needed, set them to 0

void ExactU2(double x, double y, double *values)
{
  values[0] = -3*Pi*sin(Pi*x)*sin(Pi*x)*cos(Pi*x)*sin(Pi*y)*sin(Pi*y);
  values[1] = -6* Pi*Pi*sin(Pi*x)*cos(Pi*x)*cos(Pi*x)*sin(Pi*y)*sin(Pi*y)
    +3*Pi*Pi*sin(Pi*x)*sin(Pi*x)*sin(Pi*x)*sin(Pi*y)*sin(Pi*y);
  values[2]= -6*Pi*Pi*sin(Pi*x)*sin(Pi*x)*cos(Pi*x)*sin(Pi*y)*cos(Pi*y);
  values[3] = 0;
}

// pressure p
// values[0] = p
// values[1] = \partial p/ \partial x
// values[2] = \partial p/ \partial y
// values[3] = Laplacian of p
// if the functions are not known or if they are not needed, set them to 0

void ExactP(double x, double y, double *values)
{
  values[0] = cos(Pi*x)+cos(Pi*y);
  values[1] = -Pi*sin(Pi*x);
  values[2] = -Pi*sin(Pi*y);
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================

// type of the boundary condition
// possible types:
// DIRICHLET, NEUMANN, SLIP_FRICTION_PENETRATION_RESISTANCE
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
  // no-slip on the whole boundary
  // slip on the whole boundary
  switch(BdComp)
  {
    case 0: cond =  SLIP_FRICTION_PENETRATION_RESISTANCE;
      break;
    case 1: cond =  SLIP_FRICTION_PENETRATION_RESISTANCE;
      break;
    case 2: cond =  SLIP_FRICTION_PENETRATION_RESISTANCE;
      break;
    case 3: cond =  SLIP_FRICTION_PENETRATION_RESISTANCE;
      break;
  }
  TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
}

// boundary values of u1
// counting the boundary parts starts at 0
void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

// boundary values of u2
void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

// ========================================================================
// coefficients for the Navier--Stokes equations
// viscosity \nu and right hand side f 
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  // set nu
  static double nu=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y;
  double u1,u2, u1x, u1y, u2x, u2y, u1xx, u1yy, u2xx, u2yy, px, py;
  double sinx, siny, cosx, cosy;

  // the coefficients are needed for a set of quadrature points
  // loop over all points
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    // coordinate of the quadrature point
    x = X[i];
    y = Y[i];

    // coeff[0] is the viscosity
    coeff[0] = nu;
    sinx = sin(Pi*x);
    siny = sin(Pi*y);
    cosx = cos(Pi*x);
    cosy = cos(Pi*y);
    px= -Pi*sinx; 
    py= -Pi*siny;
    u1  = 2*Pi*siny*cosy*sinx*sinx*sinx;
    u1x = 6*Pi*Pi*siny*cosy*sinx*sinx*cosx;
    u1y = 2*Pi*Pi*sinx*sinx*sinx*(cosy-siny)*(cosy+siny);
    u2 =  -3*Pi*sinx*sinx*cosx*siny*siny;
    u2x = 3*Pi*Pi*sinx*siny*siny*(-2*cosx*cosx+sinx*sinx);
    u2y = -6*Pi*Pi*siny*cosy*sinx*sinx*cosx;
    u1xx= -6*Pi*Pi*Pi*cosy*sinx*siny*(-2*cosx*cosx+sinx*sinx);
    u1yy= -8*Pi*Pi*Pi*siny*cosy*sinx*sinx*sinx;
    u2xx= 3*Pi*Pi*Pi*cosx*siny*siny*(-2*cosx*cosx+7*sinx*sinx);
    u2yy= -6*Pi*Pi*Pi*sinx*sinx*cosx*(cosy-siny)*(cosy+siny);

// coeff[1] is the first component of the right hand side 
    coeff[1] = -nu*(u1xx+u1yy)+ u1*u1x+u2*u1y+px;
    coeff[2]=  -nu*(u2xx+u2yy)+ u1*u2x+u2*u2y+py;           
  }
}

