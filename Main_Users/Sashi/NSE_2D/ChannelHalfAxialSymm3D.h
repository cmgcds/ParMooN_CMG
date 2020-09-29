// Navier-Stokes problem 
// for example with analytical solution
// comparison of strong no-slip and weak no-slip
// strong no-slip case


#include <math.h>
#include <Constants.h>
void ExampleFile()
{
  OutPut("Example: ChannelHalfPolar.h" << endl) ;
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
  values[0] =1-y*y/0.25;
  values[1] =0;
  values[2] =-8*y;
  values[3] =0;
}

// second component of the velocity u2
// values[0] = u2
// values[1] = \partial u2/ \partial x
// values[2] = \partial u2/ \partial y
// values[3] = Laplacian of u2
// if the functions are not known or if they are not needed, set them to 0

void ExactU2(double x, double y, double *values)
{
  values[0] =0;
  values[1] =0;
  values[2] =0;
  values[3] =0;
}

// pressure p
// values[0] = p
// values[1] = \partial p/ \partial x
// values[2] = \partial p/ \partial y
// values[3] = Laplacian of p
// if the functions are not known or if they are not needed, set them to 0

void ExactP(double x, double y, double *values)
{
  static double nu=1/TDatabase::ParamDB->RE_NR;

  values[0] =-16*nu*x;
  values[1] =-16*nu;
  values[2] =0;
  values[3] =0;
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
  switch(BdComp)
  {
     case 0:  cond =  SLIP_FRICTION_PENETRATION_RESISTANCE;
       TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
      break;
     case 1: cond = NEUMANN;
        TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;      
        break;
     case 2: cond = DIRICHLET;
        break;
     case 3: cond = DIRICHLET;
        break;
  }
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
    case 3: value=1-(1-Param)*(1-Param);
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
  double sinx, siny, cosx, cosy, delta=TDatabase::ParamDB->P5;


  // the coefficients are needed for a set of quadrature points
  // loop over all points
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    // coeff[0] is the viscosity
    coeff[0] = nu;
    coeff[1] = 0;
    coeff[2]= 0;
   
  }
}

