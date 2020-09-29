// Navier-Stokes problem 
// backward facing step flow

#include <math.h>
#include <Constants.h>

// ========================================================================
// exact solution
// ========================================================================
// first component of the velocity u1
// values[0] = u1
// values[1] = \partial u1/ \partial x
// values[2] = \partial u1/ \partial y
// values[3] = Laplacian of u1
// if the functions are not known or if they are not needed, set them to 0

void ExampleFile()
{
  OutPut("Example: Channel30_diri.h" << endl) ;
}

void ExactU1(double x, double y, double *values)
{
  
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
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
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
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
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
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
  switch(BdComp)
  {
    case 0: 
      cond = DIRICHLET;
      break;
    case 1: 
      cond = NEUMANN;
      TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
      break;
    case 2: 
      cond = DIRICHLET;
      break;
    case 3: 
      cond = DIRICHLET;
      break;
  }
}

// boundary values of u1
// counting the boundary parts starts at 0
void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: 
      value=0;
      break;
    case 1: 
      value=0;
      break;
    case 2: 
      value=0;
      break;
    case 3: 
      // inflow boundary
      if (Param>=0.5)
        value=0;
      else
        value = 24*Param*(0.5-Param);
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
  double *coeff;

  // the coefficients are needed for a set of quadrature points
  // loop over all points
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    // coeff[0] is the viscosity
    coeff[0] = nu;
    
    // coeff[1] is the first component of the right hand side 
    coeff[1] = 0;
    coeff[2]=  0;           
  }
}

void InitialU1(double x, double y, double *values)
{
    return;
}
void InitialU2(double x, double y, double *values)
{
    return;
}
void InitialP(double x, double y, double *values)
{
    return;
}
