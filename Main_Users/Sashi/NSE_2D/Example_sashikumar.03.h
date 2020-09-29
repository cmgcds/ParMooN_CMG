// Navier-Stokes problem 
// for example with analytical solution
// comparison of strong no-slip and weak no-slip
// strong no-slip case


#include <math.h>
#include <Constants.h>
void ExampleFile()
{
  OutPut("Example: Example_sashikumar.01.h" << endl) ;
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
  double t1, t2, t3, t4, t5, t6, t7, t9, t10, t13, t15;

  t1 = x*x;
  t2 = 1.0-x;
  t3 = t2*t2;
  t4 = t1*t3;
  t5 = 1.0-y;
  t6 = t5*t5;
  t7 = y*t6;
  t9 = y*y;
  t10 = t9*t5;
  t13 = x*t3;
  t15 = t1*t2;
  values[0] = 2.0*t4*(t7-t10);
  values[1] = 4.0*(t13-t15)*(t7-t10);
  values[2] = 2.0*t4*(t6-4.0*y*t5+t9);
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
  double t1, t2, t3, t4, t5, t6, t7, t9, t10, t20, t23;

  t1 = 1.0-x;
  t2 = t1*t1;
  t3 = x*t2;
  t4 = y*y;
  t5 = 1.0-y;
  t6 = t5*t5;
  t7 = t4*t6;
  t9 = x*x;
  t10 = t9*t1;
  t20 = y*t6;
  t23 = t4*t5;
  values[0] = -2.0*t7*(t3-t10);
  values[1] = -2.0*t7*(t2-4.0*x*t1+t9);
  values[2] = 4.0*(t23-t20)*(t3-t10);
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
  values[0] = x*x*x+y*y*y-0.5;
  values[1] = 3.0*x*x;
  values[2] = 3.0*y*y;
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
    case 0: cond = DIRICHLET;
      break;
    case 1: cond = DIRICHLET;
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
  double sinx, siny, cosx, cosy, delta=TDatabase::ParamDB->P5;
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t15;


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
    
    t1 = x*x;
    t2 = 1.0-x;
    t3 = t2*t2;
    t4 = t1*t3;
    t5 = 1.0-y;
    t6 = t5*t5;
    t7 = y*t6;
    t8 = x*t2;
    t9 = y*y;
    t10 = t9*t5;
    t11 = t9*t6;
    t12 = y*t5;
    t13 = x*t3;
    t15 = t1*t2;

    u1  = 2.0*t4*(t7-t10);
    u1x = 4.0*(t13-t15)*(t7-t10);
    u1y = 2.0*t4*(t6-4.0*t12+t9);
    u2 = -2.0*t11*(t13-t15);
    u2x = -2.0*t11*(t3-4.0*t8+t1);
    u2y = 4.0*(t10-t7)*(t13-t15);
    u1xx= 4.0*(t7-t10)*(t3-4.0*t8+t1);
    u1yy= 12.0*t4*(2.0*y-1);
    u2xx= -12.0*t11*(2.0*x-1);
    u2yy=-4.0*(t6-4.0*t12+t9)*(t13-t15);
    px = 3.0*t1;
    py = 3.0*t9;	

// coeff[1] is the first component of the right hand side
    if (TDatabase::ParamDB->STOKES_PROBLEM )
    {
      coeff[1] = -nu*(u1xx+u1yy)+ px;
      coeff[2]=  -nu*(u2xx+u2yy)+ py;  
    }
    else
    {       
      coeff[1] = -nu*(u1xx+u1yy)+ u1*u1x+u2*u1y+px;
      coeff[2]=  -nu*(u2xx+u2yy)+ u1*u2x+u2*u2y+py;  
    }
    // coeff[3] = x*(1-x)*y*(1-y)-delta*delta*(-2*y*(1-y)-2*x*(1-x))/24;
  }
}
