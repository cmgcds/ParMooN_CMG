// Navier-Stokes problem 
// backward facing step flow

#include <math.h>
#include <Constants.h>

#define __BW_FACING_STEP__

void ExampleFile()
{
  OutPut("Example: Calotte_test.h, angle: "
         << TDatabase::ParamDB->P9 << endl);
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

    cond = DIRICHLET;
}

// boundary values of u1
// counting the boundary parts starts at 0
void U1BoundValue(int BdComp, double Param, double &value)
{
    value=0;
}

// boundary values of u2
void U2BoundValue(int BdComp, double Param, double &value)
{
    value=0;
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
  double *coeff, angle, angle1, r, x, y, phi;

  // current angle of the rotor
  // note: the rotor has two blades
  angle = -2*Pi*TDatabase::ParamDB->P9;
  // normalize to (-Pi,Pi]
  while (angle<=-Pi)
  {
      angle += 2*Pi;
  }
  if (angle<=0)
  {
      angle1 = angle + Pi;
  }
  else
  {
      angle1 = angle - Pi;
  }

  // the coefficients are needed for a set of quadrature points
  // loop over all points
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];
    
    // coeff[0] is the viscosity
    coeff[0] = nu;
    
    // coeff[1] is the first component of the right hand side 
    coeff[1] = 0;
    coeff[2]=  0;           
    // compute angle of the point 
    r = sqrt(x*x+y*y);
    if (r<=0.5)
    {
	// -pi < angle_1 <= pi 
	phi = atan2(y,x);
	// two blades, the jump of angle_1 at the negative x-axis has
	// to be taken into account
	// change right hand side
	// driving force b=2*Pi*P7*(y,-x,0)
	// the terms on the rhs are (b*\nabla)b
	if ((fabs(angle-phi)<0.175)|| (fabs(angle1-phi)<0.175)
	    || (fabs(angle-phi+2*Pi)<0.175)|| (fabs(angle1-phi+2*Pi)<0.175) 
	    || (fabs(angle-phi-2*Pi)<0.175)|| (fabs(angle1-phi-2*Pi)<0.175) )
	{
	    //OutPut(x << " " << y << " " << phi << " " << angle << " " << angle1 << endl);
	    //OutPut("rot " << x << " " << y << " " << z << endl);
	    coeff[1] = y*2*Pi;
	    coeff[2] = -x*2*Pi;
	}
    }
  }
}



// ========================================================================
// compute velocity profile
// ========================================================================
void ComputeVelocityProfile(TFEFunction2D *u1,int level)
{ 
  double x[3],z , y, h, value[3], bottom, ceiling;
  int step,xstp,i,j;

  x[0]=7;
  x[1]=15;

  bottom = -0.5;
  ceiling = 0.5;
   
  step = 200;
  h = (ceiling-bottom)/step;
  
  for (j=0;j<2;j++){
  for (i=0;i<=step;i++)
    {
      y = bottom + i*h;
      u1->FindGradient(x[j],y,value);
      OutPut("level " << level << " " << x[j] << " " << y
             << " " << value[0] << endl);
    }
    }

  bottom = 5.70;
  ceiling = 6.20;

  xstp=200;
  y=-0.499;
  h=(ceiling-bottom)/xstp;

  for (i=0;i<=xstp;i++){
    z = bottom + i*h;
    u1->FindGradient(z,y,value);
    OutPut("level " << level << " " << y << " " << z << " " << value[0] << endl);
  }
  
  xstp=800;
  y=-0.49;
  h=(ceiling-bottom)/xstp;
   
  for (i=0;i<=xstp;i++){
    z = bottom + i*h;
    u1->FindGradient(z,y,value);
    OutPut("level " << level << " " << y << " " << z << " " << value[0] << endl);
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

