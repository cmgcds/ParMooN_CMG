// ==========================================================================
// instationary problem
// ==========================================================================
#include <TimeConvDiff2D.h>
//===========================================================================
// example file
// =========================================================================

//==========================================================================
//a zero field with non-zero Direchlet boundary condition
//of the unit square (0,1)X(0,1)
//==========================================================================
#define __MOVINGMESH__ 

void ExampleFile()
{
  OutPut("Example: TimeDomianNoSource.h" << endl); 
  OutPut("__MOVINGMESH__" << endl); 
}
// exact solution
void Exact(double x, double y, double *values)
{    
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;

//   values[0] = exp(t)*(sin(2*Pi*x)*sin(2*Pi*y));
//   values[1] = exp(t)*2*Pi*cos(2*Pi*x)*sin(2*Pi*y);
//   values[2] = exp(t)*2*Pi*sin(2*Pi*x)*cos(2*Pi*y);
//   values[3] = 0;
  
// values[0] = 0;
// values[1] = 0;
// values[2] = 0;
// values[3] = 0;
  values[0] = 16.*(1.+0.5* sin(5.*Pi*t))*x*(1-x)*y*(1-y);
  values[1] = 16.*(1.+0.5* sin(5.*Pi*t))*(1-2.*x)*y*(1-y);
  values[2] = 16.*(1.+0.5* sin(5.*Pi*t))*x*(1-x)*(1-2.*y);
  values[3] = -32.*(1.+ 0.5 * sin(5.*Pi*t)) *(x - (x*x) + y - (y*y));  
  
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double t, BoundCond &cond)
{
 cond = DIRICHLET;
//   cond = NEUMANN;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{  
 value = 0;
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
//  values[0] = 16.*x*(1.-x)*y*(1.-y);
//   values[0] = 0.;
//  values[0] = 1. - sqrt(x*x + y*y);
  double t;

  t = TDatabase::TimeDB->CURRENTTIME;
//   values[0] = exp(t)*(sin(2*Pi*x)*sin(2*Pi*y)); 
  
  values[0] = 16.*x*(1. - x)*y*( 1. - y);  
  
}

void ModifyCoords(double x, double y, double &X, double &Y)
{
 double t = TDatabase::TimeDB->CURRENTTIME; 
 double xc, yc, r, tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH; 
 
  X = x*(2. - cos(10.*Pi*t));
  Y = y*(2. - cos(10.*Pi*t));
  
//     X = x;
//     Y = y;
  
//  
// //   X = x*cos(20.*t) - y*sin(20.*t);
// //   Y = x*sin(20.*t) + y*cos(20.*t);
 
//   xc = x-0.5;
//   yc = y-0.5;
//   r = sqrt(xc*xc + yc*yc);
//   if(r<0.5)
//   {
//    X +=  (0.5 - r)*(cos((Pi/2.)*t))*tau*xc ;
//    Y +=  (0.5 - r)*(cos((Pi/2.)*t))*tau*yc ;
//   }
//   
 
// cout << x << " ............ " << y << endl;
}


void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1./TDatabase::ParamDB->RE_NR;
  double a=0, b=0, c=0;
  int i;
  double *coeff, *param;
  double x, y;
  double t = TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps;
    coeff[1] = a;
    coeff[2] = b;
    coeff[3] = c;

//     coeff[4] = exp(t)*(sin(2*Pi*x)*sin(2*Pi*y))
// 	- eps * exp(t)*4*Pi*Pi*(-sin(2*Pi*x)*sin(2*Pi*y)-sin(2*Pi*x)*sin(2*Pi*y))
//        + a * exp(t)*2*Pi*cos(2*Pi*x)*sin(2*Pi*y)
//        + b * exp(t)*2*Pi*sin(2*Pi*x)*cos(2*Pi*y)
// 	+ c *  exp(t)*(sin(2*Pi*x)*sin(2*Pi*y));
    
       coeff[4] = 40.* Pi * cos(5.* Pi * t)*x*(1-x)*y*(1-y)  + eps * 32.*(1.+ 0.5 * sin(5.*Pi*t)) *(x - (x*x) + y - (y*y));    
    
  }
} // void BilinearCoeffs(in


// kind of boundary condition (for FE space needed)
void GridBoundCondition(int BdComp, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void GridBoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

 
