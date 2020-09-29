#include <math.h>

#define PI 3.14159265

void ExactC_Surf(double x, double y, double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  
  values[0] = exp(-6*t)*x*y+0.5; 
}

/// ========================================================================
/// initial solution
/// ========================================================================
void InitialC_Surf(double x, double y, double z, double *values)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
//   values[0] = 0.5;
//   values[0] = 0.5 + 0.5*z;
  values[0] = exp(-6*t)*x*y+0.5; 
}

/// ========================================================================
/// boundary conditions on surface
/// ========================================================================
/// kind of boundary condition (for FE space needed)
void SurfaceBoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}


/// ========================================================================
/// coefficients 
/// ========================================================================
void SurfCoeffs(int n_points, double *X, double *Y, double *Z,
                double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->PE_NR;
  int i;
  double t=TDatabase::TimeDB->CURRENTTIME;
  double *coeff, x, y, z;
  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    x = X[i];
    y = Y[i];
    z = Z[i];
    coeff[0] = eps;
//     coeff[1] = 0; // f1
//     coeff[2] = 0; // f2
//     coeff[3] = ratio_rho*g; // f3
  }  
}



//==================================================================
// domain velocity ( for testing surface equations )
//==================================================================
void ExactG1(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  double r = sqrt(x*x+y*y+z*z);
  double rho = sqrt(x*x+y*y);
  double omega = 2*PI*cos(4*PI*t);
  double fact = 0.1*cos(2*PI*z)/**cos(2*PI*y)*cos(2*PI*x)*/;

//   values[0] = x/(r*r*r);
  
  values[0] = omega*y;

//   values[0] = 0;

//   fact = (1.+0.25*sin(t));
//   fact *= fact;
//   values[0] = 0.25*x*x*cos(t) / fact;
  
//   values[0] = x / r * fact * sin(t);
}

void ExactG2(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  double r = sqrt(x*x+y*y+z*z);
  double rho = sqrt(x*x+y*y);
  double omega = 2*PI*cos(4*PI*t);
  double fact = 0.1*cos(2*PI*z)/**cos(2*PI*y)*cos(2*PI*x)*/;

//   values[0] = y/(r*r*r);;
  
  values[0] = -omega*x;

//   values[0] = 0;

//   values[0] = y / r * fact * sin(t);
}

void ExactG3(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  double r = sqrt(x*x+y*y+z*z);
  double fact = 0.1*cos(2*PI*z)/**cos(2*PI*y)*cos(2*PI*x)*/;

//   values[0] = z/(r*r*r);  
  
  values[0] = 0;
//   values[0] = z / r * fact * sin(t);
}

