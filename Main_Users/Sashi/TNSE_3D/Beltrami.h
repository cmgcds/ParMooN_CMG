// Navier-Stokes problem
//
//
// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  OutPut("Example: Beltrami.h; a "<< TDatabase::ParamDB->P5 << " d " <<  
         TDatabase::ParamDB->P6 << endl);
}

// ========================================================================
// exact solution
// ========================================================================
void InitialU1(double x, double y, double z, double *values)
{
  double a = TDatabase::ParamDB->P5;
  double d = TDatabase::ParamDB->P6;
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = -a*(exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y))*exp(-eps*d*d*t);
}

void InitialU2(double x, double y, double z, double *values)
{
  double a = TDatabase::ParamDB->P5;
  double d = TDatabase::ParamDB->P6;
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = -a*(exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z))*exp(-eps*d*d*t);
}

void InitialU3(double x, double y, double z, double *values)
{
  double a = TDatabase::ParamDB->P5;
  double d = TDatabase::ParamDB->P6;
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = -a*(exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x))*exp(-eps*d*d*t);
}

void InitialP(double x, double y,  double z, double *values)
{
  double a = TDatabase::ParamDB->P5;
  double d = TDatabase::ParamDB->P6;
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double t=TDatabase::TimeDB->CURRENTTIME;
  
  double trans = -0.75*exp(-2*eps*d*d*t)*a*(-d*d*exp(2*a)*sin(a+d)-a*a*cos(a-d)+exp(a)*cos(2*a+d)*d*d+d*d*exp(2*a)*cos(a+d)+d*d*exp(2*a)*sin(2*a+2*d)-d*d*exp(a)*sin(d)+a*a*exp(a)*cos(a+2*d)+a*a*exp(a)*cos(2*a+d)-cos(a-d)*d*d-cos(a+d)*d*d+a*a*exp(2*a)*cos(a-d)+3*d*d*exp(a)*sin(a)-d*d*exp(a)*sin(2*a+d)-a*a*exp(2*a)*sin(2*a+2*d)-exp(a)*sin(a+2*d)*d*d+a*a*exp(2*a)*sin(a-d)+a*a*exp(a)*sin(2*a+d)+exp(a)*cos(a+2*d)*d*d+a*a*exp(a)*sin(a+2*d)-3*a*a*exp(a)*sin(a)+a*a*exp(2*a)*sin(a+d)-sin(a-d)*d*d-a*a*sin(a+d)+sin(a+d)*d*d+a*a*sin(a-d)+a*a*exp(a)*sin(d)-d*d*exp(2*a)*sin(a-d)+d*d*d*d-exp(2*a)*d*d*cos(2*a+2*d)-2*a*d*exp(2*a)*cos(a+d)-exp(2*a)*d*d*d*d+exp(2*a)*a*a*a*a-2*a*cos(a-d)*d-a*a*exp(2*a)*cos(2*a+2*d)-a*a*cos(a+d)-d*d*exp(2*a)-d*d*exp(a)*cos(a)-exp(a)*cos(d)*d*d+2*a*d*exp(2*a)*cos(a-d)+2*a*d*exp(a)*cos(d)-2*a*exp(2*a)*d+d*d*exp(2*a)*cos(a-d)-2*a*d*exp(a)*cos(a+2*d)-2*a*exp(a)*cos(2*a+d)*d+2*a*cos(a+d)*d+a*a*exp(2*a)*cos(a+d)+2*a*d*exp(2*a)*cos(2*a+2*d)+2*a*exp(a)*cos(a)*d-a*a*exp(2*a)-a*a*exp(a)*cos(d)-a*a*exp(a)*cos(a)+2*d*d-a*a*a*a+2*a*a)/(a-d)/(a+d)/(a*a+d*d);

  values[0] = -a*a/2*(exp(2*a*x)+exp(2*a*y)+exp(2*a*z)
		      +2*sin(a*x+d*y)*cos(a*z+d*x)*exp(a*(y+z))
		      +2*sin(a*y+d*z)*cos(a*x+d*y)*exp(a*(z+x))
		      +2*sin(a*z+d*x)*cos(a*y+d*z)*exp(a*(x+y)))
      *exp(-2*eps*d*d*t) - trans;
}

void ExactU1(double x, double y,  double z, double *values)
{
  double a = TDatabase::ParamDB->P5;
  double d = TDatabase::ParamDB->P6;
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = -a*(exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y))*exp(-eps*d*d*t);
  values[1] = -a*(a*exp(a*x)*sin(a*y+d*z)-a*exp(a*z)*sin(a*x+d*y))*exp(-eps*d*d*t);
  values[2] = -a*(a*exp(a*x)*cos(a*y+d*z)-d*exp(a*z)*sin(a*x+d*y))*exp(-eps*d*d*t);
  values[3] = -a*(d*exp(a*x)*cos(a*y+d*z)+a*exp(a*z)*cos(a*x+d*y))*exp(-eps*d*d*t);
  values[4] = 0; 
}

void ExactU2(double x, double y,  double z, double *values)
{
  double a = TDatabase::ParamDB->P5;
  double d = TDatabase::ParamDB->P6;
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = -a*(exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z))*exp(-eps*d*d*t);
  values[1] = -a*(d*exp(a*y)*cos(a*z+d*x)+a*exp(a*x)*cos(a*y+d*z))*exp(-eps*d*d*t);
  values[2] = -a*(a*exp(a*y)*sin(a*z+d*x)-a*exp(a*x)*sin(a*y+d*z))*exp(-eps*d*d*t);
  values[3] = -a*(a*exp(a*y)*cos(a*z+d*x)-d*exp(a*x)*sin(a*y+d*z))*exp(-eps*d*d*t);
  values[4] = 0; 
}

void ExactU3(double x, double y,  double z, double *values)
{
  double a = TDatabase::ParamDB->P5;
  double d = TDatabase::ParamDB->P6;
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double t=TDatabase::TimeDB->CURRENTTIME;

  values[0] = -a*(exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x))*exp(-eps*d*d*t);
  values[1] = -a*(a*exp(a*z)*cos(a*x+d*y)-d*exp(a*y)*sin(a*z+d*x))*exp(-eps*d*d*t);
  values[2] = -a*(d*exp(a*z)*cos(a*x+d*y)+a*exp(a*y)*cos(a*z+d*x))*exp(-eps*d*d*t);
  values[3] = -a*(a*exp(a*z)*sin(a*x+d*y)-a*exp(a*y)*sin(a*z+d*x))*exp(-eps*d*d*t);

  values[4] = 0; 

}

void ExactP(double x, double y,  double z, double *values)
{
  double a = TDatabase::ParamDB->P5;
  double d = TDatabase::ParamDB->P6;
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double t=TDatabase::TimeDB->CURRENTTIME;
  
  double trans = -0.75*exp(-2*eps*d*d*t)*a*(-d*d*exp(2*a)*sin(a+d)-a*a*cos(a-d)+exp(a)*cos(2*a+d)*d*d+d*d*exp(2*a)*cos(a+d)+d*d*exp(2*a)*sin(2*a+2*d)-d*d*exp(a)*sin(d)+a*a*exp(a)*cos(a+2*d)+a*a*exp(a)*cos(2*a+d)-cos(a-d)*d*d-cos(a+d)*d*d+a*a*exp(2*a)*cos(a-d)+3*d*d*exp(a)*sin(a)-d*d*exp(a)*sin(2*a+d)-a*a*exp(2*a)*sin(2*a+2*d)-exp(a)*sin(a+2*d)*d*d+a*a*exp(2*a)*sin(a-d)+a*a*exp(a)*sin(2*a+d)+exp(a)*cos(a+2*d)*d*d+a*a*exp(a)*sin(a+2*d)-3*a*a*exp(a)*sin(a)+a*a*exp(2*a)*sin(a+d)-sin(a-d)*d*d-a*a*sin(a+d)+sin(a+d)*d*d+a*a*sin(a-d)+a*a*exp(a)*sin(d)-d*d*exp(2*a)*sin(a-d)+d*d*d*d-exp(2*a)*d*d*cos(2*a+2*d)-2*a*d*exp(2*a)*cos(a+d)-exp(2*a)*d*d*d*d+exp(2*a)*a*a*a*a-2*a*cos(a-d)*d-a*a*exp(2*a)*cos(2*a+2*d)-a*a*cos(a+d)-d*d*exp(2*a)-d*d*exp(a)*cos(a)-exp(a)*cos(d)*d*d+2*a*d*exp(2*a)*cos(a-d)+2*a*d*exp(a)*cos(d)-2*a*exp(2*a)*d+d*d*exp(2*a)*cos(a-d)-2*a*d*exp(a)*cos(a+2*d)-2*a*exp(a)*cos(2*a+d)*d+2*a*cos(a+d)*d+a*a*exp(2*a)*cos(a+d)+2*a*d*exp(2*a)*cos(2*a+2*d)+2*a*exp(a)*cos(a)*d-a*a*exp(2*a)-a*a*exp(a)*cos(d)-a*a*exp(a)*cos(a)+2*d*d-a*a*a*a+2*a*a)/(a-d)/(a+d)/(a*a+d*d);

  values[0] = -a*a/2*(exp(2*a*x)+exp(2*a*y)+exp(2*a*z)
		      +2*sin(a*x+d*y)*cos(a*z+d*x)*exp(a*(y+z))
		      +2*sin(a*y+d*z)*cos(a*x+d*y)*exp(a*(z+x))
		      +2*sin(a*z+d*x)*cos(a*y+d*z)*exp(a*(x+y)))
      *exp(-2*eps*d*d*t) - trans;
  values[1] = -1/2*a*a*(2*a*exp(2*a*x)+2*cos(a*x+d*y)*a*cos(a*z+d*x)*exp(a*(y+z))-2*sin(a*x+d*y)*sin(a*z+d*x)*d*exp(a*(y+z))-2*sin(a*y+d*z)*sin(a*x+d*y)*a*exp(a*(z+x))+2*sin(a*y+d*z)*cos(a*x+d*y)*a*exp(a*(z+x))+2*cos(a*z+d*x)*d*cos(a*y+d*z)*exp(a*(x+y))+2*sin(a*z+d*x)*cos(a*y+d*z)*a*exp(a*(x+y)))*exp(-2*eps*d*d*t);
  values[2] = -1/2*a*a*(2*a*exp(2*a*y)+2*cos(a*x+d*y)*d*cos(a*z+d*x)*exp(a*(y+z))+2*sin(a*x+d*y)*cos(a*z+d*x)*a*exp(a*(y+z))+2*cos(a*y+d*z)*a*cos(a*x+d*y)*exp(a*(z+x))-2*sin(a*y+d*z)*sin(a*x+d*y)*d*exp(a*(z+x))-2*sin(a*z+d*x)*sin(a*y+d*z)*a*exp(a*(x+y))+2*sin(a*z+d*x)*cos(a*y+d*z)*a*exp(a*(x+y)))*exp(-2*eps*d*d*t);
  values[3] = -1/2*a*a*(2*a*exp(2*a*z)-2*sin(a*x+d*y)*sin(a*z+d*x)*a*exp(a*(y+z))+2*sin(a*x+d*y)*cos(a*z+d*x)*a*exp(a*(y+z))+2*cos(a*y+d*z)*d*cos(a*x+d*y)*exp(a*(z+x))+2*sin(a*y+d*z)*cos(a*x+d*y)*a*exp(a*(z+x))+2*cos(a*z+d*x)*a*cos(a*y+d*z)*exp(a*(x+y))-2*sin(a*z+d*x)*sin(a*y+d*z)*d*exp(a*(x+y)))*exp(-2*eps*d*d*t);
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  double a = TDatabase::ParamDB->P5;
  double d = TDatabase::ParamDB->P6;
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double t=TDatabase::TimeDB->CURRENTTIME;

  value = -a*(exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y))*exp(-eps*d*d*t);
  
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  double a = TDatabase::ParamDB->P5;
  double d = TDatabase::ParamDB->P6;
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double t=TDatabase::TimeDB->CURRENTTIME;

  value = -a*(exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z))*exp(-eps*d*d*t);   
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
  double a = TDatabase::ParamDB->P5;
  double d = TDatabase::ParamDB->P6;
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double t=TDatabase::TimeDB->CURRENTTIME;

  value = -a*(exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x))*exp(-eps*d*d*t);
  
}

// value of boundary condition
void U1BoundValue_diff(double x, double y, double z, double &value)
{
  double a = TDatabase::ParamDB->P5;
  double d = TDatabase::ParamDB->P6;
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double t=TDatabase::TimeDB->CURRENTTIME;

  value = eps*d*d*a*(exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y))*exp(-eps*d*d*t);
  
}

// value of boundary condition
void U2BoundValue_diff(double x, double y, double z, double &value)
{
  double a = TDatabase::ParamDB->P5;
  double d = TDatabase::ParamDB->P6;
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double t=TDatabase::TimeDB->CURRENTTIME;

  value = eps*d*d*a*(exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z))*exp(-eps*d*d*t);   
}

// value of boundary condition
void U3BoundValue_diff(double x, double y, double z, double &value)
{
  double a = TDatabase::ParamDB->P5;
  double d = TDatabase::ParamDB->P6;
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double t=TDatabase::TimeDB->CURRENTTIME;

  value = eps*d*d*a*(exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x))*exp(-eps*d*d*t);
  
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y, z, u1, u2, u3, ux, uy, uz;

  if (TDatabase::ParamDB->STOKES_PROBLEM)
  {
      OutPut("Stokes case not implemented !!!"<<endl);
      exit(4711);
  }
  else
  {
    for(i=0;i<n_points;i++)
    {
      coeff = coeffs[i];
      
      coeff[0] = eps;
      coeff[1] = 0;
      coeff[2] = 0;
      coeff[3] = 0;
      
      coeff[4] = 0;
      coeff[5] = 0;
      coeff[6] = 0;
      
    }
  }
    
}
