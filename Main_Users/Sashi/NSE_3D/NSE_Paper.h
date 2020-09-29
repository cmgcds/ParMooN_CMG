// Navier-Stokes problem, solution in ansatz space
// velocity pw quadratic, pressure linear
//
//

void ExampleFile()
{
  OutPut("Example: NSE_Paper.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = sin(Pi*x)*sin(Pi*y)*sin(Pi*z)+x*x*x*x*cos(Pi*y);
  values[1] = cos(Pi*x)*Pi*sin(Pi*y)*sin(Pi*z)+4*x*x*x*cos(Pi*y);
  values[2] = sin(Pi*x)*cos(Pi*y)*Pi*sin(Pi*z)-x*x*x*x*sin(Pi*y)*Pi;
  values[3] = sin(Pi*x)*sin(Pi*y)*cos(Pi*z)*Pi;
  values[4] = -3*sin(Pi*x)*Pi*Pi*sin(Pi*y)*sin(Pi*z)
                       +12*x*x*cos(Pi*y)-x*x*x*x*cos(Pi*y)*Pi*Pi;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] =  cos(Pi*x)*cos(Pi*y)*cos(Pi*z)-3*y*y*y*z;
  values[1] = -Pi*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[2] = -Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z)-9*y*y*z;
  values[3] = -Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z)-3*y*y*y;
  values[4] = -3*cos(Pi*x)*Pi*Pi*cos(Pi*y)*cos(Pi*z)-18*y*z;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = cos(Pi*x)*sin(Pi*y)*cos(Pi*z)+cos(Pi*x)*sin(Pi*y)*sin(Pi*z)
    -4*x*x*x*cos(Pi*y)*z+4.5*y*y*z*z;
  values[1] = -sin(Pi*x)*sin(Pi*y)*cos(Pi*z)*Pi
    -12*x*x*cos(Pi*y)*z-sin(Pi*z)*sin(Pi*x)*Pi*sin(Pi*y);
  values[2] = cos(Pi*z)*cos(Pi*x)*cos(Pi*y)*Pi
    +4*x*x*x*sin(Pi*y)*Pi*z+cos(Pi*x)*cos(Pi*y)*sin(Pi*z)*Pi+9*y*z*z;
  values[3] = -cos(Pi*x)*Pi*sin(Pi*y)*sin(Pi*z)
    -4*x*x*x*cos(Pi*y)+cos(Pi*x)*sin(Pi*y)*Pi*cos(Pi*z)+9*y*y*z;
  values[4] = -3*cos(Pi*x)*Pi*Pi*sin(Pi*y)*cos(Pi*z)
                       -24*x*cos(Pi*y)*z
                       -3*sin(Pi*z)*cos(Pi*x)*Pi*Pi*sin(Pi*y)
                       +4*x*x*x*cos(Pi*y)*Pi*Pi*z
                       +9*z*z+9*y*y;
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 3*x-sin(y+4*z)-1.5-0.25*sin(5.0)+0.25*sin(4.0)+0.25*sin(1.0);
  values[1] = 3;
  values[2] = -cos(y+4*z);
  values[3] = -4*cos(y+4*z);
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
  value = sin(Pi*x)*sin(Pi*y)*sin(Pi*z)+x*x*x*x*cos(Pi*y);
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  value = cos(Pi*x)*cos(Pi*y)*cos(Pi*z)-3*y*y*y*z;
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
  value = cos(Pi*x)*sin(Pi*y)*cos(Pi*z)+cos(Pi*x)*sin(Pi*y)*sin(Pi*z)
    -4*x*x*x*cos(Pi*y)*z+4.5*y*y*z*z;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y, z; 
  double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
  double t13,t15,t18,t20,t22,t23,t25,t26,t28,t29,t31;
  double t34,t36,t37,t42,t43,t46,t48,t49,t50,t52,t53;
  double t59,t60,t62,t64,t69,t73,t88,t91,t92,t95,t99;
  double t103,t113,t116,t117,t122,t123,t125,t154,t160;
  double t203;

  if (TDatabase::ParamDB->STOKES_PROBLEM)
  {
    for(i=0;i<n_points;i++)
    {
      coeff = coeffs[i];
      
      x = X[i];
      y = Y[i];
      z = Z[i];
      coeff[0] = eps;
      coeff[1] = -eps*(-3*sin(Pi*x)*Pi*Pi*sin(Pi*y)*sin(Pi*z)
	+12*x*x*cos(Pi*y)-x*x*x*x*cos(Pi*y)*Pi*Pi)+3.0; // f1
      coeff[2] = -eps*(-3*cos(Pi*x)*Pi*Pi*cos(Pi*y)*cos(Pi*z)-18*y*z) 
        - cos(y+4*z); // f2
      coeff[3] = -eps*(-3*cos(Pi*x)*Pi*Pi*sin(Pi*y)*cos(Pi*z)
                       -24*x*cos(Pi*y)*z
                       -3*sin(Pi*z)*cos(Pi*x)*Pi*Pi*sin(Pi*y)
                       +4*x*x*x*cos(Pi*y)*Pi*Pi*z
                       +9*z*z+9*y*y) 
        - 4*cos(y+4*z);   // f3
    }
  }
  else
  {
    for(i=0;i<n_points;i++)
    {

      coeff = coeffs[i];
      
      x = X[i];
      y = Y[i];
      z = Z[i];
      
      coeff[0] = eps;
      t1 = Pi*x;
      t2 = sin(t1);
      t4 = Pi*Pi;
      t5 = Pi*y;
      t6 = sin(t5);
      t7 = t4*t6;
      t8 = Pi*z;
      t9 = sin(t8);
      t10 = t7*t9;
      t13 = x*x;
      t15 = cos(t5);
      t18 = t13*t13;
      t20 = t15*t4;
      t22 = cos(t1);
      t23 = t2*t22;
      t25 = t15*t15;
      t26 = Pi*t25;
      t28 = t2*t6;
      t29 = t13*x;
      t31 = t9*t29*t15;
      t34 = t18*t15;
      t36 = Pi*t6;
      t37 = t36*t9;
      t42 = t22*t15;
      t43 = cos(t8);
      t46 = t18*t6*Pi;
      t48 = y*y;
      t49 = t48*y;
      t50 = t49*z;
      t52 = t15*Pi;
      t53 = t52*t9;
      t59 = Pi*t22;
      t60 = t59*t9;
      t62 = t28*t43;
      t64 = t15*z;
      t69 = z*z;
      
      coeff[1] = 3.0*eps*t2*t10-12.0*eps*t13*t15+eps*t18*t20+t23*Pi
	-t23*t26+4.0*t28*t31+t34*t22*t37+4.0*t18*t29*t25-t42*t43*t46
	-3.0*t50*t2*t53+3.0*t50*t46+t2*t43*t60-4.0*t62*Pi*t29*t64
	+4.5*t62*Pi*t48*t69+3.0; // f1
      
      t73 = eps*t22;
      t88 = t43*t48*z;
      t91 = t50*t22;
      t92 = t36*t43;
      t95 = t48*t48;
      t99 = t22*t6;
      t103 = t22*t22;
      t113 = t29*t15;
      t116 = t48*t69;
      t117 = t116*t22;
      t122 = cos(y+4.0*z);
      
      coeff[2] = 3.0*t73*t20*t43+18.0*eps*y*z-t52*t43*t6*t9
	-t2*Pi*t25*t43*t18-9.0*t42*t88+3.0*t91*t92+13.5*t95*y*t69
	-3.0*t99*t43*t49-t103*t6*t52-3.0*t99*t9*t49+4.0*t29*t25*z*t60
	+12.0*t113*t50-4.5*t117*t53-t122; // f2
      
      t123 = t43*t43;
      t125 = t9*t43;
      t154 = t34*t2;
      t160 = -t26*t123-t125*Pi+t103*t123*Pi+4.0*t18*t13*t25*z
	-54.0*t113*t116+24.0*eps*x*t64+t125*t26+t125*Pi*t103
	-4.0*t122-9.0*eps*t48+13.5*t95*t69*z-9.0*eps*t69
	+3.0*t73*t10+3.0*t73*t7*t43-t154*t92-4.0*eps*t29*t20*z;
      
      t203 = 9.0*t42*t43*y*t69-t154*t37-3.0*t91*t53-3.0*t91*t52*t43
	-4.0*t99*t43*t29*t15-12.0*t49*t69*t29*t6*Pi-4.0*t99*t31
	-12.0*t28*t9*t13*t15*z+9.0*t99*t88-4.5*t117*t37
	+4.0*t113*z*t59*t6*t9+9.0*t99*t9*t48*z+t26+4.5*t117*t92+Pi*t123-Pi;
      
      coeff[3] = t160+t203; // f3
     
    }
  }
    
}
