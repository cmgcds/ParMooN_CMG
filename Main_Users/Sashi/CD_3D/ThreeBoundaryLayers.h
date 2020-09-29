// ======================================================================
// three boundary layer problem
// ======================================================================
#include <ConvDiff3D.h>

void ExampleFile()
{
  OutPut("Example: ThreeBoundaryLayers.h" << endl) ;
}

// exact solution
void Exact(double x, double y, double z, double *values)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;

  double t2, t4, t5, t6, t9, t10, t11, t12, t16, t17, t30, t31;
  t2 = 1/eps;
  t4 = exp(2.0*(x-1.0)*t2);
  t5 = x-t4;
  t6 = y*y;
  t9 = exp(3.0*(y-1.0)*t2);
  t10 = t6-t9;
  t11 = t5*t10;
  t12 = z*z;
  t16 = exp(4.0*(z-1.0)*t2);
  t17 = t12*z-t16;
  t30 = eps*eps;
  t31 = 1/t30;
  values[0] = t11*t17;
  values[1] = (1.0-2.0*t2*t4)*t10*t17;
  values[2] = t5*(2.0*y-3.0*t2*t9)*t17;
  values[3] = t11*(3.0*t12-4.0*t2*t16);
  values[4] = -4.0*t31*t4*t10*t17+
    t5*(2.0-9.0*t31*t9)*t17+
    t11*(6.0*z-16.0*t31*t16);
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;

  double t2, t6, t12;
  t2 = 1/eps;
  t6 = y*y;
  t12 = z*z;

  value = (x-exp(2.0*(x-1.0)*t2))*
    (t6-exp(3.0*(y-1.0)*t2))*
    (t12*z-exp(4.0*(z-1.0)*t2));
}

void BilinearCoeffs(int n_points, double *X, double *Y, double *Z,
        double **parameters, double **coeffs)
{
  static double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double t1, t2, t4, t6, t8, t11, t12, t13, t17, t18, t21, t26, t44;
  double *coeff, *param;
  double x, y, z;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];
    x=X[i];
    y=Y[i];
    z=Z[i];

    t1 = eps*eps;
    t2 = 1/t1;
    t4 = 1/eps;
    t6 = exp(2.0*(x-1.0)*t4);
    t8 = y*y;
    t11 = exp(3.0*(y-1.0)*t4);
    t12 = t8-t11;
    t13 = z*z;
    t17 = exp(4.0*(z-1.0)*t4);
    t18 = t13*z-t17;
    t21 = x-t6;
    t26 = t21*t12;
    t44 = -eps*(-4.0*t2*t6*t12*t18+t21*(2.0-9.0*t2*t11)*t18+
                t26*(6.0*z-16.0*t2*t17))
         +2.0*(1.0-2.0*t4*t6)*t12*t18
         +3.0*t21*(2.0*y-3.0*t4*t11)*t18
         +4.0*t26*(3.0*t13-4.0*t4*t17)
         +t26*t18;

    coeff[0] = eps;
    coeff[1] = 2;
    coeff[2] = 3;
    coeff[3] = 4;
    coeff[4] = 1;
    coeff[5] = t44;
  }
}
