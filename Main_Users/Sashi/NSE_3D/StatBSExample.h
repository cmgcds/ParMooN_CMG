// Time Dependent Navier-Stokes problem
// for paper with Traian Iliescu and Bill Layton
//
#include <TNSE3D_Routines.h>

void ExampleFile()
{
  OutPut("Example: StatBSExample.h" << endl) ;
}

void InitialU1(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  fac = sin(t*Pi/2);
  fac = 1.0;

  values[0] = sin(Pi*x)*sin(Pi*y)*sin(Pi*z)+x*x*x*x*cos(Pi*y);
  values[0] *=fac ;
}

void InitialU2(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  fac = sin(t*Pi/2);
  fac = 1.0;

  values[0] =  cos(Pi*x)*cos(Pi*y)*cos(Pi*z)-3*y*y*y*z;
  values[0] *= fac;
}

void InitialU3(double x, double y, double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  fac = sin(t*Pi/2);
  fac = 1.0;

  values[0] = cos(Pi*x)*sin(Pi*y)*cos(Pi*z)+cos(Pi*x)*sin(Pi*y)*sin(Pi*z)
    -4*x*x*x*cos(Pi*y)*z+4.5*y*y*z*z;
  values[0] *= fac;
}
void InitialP(double x, double y, double z, double *values)
{
  values[0] = 0;
}
// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  fac = sin(t*Pi/2);
  fac = 1.0;

  values[0] = sin(Pi*x)*sin(Pi*y)*sin(Pi*z)+x*x*x*x*cos(Pi*y);
  values[1] = cos(Pi*x)*Pi*sin(Pi*y)*sin(Pi*z)+4*x*x*x*cos(Pi*y);
  values[2] = sin(Pi*x)*cos(Pi*y)*Pi*sin(Pi*z)-x*x*x*x*sin(Pi*y)*Pi;
  values[3] = sin(Pi*x)*sin(Pi*y)*cos(Pi*z)*Pi;
  values[4] = -3*sin(Pi*x)*Pi*Pi*sin(Pi*y)*sin(Pi*z)
                       +12*x*x*cos(Pi*y)-x*x*x*x*cos(Pi*y)*Pi*Pi;
  values[0] *= fac;
  values[1] *= fac;
  values[2] *= fac;
  values[3] *= fac;
  values[4] *= fac;
}

void ExactU2(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  fac = sin(t*Pi/2);
  fac = 1.0;

  values[0] =  cos(Pi*x)*cos(Pi*y)*cos(Pi*z)-3*y*y*y*z;
  values[1] = -Pi*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
  values[2] = -Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z)-9*y*y*z;
  values[3] = -Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z)-3*y*y*y;
  values[4] = -3*cos(Pi*x)*Pi*Pi*cos(Pi*y)*cos(Pi*z)-18*y*z;
  values[0] *= fac;
  values[1] *= fac;
  values[2] *= fac;
  values[3] *= fac;
  values[4] *= fac;
}

void ExactU3(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  fac = sin(t*Pi/2);
  fac = 1.0;

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
  values[0] *= fac;
  values[1] *= fac;
  values[2] *= fac;
  values[3] *= fac;
  values[4] *= fac;
}

void ExactP(double x, double y,  double z, double *values)
{

  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  fac = sin(t*Pi/2);
  fac = 1.0;

  values[0] = 3*x-sin(y+4*z)-1.5-sin(1.0)*cos(1.0)
    -4*cos(1.0)*cos(1.0)*cos(1.0)*cos(1.0)*sin(1.0)
    +sin(1.0)*cos(1.0)*cos(1.0)-2*sin(1.0)*sin(1.0)*sin(1.0)
    +2*sin(1.0)*cos(1.0)*cos(1.0)*cos(1.0)+2*sin(1.0);
  values[1] = 3;
  values[2] = -cos(y+4*z);
  values[3] = -4*cos(y+4*z);
  values[4] = 0;
  values[0] *= fac;
  values[1] *= fac;
  values[2] *= fac;
  values[3] *= fac;
  values[4] *= fac;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void U1BoundValue(double x, double y, double z, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  value = sin(Pi*x)*sin(Pi*y)*sin(Pi*z)+x*x*x*x*cos(Pi*y);
  //value *=sin(t*Pi/2);
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  value = cos(Pi*x)*cos(Pi*y)*cos(Pi*z)-3*y*y*y*z;
  //value *=sin(t*Pi/2);
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  value = cos(Pi*x)*sin(Pi*y)*cos(Pi*z)+cos(Pi*x)*sin(Pi*y)*sin(Pi*z)
    -4*x*x*x*cos(Pi*y)*z+4.5*y*y*z*z;
  //value *=sin(t*Pi/2);
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  double eps = 1/TDatabase::ParamDB->RE_NR;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double gamma =  TDatabase::ParamDB->GAUSSIAN_GAMMA;
  int i;
  double *coeff, x, y, z, u1, u2, u3, p_x, p_y, p_z;
  double u1_x, u1_y, u1_z, u2_x, u2_y, u2_z, u3_x, u3_y, u3_z; 
  double u1_xx, u1_yy, u1_zz, u2_xx, u2_yy, u2_zz, u3_xx, u3_yy, u3_zz; 
  double u1_xy, u1_xz, u1_yz, u2_xy, u2_xz, u2_yz, u3_xy, u3_xz, u3_yz; 
  double u1_t, u2_t, u3_t, fact, fact_t, c, tmp, delta;
  double div_1, div_2, div_3, taylor_1, taylor_2, taylor_3;
  double  nu_T,  nu_T_x,  nu_T_y,  nu_T_z, grad_u[9];

  // compute characteristic filter width
  // it is assumed that the filter width does not depend on h
  // such that tmp=1 is just a dummy input
  tmp = 1;
  delta =  CharacteristicFilterWidth(tmp);
  c = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT*delta*delta;
  fact = sin(t*Pi/2);
  fact_t = Pi*cos(t*Pi/2)/2.0;
  fact = 1.0;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    // coordinates
    x = X[i];
    y = Y[i];
    z = Z[i];
    
    //functions
    u1 = sin(Pi*x)*sin(Pi*y)*sin(Pi*z)+x*x*x*x*cos(Pi*y);
    u2 =  cos(Pi*x)*cos(Pi*y)*cos(Pi*z)-3*y*y*y*z;
    u3 = cos(Pi*x)*sin(Pi*y)*cos(Pi*z)+cos(Pi*x)*sin(Pi*y)*sin(Pi*z)
      -4*x*x*x*cos(Pi*y)*z+4.5*y*y*z*z;
    u1 *= fact;
    u2 *= fact;
    u3 *= fact;

    u1_t = sin(Pi*x)*sin(Pi*y)*sin(Pi*z)+x*x*x*x*cos(Pi*y);
    u1_t *= fact_t;
    u1_x =  cos(Pi*x)*Pi*sin(Pi*y)*sin(Pi*z)+4*x*x*x*cos(Pi*y);
    u1_y = sin(Pi*x)*cos(Pi*y)*Pi*sin(Pi*z)-x*x*x*x*sin(Pi*y)*Pi;
    u1_z = sin(Pi*x)*sin(Pi*y)*cos(Pi*z)*Pi;
    u1_x *= fact;
    u1_y *= fact;
    u1_z *= fact;

    u2_t =  cos(Pi*x)*cos(Pi*y)*cos(Pi*z)-3*y*y*y*z;
    u2_t *= fact_t;
    u2_x = -Pi*sin(Pi*x)*cos(Pi*y)*cos(Pi*z);
    u2_y = -Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z)-9*y*y*z;
    u2_z = -Pi*cos(Pi*x)*cos(Pi*y)*sin(Pi*z)-3*y*y*y;      
    u2_x *= fact;
    u2_y *= fact;
    u2_z *= fact;

    u3_t =  cos(Pi*x)*sin(Pi*y)*cos(Pi*z)+cos(Pi*x)*sin(Pi*y)*sin(Pi*z)
      -4*x*x*x*cos(Pi*y)*z+4.5*y*y*z*z;
    u3_t *=  fact_t;
    u3_x = -sin(Pi*x)*sin(Pi*y)*cos(Pi*z)*Pi
      -12*x*x*cos(Pi*y)*z-sin(Pi*z)*sin(Pi*x)*Pi*sin(Pi*y);
    u3_y = cos(Pi*z)*cos(Pi*x)*cos(Pi*y)*Pi
      +4*x*x*x*sin(Pi*y)*Pi*z+cos(Pi*x)*cos(Pi*y)*sin(Pi*z)*Pi+9*y*z*z;
    u3_z = -cos(Pi*x)*Pi*sin(Pi*y)*sin(Pi*z)
      -4*x*x*x*cos(Pi*y)+cos(Pi*x)*sin(Pi*y)*Pi*cos(Pi*z)+9*y*y*z;
    u3_x *= fact;
    u3_y *= fact;
    u3_z *= fact;
 
    u1_xx = -Pi*sin(Pi*x)*Pi*sin(Pi*y)*sin(Pi*z)+12*x*x*cos(Pi*y);
    u1_yy = -Pi*sin(Pi*x)*sin(Pi*y)*Pi*sin(Pi*z)-x*x*x*x*cos(Pi*y)*Pi*Pi; 
    u1_zz = -sin(Pi*x)*sin(Pi*y)*sin(Pi*z)*Pi*Pi;
    u1_xy = cos(Pi*x)*Pi*Pi*cos(Pi*y)*sin(Pi*z)-4*x*x*x*sin(Pi*y)*Pi;
    u1_xz = cos(Pi*x)*Pi*sin(Pi*y)*cos(Pi*z)*Pi;
    u1_yz = sin(Pi*x)*cos(Pi*y)*Pi*cos(Pi*z)*Pi;
    u1_xx *= fact;
    u1_yy *= fact;
    u1_zz *= fact;
    u1_xy *= fact;
    u1_xz *= fact;
    u1_yz *= fact;

    u2_xx = -Pi*Pi*cos(Pi*x)*cos(Pi*y)*cos(Pi*z);
    u2_yy = -Pi*cos(Pi*x)*Pi*cos(Pi*y)*cos(Pi*z)-18*y*z;
    u2_zz = -Pi*cos(Pi*x)*cos(Pi*y)*Pi*cos(Pi*z);
    u2_xy = Pi*sin(Pi*x)*Pi*sin(Pi*y)*cos(Pi*z);
    u2_xz = Pi*sin(Pi*x)*cos(Pi*y)*Pi*sin(Pi*z);
    u2_yz = Pi*cos(Pi*x)*sin(Pi*y)*Pi*sin(Pi*z)-9*y*y;
    u2_xx *= fact;
    u2_yy *= fact;
    u2_zz *= fact;
    u2_xy *= fact;
    u2_xz *= fact;
    u2_yz *= fact;

    u3_xx = -Pi*cos(Pi*x)*sin(Pi*y)*cos(Pi*z)*Pi
      -24*x*cos(Pi*y)*z-sin(Pi*z)*Pi*cos(Pi*x)*Pi*sin(Pi*y);
    u3_yy = -cos(Pi*z)*cos(Pi*x)*Pi*sin(Pi*y)*Pi
      +4*x*x*x*Pi*cos(Pi*y)*Pi*z-cos(Pi*x)*Pi*sin(Pi*y)*sin(Pi*z)*Pi+9*z*z;
    u3_zz =  -cos(Pi*x)*Pi*sin(Pi*y)*Pi*cos(Pi*z)
      -cos(Pi*x)*sin(Pi*y)*Pi*Pi*sin(Pi*z)+9*y*y;
    u3_xy =  -sin(Pi*x)*Pi*cos(Pi*y)*cos(Pi*z)*Pi
      +12*x*x*Pi*sin(Pi*y)*z-sin(Pi*z)*sin(Pi*x)*Pi*Pi*cos(Pi*y);
    u3_xz = sin(Pi*x)*sin(Pi*y)*Pi*sin(Pi*z)*Pi
      -12*x*x*cos(Pi*y)-Pi*cos(Pi*z)*sin(Pi*x)*Pi*sin(Pi*y);
    u3_yz =  -Pi*sin(Pi*z)*cos(Pi*x)*cos(Pi*y)*Pi
      +4*x*x*x*sin(Pi*y)*Pi+cos(Pi*x)*cos(Pi*y)*Pi*cos(Pi*z)*Pi+18*y*z;
    u3_xx *= fact;
    u3_yy *= fact;
    u3_zz *= fact;
    u3_xy *= fact;
    u3_xz *= fact;
    u3_yz *= fact;

    
// turbulent viscosity (Smagorinsky term)    
    grad_u[0] = u1_x;
    grad_u[1] = u2_x;
    grad_u[2] = u3_x;
    grad_u[3] = u1_y;
    grad_u[4] = u2_y;
    grad_u[5] = u3_y;
    grad_u[6] = u1_z;
    grad_u[7] = u2_z;
    grad_u[8] = u3_z;
    nu_T = TurbulentViscosity3D(delta, grad_u, grad_u, grad_u);
    if (nu_T == 0)
    {
      nu_T_x = 0;
      nu_T_y = 0;
      nu_T_z = 0;
    }
    else
    {
      nu_T_x = u1_x*u1_xx+u1_y*u1_xy+ u1_z*u1_xz;
      nu_T_x += u2_x*u2_xx+u2_y*u2_xy + u2_z*u2_xz;
      nu_T_x += u3_x*u3_xx+u3_y*u3_xy + u3_z*u3_xz;
      nu_T_x *= c*c/nu_T;
      nu_T_y = u1_x*u1_xy+u1_y*u1_yy+ u1_z*u1_yz;
      nu_T_y += u2_x*u2_xy+u2_y*u2_yy + u2_z*u2_yz;
      nu_T_y += u3_x*u3_xy+u3_y*u3_yy + u3_z*u3_yz;
      nu_T_y *= c*c/nu_T;
      nu_T_z = u1_x*u1_xz+u1_y*u1_yz+ u1_z*u1_zz;
      nu_T_z += u2_x*u2_xz+u2_y*u2_yz + u2_z*u2_zz;
      nu_T_z += u3_x*u3_xz+u3_y*u3_yz + u3_z*u3_zz;
      nu_T_z *= c*c/nu_T;
    }

    taylor_1 =  2 *u1_x*u1_xx + 2 * u1_y*u1_xy + 2*u1_z*u1_xz;
    taylor_1 += u1_xy * u2_x + u1_x * u2_xy + u1_yy * u2_y + u1_y * u2_yy 
      + u1_yz*u2_z + u1_z * u2_yz;
    taylor_1 += u1_xz*u3_x + u1_x * u3_xz + u1_yz*u3_y + u1_y * u3_yz 
      + u1_zz*u3_z + u1_z * u3_zz;
    taylor_2 =  u1_xx * u2_x + u1_x * u2_xx +  u1_xy * u2_y+  u1_y * u2_xy
      + u1_xz * u2_z + u1_z * u2_xz;
    taylor_2 += 2 *u2_x*u2_xy + 2*u2_y*u2_yy + 2* u2_z*u2_yz;
    taylor_2 += u2_xz * u3_x + u2_x *u3_xz + u2_yz*u3_y + u2_y*u3_yz
      + u2_zz*u3_z + u2_z * u3_zz;
    taylor_3 = u1_xx+u3_x + u1_x *u3_xx+ u1_xy*u3_y + u1_y * u3_xy
      + u1_xz *u3_z + u1_z * u3_xz;
    taylor_3 += u2_xy*u3_x + u2_x * u3_xy + u2_yy * u3_y + u2_y * u3_yy
      + u2_yz * u3_z + u2_z * u3_yz;
    taylor_3 += 2 * u3_x * u3_xz + 2* u3_y * u3_yz + 2 * u3_z * u3_zz;
    
// divergence of the Laplacian
    div_1 = u1_xx + u1_yy + u1_zz;
    div_2 = u2_xx + u2_yy + u2_zz;
    div_3 = u3_xx + u3_yy + u3_zz;

    // pressure
    p_x = 3;
    p_y = -cos(y+4*z);
    p_z = -4*cos(y+4*z);

    coeff[0] = eps;

    taylor_1 = taylor_2 = taylor_3 = 0.0;
    //nu_T_x = nu_T_y = nu_T_z = 0.0;
    coeff[1] = u1_t - (eps + nu_T) *div_1 
      -(u1_x*nu_T_x + u1_y*nu_T_y + u1_z * nu_T_z)
      + u1 * u1_x + u2 * u1_y + u3 * u1_z +
      p_x + delta * delta /(2*gamma) * taylor_1;
    coeff[2] = u2_t -  (eps + nu_T)*div_2 
      -(u2_x*nu_T_x + u2_y*nu_T_y + u2_z * nu_T_z)
      + u1 * u2_x + u2 * u2_y + u3 * u2_z
      + p_y + delta * delta /(2*gamma) * taylor_2;
    coeff[3] = u3_t -  (eps + nu_T)*div_3 
      -(u3_x*nu_T_x + u3_y*nu_T_y + u3_z * nu_T_z)
      + u1 * u3_x + u2 * u3_y + u3 * u3_z
      + p_z + delta * delta /(2*gamma) * taylor_3;

    coeff[1] = - eps *div_1 + u1 * u1_x + u2 * u1_y + u3 * u1_z + p_x;
    coeff[1] += -nu_T *div_1 -(u1_x*nu_T_x + u1_y*nu_T_y + u1_z * nu_T_z);  
    coeff[2] = - eps *div_2 + u1 * u2_x + u2 * u2_y + u3 * u2_z + p_y;
    coeff[2] += - nu_T*div_2 -(u2_x*nu_T_x + u2_y*nu_T_y + u2_z * nu_T_z);  
    coeff[3] = - eps *div_3 + u1 * u3_x + u2 * u3_y + u3 * u3_z + p_z;
    coeff[3] += -nu_T*div_3 -(u3_x*nu_T_x + u3_y*nu_T_y + u3_z * nu_T_z);  
  }
}

