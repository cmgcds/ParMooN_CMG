// Navier-Stokes problem, solution in ansatz space
// velocity pw quadratic, pressure linear
// 
// u(x,y) = t(x^2+y^2+z^2+y^5, x^2+2xy+13+3*z^4, -2xz+5y^2-x^4y)^T
// p(x,y) = 2t*(3x-2y+7z-4)

#include <TNSE3D_Routines.h>

void ExampleFile()
{
  OutPut("Example: Polynom_Taylor.h" << endl);
}

// ========================================================================
// initial condition
// ========================================================================
void InitialU1(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*(x*x+y*y+z*z+y*y*y*y*y);
}

void InitialU2(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*(x*x+2*x*z+13+3*z*z*z*z);
}

void InitialU3(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  values[0] = t*(-2*x*z+5*y*y-x*x*x*x*y);
}

void InitialP(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  values[0] = 2*t*(3*x-2*y+7*z-4);
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  values[0] = x*x+y*y+z*z+y*y*y*y*y;
  values[1] = 2*x;
  values[2] = 2*y+5*y*y*y*y;
  values[3] = 2*z;
  values[4] = 6;
  values[0] *= t;
  values[1] *= t;
  values[2] *= t;
  values[3] *= t;
  values[4] *= t;
}

void ExactU2(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  values[0] = x*x+2*x*z+13+3*z*z*z*z;
  values[1] = 2*x+2*z;
  values[2] = 0;
  values[3] = 2*x+12*z*z*z;
  values[4] = 2;
  values[0] *= t;
  values[1] *= t;
  values[2] *= t;
  values[3] *= t;
  values[4] *= t;
}

void ExactU3(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  values[0] = -2*x*z+5*y*y-x*x*x*x*y;
  values[1] = -2*z-4*x*x*x*y;
  values[2] = 10*y-x*x*x*x;
  values[3] = -2*x;
  values[4] = 10;
  values[0] *= t;
  values[1] *= t;
  values[2] *= t;
  values[3] *= t;
  values[4] *= t;
}

void ExactP(double x, double y,  double z, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  values[0] = 3*x-2*y+7*z-4;
  values[1] = 3;
  values[2] = -2;
  values[3] = 7;
  values[4] = 0;
  values[0] *= 2*t;
  values[1] *= 2*t;
  values[2] *= 2*t;
  values[3] *= 2*t;
  values[4] *= 2*t;
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
  value = x*x+y*y+z*z+y*y*y*y*y;
  value *= t;
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  value = x*x+2*x*z+13+3*z*z*z*z;
  value *= t;
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  value = -2*x*z+5*y*y-x*x*x*x*y;
  value *= t;
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
  fact = t;
  fact_t = 1;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    // coordinates
    x = X[i];
    y = Y[i];
    z = Z[i];
    
    //functions
    u1 = x*x+y*y+z*z+y*y*y*y*y;
    u2 = x*x+2*x*z+13+3*z*z*z*z ;
    u3 = -2*x*z+5*y*y-x*x*x*x*y;
    u1 *= fact;
    u2 *= fact;
    u3 *= fact;

    u1_t = x*x+y*y+z*z+y*y*y*y*y;
    u1_t *= fact_t;
    u1_x =  2*x;
    u1_y =  2*y+5*y*y*y*y;
    u1_z =  2*z;
    u1_x *= fact;
    u1_y *= fact;
    u1_z *= fact;

    u2_t =  x*x+2*x*z+13+3*z*z*z*z;
    u2_t *= fact_t;
    u2_x = 2*x+2*z;
    u2_y = 0;
    u2_z = 2*x+12*z*z*z;      
    u2_x *= fact;
    u2_y *= fact;
    u2_z *= fact;

    u3_t =  -2*x*z+5*y*y-x*x*x*x*y;
    u3_t *=  fact_t;
    u3_x = -2*z-4*x*x*x*y;
    u3_y = 10*y-x*x*x*x;
    u3_z = -2*x;
    u3_x *= fact;
    u3_y *= fact;
    u3_z *= fact;
 
    u1_xx = 2;
    u1_yy = 2+20*y*y*y; 
    u1_zz = 2;
    u1_xy = 0;
    u1_xz = 0;
    u1_yz = 0;
    u1_xx *= fact;
    u1_yy *= fact;
    u1_zz *= fact;
    u1_xy *= fact;
    u1_xz *= fact;
    u1_yz *= fact;

    u2_xx = 2;
    u2_yy = 0;
    u2_zz = 36*z*z;
    u2_xy = 0;
    u2_xz = 2;
    u2_yz = 0;
    u2_xx *= fact;
    u2_yy *= fact;
    u2_zz *= fact;
    u2_xy *= fact;
    u2_xz *= fact;
    u2_yz *= fact;

    u3_xx = -12*x*x*y;
    u3_yy = 10;
    u3_zz =  0;
    u3_xy = -4*x*x*x;
    u3_xz = -2;
    u3_yz =  0;
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
    p_x = 3*2*t;
    p_y = -2*2*t;
    p_z = 7*2*t;

    coeff[0] = eps;

    //taylor_1 = taylor_2 = taylor_3 = 0.0;
    //nu_T= nu_T_x = nu_T_y = nu_T_z = 0.0;
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

  }
}

