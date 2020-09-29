// time dependent Navier-Stokes problem with Smagorinsky
// turbulent viscosiy 
//
// Tafti, Computers & Fluids (25), p. 647 - 665, 1996
// 
// u_1(x,y) = - cos(n pi x) sin(n pi y) exp(-2 (n pi)^2 t/T)
// u_2(x,y) = sin(n pi x) cos(n pi y) exp(-2 (n pi)^2 t/T)
// p(x,y) = -1/4 (cos(2n Pi x) + cos(2n Pi y))  exp(-4 (n pi)^2 t/T)
//
// n - number of oscillations
//
// rhs for Tazlor polynomial LES model

#include <TNSE2D_Routines.h>

#define N_OSCILLATIONS    4
#define SINCOSEXP_ENDTIME    1000.0

void ExampleFile()
{
  OutPut("Example: SinCosExpT_Smago.h : N_OSCILLATIONS " << N_OSCILLATIONS ) 
  OutPut(" SINCOSEXP_ENDTIME "<< SINCOSEXP_ENDTIME << endl);
}

// ========================================================================
// initial solution
// ========================================================================

void InitialU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
 
  fac = exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*t/SINCOSEXP_ENDTIME);
  values[0] = -cos(N_OSCILLATIONS*pi*x)*sin(N_OSCILLATIONS*pi*y)*fac;
}

void InitialU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
 
  fac = exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*t/SINCOSEXP_ENDTIME);
  values[0] = sin(N_OSCILLATIONS*pi*x)*cos(N_OSCILLATIONS*pi*y)*fac;
}

void InitialP(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
 
  fac = exp(-4*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*t/SINCOSEXP_ENDTIME);
  values[0] = -(cos(2*N_OSCILLATIONS*pi*x)+cos(2*N_OSCILLATIONS*pi*y))*fac/4.0;
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
 
  fac = exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*t/SINCOSEXP_ENDTIME);
  values[0] = -cos(N_OSCILLATIONS*pi*x)*sin(N_OSCILLATIONS*pi*y)*fac;
  values[1] = N_OSCILLATIONS*pi*sin(N_OSCILLATIONS*pi*x)*sin(N_OSCILLATIONS*pi*y)*fac;
  values[2] = -N_OSCILLATIONS*pi*cos(N_OSCILLATIONS*pi*x)*cos(N_OSCILLATIONS*pi*y)*fac;
  //values[1] =  values[0] =  values[2] = 0;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
 
  fac = exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*t/SINCOSEXP_ENDTIME);
  values[0] = sin(N_OSCILLATIONS*pi*x)*cos(N_OSCILLATIONS*pi*y)*fac;
  values[1] = N_OSCILLATIONS*pi*cos(N_OSCILLATIONS*pi*x)*cos(N_OSCILLATIONS*pi*y)*fac;
  values[2] = -N_OSCILLATIONS*pi* sin(N_OSCILLATIONS*pi*x)*sin(N_OSCILLATIONS*pi*y)*fac;
  //values[1] =  values[0] =  values[2] = 0;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double pi = 3.14159265358979;
 
  fac = exp(-4*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*t/SINCOSEXP_ENDTIME);
  values[0] = -(cos(2*N_OSCILLATIONS*pi*x)+cos(2*N_OSCILLATIONS*pi*y))*fac/4.0;
  values[1] = N_OSCILLATIONS*pi*sin(2*N_OSCILLATIONS*pi*x)*fac/2.0;
  values[2] = N_OSCILLATIONS*pi*sin(2*N_OSCILLATIONS*pi*y)*fac/2.0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  double fac;
  double pi = 3.14159265358979;

  switch(BdComp)
  {
  case 0: 
  case 2:
    value=0;
    break;
    case 1: 
      fac = exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*t/SINCOSEXP_ENDTIME);
      value = -cos(N_OSCILLATIONS*pi)*sin(N_OSCILLATIONS*pi*Param)*fac;
      break;
    case 3: 
      fac = exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*t/SINCOSEXP_ENDTIME);
      value = -sin(N_OSCILLATIONS*pi*(1-Param))*fac;
      break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  double fac;
  double pi = 3.14159265358979;
  value = 0;
  switch(BdComp)
  {
  case 1: 
  case 3:
    value=0;
    break;
  case 0: 
    fac = exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*t/SINCOSEXP_ENDTIME);
    value = sin(N_OSCILLATIONS*pi*Param)*fac;
    break;
  case 2: 
    fac = exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*t/SINCOSEXP_ENDTIME);
    value = sin(N_OSCILLATIONS*pi*(1-Param))*cos(N_OSCILLATIONS*pi)*fac;
    break;
  default: cout << "wrong boundary part number" << endl;
    break;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff, x, y; 
  static double a=1;
  static double nu=1/TDatabase::ParamDB->RE_NR;
  double gamma =  TDatabase::ParamDB->GAUSSIAN_GAMMA;
  double pi = 3.14159265358979;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double scale1, scale2, tmp, delta, c;
  double facu, facu_t, facp;
  double u1, u2, p, u1_t, u2_t, grad_u[4], u1_x, u1_y;
  double u2_x, u2_y, p_x, p_y, du_11, du_12, du_21, du_22;
  double u1_xx, u1_xy, u1_yy, u2_xx, u2_xy, u2_yy;
  double div_1, div_2, nu_T,nu_T_x,nu_T_y;

  // compute characteristic filter width
  // it is assumed that the filter width does not depend on h
  // such that tmp=1 is just a dummy input
  tmp = 1;
  delta =  CharacteristicFilterWidth(tmp);
  c = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT*delta*delta;

  scale1 = N_OSCILLATIONS*pi;
  scale2 = scale1*scale1;
  facu = exp(-2*scale2*t/SINCOSEXP_ENDTIME);
  facu_t = -2*scale2 * exp(-2*scale2*t/SINCOSEXP_ENDTIME)/ SINCOSEXP_ENDTIME;
  facp = exp(-4*scale2*t/SINCOSEXP_ENDTIME);

  //cout << c << " ";
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    // coordinates
    x = X[i];
    y = Y[i];
    // viscosity
    coeff[0] = nu;
    // build rhs
     
    u1  = -cos(scale1*x)*sin(scale1*y)*facu;
    u2 = sin(scale1*x)*cos(scale1*y)*facu;
    p = -(cos(2*scale1*x)+cos(2*scale1*y))*facp/4.0;

    u1_t =  -cos(scale1*x)*sin(scale1*y)*facu_t;
    u2_t = sin(scale1*x)*cos(scale1*y)*facu_t;
    
    grad_u[0] = u1_x = scale1*sin(scale1*x)*sin(scale1*y)*facu;
    grad_u[2] = u1_y = -scale1*cos(scale1*x)*cos(scale1*y)*facu;
    grad_u[1] = u2_x = scale1*cos(scale1*x)*cos(scale1*y)*facu;
    grad_u[3] = u2_y = -scale1*sin(scale1*x)*sin(scale1*y)*facu;

    p_x = scale1*sin(2*scale1*x)*facp/2.0;
    p_y = scale1*sin(2*scale1*y)*facp/2.0;

    du_11 = u1_x;
    du_12 = (u2_x + u1_y)/2.0;
    du_21 = du_12;
    du_22 = u2_y;

    u1_xx = scale2*cos(scale1*x)*sin(scale1*y)*facu;
    u1_xy = scale2*sin(scale1*x)*cos(scale1*y)*facu;
    u1_yy = scale2*cos(scale1*x)*sin(scale1*y)*facu;
    u2_xx = -scale2*sin(scale1*x)*cos(scale1*y)*facu;
    u2_xy = -scale2*cos(scale1*x)*sin(scale1*y)*facu;
    u2_yy = -scale2*sin(scale1*x)*cos(scale1*y)*facu;
    
// divergence of the Laplacian
    div_1 = u1_xx + u1_yy/2.0 + u2_xy/2.0;
    div_2 = u2_xx/2.0 + u1_xy/2.0 + u2_yy;
    
// turbulent viscosity (Smagorinsky term)    
    nu_T = TurbulentViscosity(delta, grad_u, grad_u, grad_u);
    if (nu_T == 0)
    {
      nu_T_x = 0;
      nu_T_y = 0;
    }
    else
    {
      nu_T_x = (2*u1_x*u1_xx+u1_y*u1_xy + u1_xy*u2_x + u1_y*u2_xx 
		+ u2_x * u2_xx + 2*u2_y*u2_xy);
      nu_T_x *= c*c/(2*nu_T);
      nu_T_y = (2*u1_x*u1_xy + u1_y*u1_yy + u1_yy*u2_x + u1_y * u2_xy
		+ u2_x * u2_xy + 2*u2_y*u2_yy);
      nu_T_y *= c*c/(2*nu_T);
    }
 
    coeff[1] = u1_t - (2*nu + nu_T) *div_1 -(du_11*nu_T_x + du_12*nu_T_y)
      + u1 * u1_x + u2 * u1_y + p_x; 
    coeff[2] = u2_t -  (2*nu + nu_T)*div_2 -(du_21*nu_T_x + du_22*nu_T_y)
      + u1 * u2_x + u2 * u2_y + p_y;

  }
}
