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

#define N_OSCILLATIONS    4
#define SINCOSEXP_ENDTIME    10000.0
#define __PAPER2__

void ExampleFile()
{
  OutPut("Example: SinCosExpT_paper2.h : N_OSCILLATIONS " << N_OSCILLATIONS ) 
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
  values[1] = 2*N_OSCILLATIONS*pi*sin(2*N_OSCILLATIONS*pi*x)*fac/4.0;
  values[2] = 2*N_OSCILLATIONS*pi*sin(2*N_OSCILLATIONS*pi*y)*fac/4.0;
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
  double pi = 3.14159265358979;
  double t = TDatabase::TimeDB->CURRENTTIME;
  double t1,t4,t6,t2,t8,t11,t14,t18, t7;
  double delta, tmp, fac,c, d11, d22, dnorm;
  double ddnorm_dx, ddnorm_dy, d11_dx, d11_dy, d22_dx, d22_dy;

  t1 = N_OSCILLATIONS*pi; 
  t4 = N_OSCILLATIONS*N_OSCILLATIONS; 
  t6 = pi*pi;
  t7 = t4*t6;

  fac = exp(-2*t7*t/SINCOSEXP_ENDTIME);

  // compute characteristic filter width
  // it is assumed that the filter width does not depend on h
  // such that tmp=1 is just a dummy input
  tmp = 1;
  delta =  CharacteristicFilterWidth(tmp);

  c = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT*delta*delta;
  //cout << c << " ";
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = nu;

    t2 = t1*x;
    t8 = t1*y;
    t11 = 1.0/SINCOSEXP_ENDTIME;
    t14 = exp(-2.0*t4*t6*t*t11);
    t18 = (-1.0+nu*SINCOSEXP_ENDTIME)*t11;

    coeff[1] = -2.0*cos(t2)*t4*t6*sin(t8)*t14*t18;
    coeff[2] = 2.0*cos(t8)*t4*t6*sin(t2)*t14*t18;

    d11 = t1*sin(t1*x)*sin(t1*y)*fac;
    d22 = -t1* sin(t1*x)*sin(t1*y)*fac;
    
    d11_dx = t7*cos(t1*x)*sin(t1*y)*fac;
    d11_dy = t7*sin(t1*x)*cos(t1*y)*fac;
    d22_dx = -t7*cos(t1*x)*sin(t1*y)*fac;
    d22_dy = -t7*sin(t1*x)*cos(t1*y)*fac;

    dnorm = sqrt(d11*d11+d22*d22);
    ddnorm_dx = (d11*d11_dx+d22*d22_dx)/dnorm;
    ddnorm_dy = (d11*d11_dy+d22*d22_dy)/dnorm;

    coeff[1] += (-ddnorm_dx*d11-dnorm*d11_dx)*c;
    coeff[2] += (-ddnorm_dy*d22-dnorm*d22_dy)*c;    
   
  }
}

// ========================================================================
// current deformation tensor 
// ========================================================================
void RhsAuxProblemPaper2(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  double t=TDatabase::TimeDB->CURRENTTIME, fac;
  double u1_x, u2_x, u1_y, u2_y;
  double pi = 3.14159265358979;
  int i;
  double *coeff, x, y; 

  fac = exp(-2*N_OSCILLATIONS*N_OSCILLATIONS*pi*pi*t/SINCOSEXP_ENDTIME);
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    u1_x = N_OSCILLATIONS*pi*sin(N_OSCILLATIONS*pi*x)*sin(N_OSCILLATIONS*pi*y)*fac;
    u1_y =  -N_OSCILLATIONS*pi*cos(N_OSCILLATIONS*pi*x)*cos(N_OSCILLATIONS*pi*y)*fac;
    u2_x = N_OSCILLATIONS*pi*cos(N_OSCILLATIONS*pi*x)*cos(N_OSCILLATIONS*pi*y)*fac;
    u2_y = -N_OSCILLATIONS*pi* sin(N_OSCILLATIONS*pi*x)*sin(N_OSCILLATIONS*pi*y)*fac;
    
// (u1_x*u1_x)+(u1_y*u1_y)
    coeff[0] = u1_x*u1_x + u1_y*u1_y;
    
// u1_x*u2_x + u1_y * u2_y
    coeff[1] = u1_x * u2_x + u1_y * u2_y;
    
// u2_x*u2_x + u2_y*u2_y
    coeff[2] =  u2_x*u2_x + u2_y*u2_y;
  }
}


