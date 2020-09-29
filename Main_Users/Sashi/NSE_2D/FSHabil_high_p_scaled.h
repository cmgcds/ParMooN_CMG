// Navier-Stokes problem from Friedhelms Habil
// 
// ========================================================================
// example file
// ========================================================================

void ExampleFile()
{
  OutPut("Example: FSHabil_high_p_scaled.h,") 
  OutPut(" solution scale by "<< TDatabase::ParamDB->P7 << endl);
}
// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t9, t10, t13, t15;
  double scale = TDatabase::ParamDB->P7;

  t1 = x*x;
  t2 = 1.0-x;
  t3 = t2*t2;
  t4 = t1*t3;
  t5 = 1.0-y;
  t6 = t5*t5;
  t7 = y*t6;
  t9 = y*y;
  t10 = t9*t5;
  t13 = x*t3;
  t15 = t1*t2;
  values[0] = 2.0*t4*t7-2.0*t4*t10;
  values[1] = 4.0*t13*t7-4.0*t15*t7-4.0*t13*t10+4.0*t15*t10;
  values[2] = 2.0*t4*t6-8.0*t4*y*t5+2.0*t4*t9;
  values[3] = 0;
  values[0] *= scale;
  values[1] *= scale;
  values[2] *= scale;
  values[3] *= scale;
}

void ExactU2(double x, double y, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t9, t10, t20, t23;
  double scale = TDatabase::ParamDB->P7;

  t1 = 1.0-x;
  t2 = t1*t1;
  t3 = x*t2;
  t4 = y*y;
  t5 = 1.0-y;
  t6 = t5*t5;
  t7 = t4*t6;
  t9 = x*x;
  t10 = t9*t1;
  t20 = y*t6;
  t23 = t4*t5;
  values[0] = -2.0*t3*t7+2.0*t10*t7;
  values[1] = -2.0*t2*t4*t6+8.0*x*t1*t7-2.0*t4*t9*t6;
  values[2] = -4.0*t3*t20+4.0*t10*t20+4.0*t3*t23-4.0*t10*t23;
  values[3] = 0;
  values[0] *= scale;
  values[1] *= scale;
  values[2] *= scale;
  values[3] *= scale;
}

void ExactP(double x, double y, double *values)
{
  double scale = TDatabase::ParamDB->P7;

  values[0] = x*x*x*x*x+y*y*y*y*y-1.0/3.0;
  values[1] = 5*x*x*x*x;
  values[2] = 5*y*y*y*y;
  values[3] = 0;
  values[0] *= scale;
  values[1] *= scale;
  values[2] *= scale;
  values[3] *= scale;
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
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value=0;
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
  static double nu=1/TDatabase::ParamDB->RE_NR;
  double t2, t3, t5, t6, t8, t9, t11, t10, t13, t14, t16;
  double t17, t21, t23, t24, t26, t29, t31, t32, t37;
  int i;
  double *coeff, x, y;
  double scale = TDatabase::ParamDB->P7;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    
    x = X[i];
    y = Y[i];

    t2 = -1.0+x;
    t3 = t2*t2;
    t5 = -1.0+y;
    t6 = t5*t5;
    t8 = x*t2;
    t9 = y*t6;
    t11 = x*x;
    t14 = y*y;
    t17 = t14*t5;
    t21 = t11*t3;
    t29 = x*t3;
    t31 = t11*t2;
    t37 = t14*t6;
    coeff[0] = nu;
    coeff[1] = -nu*(4.0*t3*y*t6+16.0*t8*t9+4.0*t11*y*t6+4.0*t3*t14*t5
                     +16.0*t8*t17+4.0*t11*t14*t5+12.0*t21*t5+12.0*t21*y)*scale
               +(2.0*t21*t9+2.0*t21*t17)*(4.0*t29*t9+4.0*t31*t9+4.0*t29*t17
                                          +4.0*t31*t17)*scale*scale
               +(-2.0*t29*t37-2.0*t31*t37)*(2.0*t21*t6+8.0*t21*y*t5
                                            +2.0*t21*t14)*scale*scale
               +5*x*x*x*x*scale;

    t2 = -1.0+x;
    t3 = y*y;
    t5 = -1.0+y;
    t6 = t5*t5;
    t10 = t2*t2;
    t11 = x*t10;
    t13 = y*t5;
    t16 = x*x;
    t17 = t16*t2;
    t23 = t16*t10;
    t24 = y*t6;
    t26 = t3*t5;
    t32 = t3*t6;
    coeff[2] = -nu*(-12.0*t2*t3*t6-12.0*x*t3*t6-4.0*t11*t6-16.0*t11*t13
                     -4.0*t11*t3-4.0*t17*t6-16.0*t17*t13-4.0*t17*t3)*scale
               +(2.0*t23*t24+2.0*t23*t26)*(-2.0*t10* t3*t6-8.0*x*t2*t32
                                           -2.0*t16*t3*t6)*scale*scale
               +(-2.0*t11*t32-2.0*t17*t32)*(-4.0*t11*t24-4.0*t11*t26
                                            -4.0*t17*t24-4.0*t17*t26)*scale*scale
               +5.0*y*y*y*y*scale;
  }
}

