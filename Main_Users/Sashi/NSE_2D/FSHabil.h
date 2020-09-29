// Navier-Stokes problem from Friedhelms Habil
// 

void ExampleFile()
{
  OutPut("Example: FSHabil.h" << endl) ;
  TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY = OSEEN_PROBLEM;
  OutPut("Example: FSHabil.h with INTERNAL_PROBLEM_IDENTITY " << 
	 TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY <<  endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t9, t10, t13, t15;

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
}

void ExactU2(double x, double y, double *values)
{
  double t1, t2, t3, t4, t5, t6, t7, t9, t10, t20, t23;

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
}

void ExactP(double x, double y, double *values)
{
  values[0] = x*x*x+y*y*y-0.5;
  values[1] = 3*x*x;
  values[2] = 3*y*y;
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
  double s2, s3, s5, s6, s10, s11, s13, s16, s17, s23, s24, s26, s32;
  double t1, t4, t7, t15, t20;
  int i;
  double *coeff, x, y;

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
                     +16.0*t8*t17+4.0*t11*t14*t5+12.0*t21*t5+12.0*t21*y)
               +3.0*t11;

    s2 = -1.0+x;
    s3 = y*y;
    s5 = -1.0+y;
    s6 = s5*s5;
    s10 = s2*s2;
    s11 = x*s10;
    s13 = y*s5;
    s16 = x*x;
    s17 = s16*s2;
    s23 = s16*s10;
    s24 = y*s6;
    s26 = s3*s5;
    s32 = s3*s6;
    coeff[2] = -nu*(-12.0*s2*s3*s6-12.0*x*s3*s6-4.0*s11*s6-16.0*s11*s13
                     -4.0*s11*s3-4.0*s17*s6-16.0*s17*s13-4.0*s17*s3)
               +3.0*s3;

    if( !(TDatabase::ParamDB->STOKES_PROBLEM) )
    {
      coeff[1] +=
           +(2.0*t21*t9+2.0*t21*t17)*(4.0*t29*t9+4.0*t31*t9+4.0*t29*t17
                                      +4.0*t31*t17)
           +(-2.0*t29*t37-2.0*t31*t37)*(2.0*t21*t6+8.0*t21*y*t5
                                        +2.0*t21*t14);
                                        
      coeff[2] +=
           +(2.0*s23*s24+2.0*s23*s26)*(-2.0*s10* s3*s6-8.0*x*s2*s32
                                       -2.0*s16*s3*s6)
           +(-2.0*s11*s32-2.0*s17*s32)*(-4.0*s11*s24-4.0*s11*s26
                                        -4.0*s17*s24-4.0*s17*s26);

      // for Oseen
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
      coeff[3] = 2.0*t4*t7-2.0*t4*t10;
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
      coeff[4] = -2.0*t3*t7+2.0*t10*t7;
    }
  }
}

