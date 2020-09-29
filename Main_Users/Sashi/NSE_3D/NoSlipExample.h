// Navier-Stokes problem, solution in ansatz space
// velocity pw quadratic, pressure linear
//
//
void ExampleFile()
{
  OutPut("Example: NoSlipExample.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  double sinx, dsinx, ddsinx, dddsinx, siny, dsiny, ddsiny, dddsiny;
  double polyz, dpolyz, ddpolyz, dddpolyz;

  sinx = sin(Pi*x)*sin(Pi*x);
  dsinx = 2*Pi*sin(Pi*x)*cos(Pi*x);
  ddsinx = 2*Pi*Pi*(cos(Pi*x)*cos(Pi*x)-sin(Pi*x)*sin(Pi*x));
  dddsinx = -8*Pi*Pi*Pi*cos(Pi*x)*sin(Pi*x);
  siny = sin(2*Pi*y)*sin(2*Pi*y);
  dsiny = 4*Pi*sin(2*Pi*y)*cos(2*Pi*y);
  ddsiny = 8*Pi*Pi*(cos(2*Pi*y)*cos(2*Pi*y) - sin(2*Pi*y)*sin(2*Pi*y));
  dddsiny = -64*Pi*Pi*Pi*sin(2*Pi*y)*cos(2*Pi*y);
  polyz = z*z*(1-z)*(1-z);
  dpolyz = 2*z*(1-z)*(1-2*z);
  ddpolyz = 2-12*z+12*z*z;
  dddpolyz = 24*z-12;

  values[0] = sinx*dsiny*polyz+sinx*siny*dpolyz;
  values[1] = dsinx*dsiny*polyz+dsinx*siny*dpolyz;
  values[2] = sinx*ddsiny*polyz+sinx*dsiny*dpolyz;
  values[3] = sinx*dsiny*dpolyz+sinx*siny*ddpolyz;
  values[4] = ddsinx*dsiny*polyz+ddsinx*siny*dpolyz+sinx*dddsiny*polyz+sinx*ddsiny*dpolyz
     + sinx*dsiny*ddpolyz+sinx*siny*dddpolyz;
}

void ExactU2(double x, double y,  double z, double *values)
{
  double sinx, dsinx, ddsinx, dddsinx, siny, dsiny, ddsiny, dddsiny;
  double polyz, dpolyz, ddpolyz, dddpolyz;

  sinx = sin(Pi*x)*sin(Pi*x);
  dsinx = 2*Pi*sin(Pi*x)*cos(Pi*x);
  ddsinx = 2*Pi*Pi*(cos(Pi*x)*cos(Pi*x)-sin(Pi*x)*sin(Pi*x));
  dddsinx = -8*Pi*Pi*Pi*cos(Pi*x)*sin(Pi*x);
  siny = sin(2*Pi*y)*sin(2*Pi*y);
  dsiny = 4*Pi*sin(2*Pi*y)*cos(2*Pi*y);
  ddsiny = 8*Pi*Pi*(cos(2*Pi*y)*cos(2*Pi*y) - sin(2*Pi*y)*sin(2*Pi*y));
  dddsiny = -64*Pi*Pi*Pi*sin(2*Pi*y)*cos(2*Pi*y);
  polyz = z*z*(1-z)*(1-z);
  dpolyz = 2*z*(1-z)*(1-2*z);
  ddpolyz = 2-12*z+12*z*z;
  dddpolyz = 24*z-12;

  values[0] = -dsinx*siny*polyz+sinx*siny*dpolyz;
  values[1] = -ddsinx*siny*polyz+dsinx*siny*dpolyz;
  values[2] = -dsinx*dsiny*polyz+sinx*dsiny*dpolyz;
  values[3] = -dsinx*siny*dpolyz+sinx*siny*ddpolyz;
  values[4] = -dddsinx*siny*polyz+ddsinx*siny*dpolyz-dsinx*ddsiny*polyz+sinx*ddsiny*dpolyz
   -dsinx*siny*ddpolyz+sinx*siny*dddpolyz;  

}

void ExactU3(double x, double y,  double z, double *values)
{
  double sinx, dsinx, ddsinx, dddsinx, siny, dsiny, ddsiny, dddsiny;
  double polyz, dpolyz, ddpolyz, dddpolyz;

  sinx = sin(Pi*x)*sin(Pi*x);
  dsinx = 2*Pi*sin(Pi*x)*cos(Pi*x);
  ddsinx = 2*Pi*Pi*(cos(Pi*x)*cos(Pi*x)-sin(Pi*x)*sin(Pi*x));
  dddsinx = -8*Pi*Pi*Pi*cos(Pi*x)*sin(Pi*x);
  siny = sin(2*Pi*y)*sin(2*Pi*y);
  dsiny = 4*Pi*sin(2*Pi*y)*cos(2*Pi*y);
  ddsiny = 8*Pi*Pi*(cos(2*Pi*y)*cos(2*Pi*y) - sin(2*Pi*y)*sin(2*Pi*y));
  dddsiny = -64*Pi*Pi*Pi*sin(2*Pi*y)*cos(2*Pi*y);
  polyz = z*z*(1-z)*(1-z);
  dpolyz = 2*z*(1-z)*(1-2*z);
  ddpolyz = 2-12*z+12*z*z;
  dddpolyz = 24*z-12;

  values[0] = -dsinx*siny*polyz-sinx*dsiny*polyz;
  values[1] = -ddsinx*siny*polyz-dsinx*dsiny*polyz;
  values[2] = -dsinx*dsiny*polyz-sinx*ddsiny*polyz;
  values[3] = -dsinx*siny*dpolyz-sinx*dsiny*dpolyz;
  values[4] = -dddsinx*siny*polyz-ddsinx*dsiny*polyz-dsinx*ddsiny*polyz-sinx*dddsiny*polyz
     -dsinx*siny*ddpolyz-sinx*dsiny*ddpolyz;

}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 3*x-sin(y+4*z)-1.5-sin(1.0)*cos(1.0)
    -4*cos(1.0)*cos(1.0)*cos(1.0)*cos(1.0)*sin(1.0)
    +sin(1.0)*cos(1.0)*cos(1.0)-2*sin(1.0)*sin(1.0)*sin(1.0)
    +2*sin(1.0)*cos(1.0)*cos(1.0)*cos(1.0)+2*sin(1.0);
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
  value = 0;
}

// value of boundary condition
void U2BoundValue(double x, double y, double z, double &value)
{
  value =  0;
}

// value of boundary condition
void U3BoundValue(double x, double y, double z, double &value)
{
  value = 0;
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
  double sinx, dsinx, ddsinx, dddsinx, siny, dsiny, ddsiny, dddsiny;
  double polyz, dpolyz, ddpolyz, dddpolyz;
  double u1, u1x, u1y, u1z, lapl_u1;
  double u2, u2x, u2y, u2z, lapl_u2;
  double u3, u3x, u3y, u3z, lapl_u3;
  double px, py, pz;

  if (TDatabase::ParamDB->STOKES_PROBLEM)
  {
    for(i=0;i<n_points;i++)
    {
      coeff = coeffs[i];
      x = X[i];
      y = Y[i];
      z = Z[i];

      sinx = sin(Pi*x)*sin(Pi*x);
      dsinx = 2*Pi*sin(Pi*x)*cos(Pi*x);
      ddsinx = 2*Pi*Pi*(cos(Pi*x)*cos(Pi*x)-sin(Pi*x)*sin(Pi*x));
      dddsinx = -8*Pi*Pi*Pi*cos(Pi*x)*sin(Pi*x);
      siny = sin(2*Pi*y)*sin(2*Pi*y);
      dsiny = 4*Pi*sin(2*Pi*y)*cos(2*Pi*y);
      ddsiny = 8*Pi*Pi*(cos(2*Pi*y)*cos(2*Pi*y) - sin(2*Pi*y)*sin(2*Pi*y));
      dddsiny = -64*Pi*Pi*Pi*sin(2*Pi*y)*cos(2*Pi*y);
      polyz = z*z*(1-z)*(1-z);
      dpolyz = 2*z*(1-z)*(1-2*z);
      ddpolyz = 2-12*z+12*z*z;
      dddpolyz = 24*z-12;

      u1 = sinx*dsiny*polyz+sinx*siny*dpolyz;
      u1x = dsinx*dsiny*polyz+dsinx*siny*dpolyz;
      u1y = sinx*ddsiny*polyz+sinx*dsiny*dpolyz;
      u1z = sinx*dsiny*dpolyz+sinx*siny*ddpolyz;
      lapl_u1 = ddsinx*dsiny*polyz+ddsinx*siny*dpolyz+sinx*dddsiny*polyz+sinx*ddsiny*dpolyz
         + sinx*dsiny*ddpolyz+sinx*siny*dddpolyz;
     
      u2 = -dsinx*siny*polyz+sinx*siny*dpolyz;
      u2x = -ddsinx*siny*polyz+dsinx*siny*dpolyz;
      u2y = -dsinx*dsiny*polyz+sinx*dsiny*dpolyz;
      u2z = -dsinx*siny*dpolyz+sinx*siny*ddpolyz;
      lapl_u2 = -dddsinx*siny*polyz+ddsinx*siny*dpolyz-dsinx*ddsiny*polyz+sinx*ddsiny*dpolyz
         -dsinx*siny*ddpolyz+sinx*siny*dddpolyz;  

      u3 = -dsinx*siny*polyz-sinx*dsiny*polyz;
      u3x = -ddsinx*siny*polyz-dsinx*dsiny*polyz;
      u3y = -dsinx*dsiny*polyz-sinx*ddsiny*polyz;
      u3z = -dsinx*siny*dpolyz-sinx*dsiny*dpolyz;
      lapl_u3 = -dddsinx*siny*polyz-ddsinx*dsiny*polyz-dsinx*ddsiny*polyz-sinx*dddsiny*polyz
         -dsinx*siny*ddpolyz-sinx*dsiny*ddpolyz;

      px = 3;
      py = -cos(y+4*z);
      pz = -4*cos(y+4*z);

      coeff[0] = eps;

      coeff[1] = -eps*lapl_u1+px; // f1
      coeff[2] = -eps*lapl_u2+py; // f2
      coeff[3] = -eps*lapl_u3+pz; // f3
 
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

      sinx = sin(Pi*x)*sin(Pi*x);
      dsinx = 2*Pi*sin(Pi*x)*cos(Pi*x);
      ddsinx = 2*Pi*Pi*(cos(Pi*x)*cos(Pi*x)-sin(Pi*x)*sin(Pi*x));
      dddsinx = -8*Pi*Pi*Pi*cos(Pi*x)*sin(Pi*x);
      siny = sin(2*Pi*y)*sin(2*Pi*y);
      dsiny = 4*Pi*sin(2*Pi*y)*cos(2*Pi*y);
      ddsiny = 8*Pi*Pi*(cos(2*Pi*y)*cos(2*Pi*y) - sin(2*Pi*y)*sin(2*Pi*y));
      dddsiny = -64*Pi*Pi*Pi*sin(2*Pi*y)*cos(2*Pi*y);
      polyz = z*z*(1-z)*(1-z);
      dpolyz = 2*z*(1-z)*(1-2*z);
      ddpolyz = 2-12*z+12*z*z;
      dddpolyz = 24*z-12;

      u1 = sinx*dsiny*polyz+sinx*siny*dpolyz;
      u1x = dsinx*dsiny*polyz+dsinx*siny*dpolyz;
      u1y = sinx*ddsiny*polyz+sinx*dsiny*dpolyz;
      u1z = sinx*dsiny*dpolyz+sinx*siny*ddpolyz;
      lapl_u1 = ddsinx*dsiny*polyz+ddsinx*siny*dpolyz+sinx*dddsiny*polyz+sinx*ddsiny*dpolyz
         + sinx*dsiny*ddpolyz+sinx*siny*dddpolyz;
     
      u2 = -dsinx*siny*polyz+sinx*siny*dpolyz;
      u2x = -ddsinx*siny*polyz+dsinx*siny*dpolyz;
      u2y = -dsinx*dsiny*polyz+sinx*dsiny*dpolyz;
      u2z = -dsinx*siny*dpolyz+sinx*siny*ddpolyz;
      lapl_u2 = -dddsinx*siny*polyz+ddsinx*siny*dpolyz-dsinx*ddsiny*polyz+sinx*ddsiny*dpolyz
         -dsinx*siny*ddpolyz+sinx*siny*dddpolyz;  

      u3 = -dsinx*siny*polyz-sinx*dsiny*polyz;
      u3x = -ddsinx*siny*polyz-dsinx*dsiny*polyz;
      u3y = -dsinx*dsiny*polyz-sinx*ddsiny*polyz;
      u3z = -dsinx*siny*dpolyz-sinx*dsiny*dpolyz;
      lapl_u3 = -dddsinx*siny*polyz-ddsinx*dsiny*polyz-dsinx*ddsiny*polyz-sinx*dddsiny*polyz
         -dsinx*siny*ddpolyz-sinx*dsiny*ddpolyz;

      px = 3;
      py = -cos(y+4*z);
      pz = -4*cos(y+4*z);

      coeff[0] = eps;

      coeff[1] = -eps*lapl_u1+u1*u1x*u2*u1y+u3*u1z+px; // f1
      coeff[2] = -eps*lapl_u2+u1*u2x*u2*u2y+u3*u2z+py; // f2
      coeff[3] = -eps*lapl_u3+u1*u3x*u2*u3y+u3*u3z+pz; // f3
 
    }
  }
    
}
void LinCoeffsOld(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y, z; 
  double sinx, dsinx, ddsinx, dddsinx, siny, dsiny, ddsiny, dddsiny;
  double polyz, dpolyz, ddpolyz, dddpolyz;
  double u1, u1x, u1y, u1z, lapl_u1;
  double u2, u2x, u2y, u2z, lapl_u2;
  double u3, u3x, u3y, u3z, lapl_u3;
  double px, py, pz;
  double t1,t2,t3,t4,t5,t8,t9,t11,t12,t13,t15,t19,t20,t21,t22,t23,t27,t28,t33,t34;
  double t35,t42,t43,t46,t47,t52,t55,t57,t68,t71,t72,t74,t83,t86,t89,t96,t99,t102;
  double t103, t7, t14, t25, t31, t39, t40, t51, t53,t26,t36,t56,t58,t69,t73, t76,t77;
  double t79,t90,t92,t95,t122,t29,t41,t44,t93,t82,t94,t114;

  for(i=0;i<n_points;i++)
  {
     coeff = coeffs[i];
     x = X[i];
     y = Y[i];
     z = Z[i];
     
     t1 = Pi*x;
     t2 = cos(t1);
     t3 = t2*t2;
     t4 = Pi*Pi;
     t5 = t4*Pi;
     t8 = 2.0*Pi*y;
     t9 = sin(t8);
     t11 = z*z;
     t12 = 1.0-z;
     t13 = t12*t12;
     t15 = cos(t8);
     t19 = sin(t1);
     t20 = t19*t19;
     t21 = t20*t9;
     t22 = t21*t11;
     t23 = t13*t15;
     t27 = t3*t4;
     t28 = t9*t9;
     t33 = t20*t28;
     t34 = z*t13;
     t35 = t34*t4;
     t42 = t11*t12;
     t43 = t42*t4;
     t46 = t15*t15;
     t47 = t20*t46;
     t52 = t23*Pi;
     t55 = t21*z;
     t57 = t12*t15*Pi;
     t68 = 8.0*t3*t5*t9*t11*t13*t15-72.0*t22*t23*t5+4.0*t27*t28*z*t13-20.0*t33
        *t35-4.0*t27*t28*t11*t12+20.0*t33*t43+16.0*t47*t35-16.0*t47*t43+8.0*t21*t52
        -32.0*t55*t57+8.0*t21*t11*t15*Pi-12.0*t33*t12+12.0*t33*z;
     t71 = 4.0*t22*t52;
     t72 = t33*t34;
     t74 = t33*t42;
     t83 = t19*t28;
     t86 = t13*t2*Pi;
     t89 = t83*t11;
     t96 = t89*t86;
     t99 = t4*t11*t13;
     t102 = t55*t52;
     t103 = t22*t57;
     coeff[0] = eps;
     coeff[1] = eps*t68+(t71+2.0*t72-2.0*t74)*(8.0*t19*t9*t11*t23*t4*t2+4.0*
                                               t83*z*t86-4.0*t89*t12*t2*Pi)+16.0*(-t96+t72-t74)*(t47*t99-t33
                                                                                                 *t99+t102-t103)+(-2.0*t96-t71)*
        (8.0*t102-8.0*t103+2.0*t33*t13-8.0*t33*z*t12+2.0*t33*t11)+3.0;
     t1 = Pi*x;
     t2 = cos(t1);
     t3 = Pi*Pi;
     t4 = t3*Pi;
     t7 = 2.0*Pi*y;
     t8 = sin(t7);
     t9 = t8*t8;
     t11 = z*z;
     t12 = 1.0-z;
     t13 = t12*t12;
     t14 = t11*t13;
     t15 = sin(t1);
     t19 = t2*t2;
     t20 = t19*t3;
     t25 = t15*t15;
     t26 = t25*t9;
     t27 = z*t13;
     t28 = t27*t3;
     t31 = t9*t11;
     t35 = t11*t12;
     t36 = t35*t3;
     t39 = cos(t7);
     t40 = t39*t39;
     t46 = t25*t40;
     t51 = t15*t9;
     t53 = t13*t2*Pi;
     t56 = t51*z;
     t58 = t12*t2*Pi;
     t69 = 24.0*t2*t4*t9*t14*t15+4.0*t20*t9*z*t13-20.0*t26*t28-4.0*t20*t31*t12
        +20.0*t26*t36-16.0*t15*t40*t4*t14*t2+16.0*t46*t28-16.0*t46*t36-4.0*t51*t53+16.0
        *t56*t58-4.0*t51*t11*t2*Pi-12.0*t26*t12+12.0*t26*z;
     t71 = t25*t8;
     t72 = t71*t11;
     t73 = t13*t39;
     t74 = t73*Pi;
     t76 = 4.0*t72*t74;
     t77 = t26*t27;
     t79 = t26*t35;
     t89 = 4.0*t56*t53;
     t90 = t51*t11;
     t92 = 4.0*t90*t58;
     t95 = t90*t53;
     t122 = cos(y+4.0*z);
     coeff[2] = eps*t69+(t76+2.0*t77-2.0*t79)*(-2.0*t20*t31*t13+2.0*t26*t14*t3
                                               +t89-t92)+16.0*(-t95+t77-t79)*(-t15*t8*t11*t73*t3*t2+t71*z*t74-t72*t12*t39*
                                                                              Pi)+(-2.0*t95-t76)*(-t89+t92+2.0*t26*t13-8.0*t26*z*t12+2.0*
                                                                                                  t26*t11)-t122;
     t1 = Pi*x;
     t2 = cos(t1);
     t3 = Pi*Pi;
     t4 = t3*Pi;
     t7 = 2.0*Pi*y;
     t8 = sin(t7);
     t9 = t8*t8;
     t11 = z*z;
     t12 = 1.0-z;
     t13 = t12*t12;
     t14 = t11*t13;
     t15 = sin(t1);
     t19 = t2*t2;
     t22 = cos(t7);
     t26 = t15*t15;
     t27 = t26*t8;
     t28 = t27*t11;
     t29 = t13*t22;
     t33 = t22*t22;
     t39 = t15*t9;
     t41 = t13*t2*Pi;
     t44 = t39*z;
     t46 = t12*t2*Pi;
     t53 = t29*Pi;
     t56 = t27*z;
     t58 = t12*t22*Pi;
     t68 = 4.0*t28*t53;
     t69 = t26*t9;
     t71 = t69*z*t13;
     t74 = t69*t11*t12;
     t82 = t14*t3;
     t83 = t69*t82;
     t89 = t15*t8*t11*t29*t3*t2;
     t93 = t39*t11;
     t94 = t93*t41;
     t114 = cos(y+4.0*z);
     coeff[3] = eps*(24.0*t2*t4*t9*t14*t15-8.0*t19*t4*t8*t14*t22+72.0*t28*t29*t4-16.0*t15*t33*t4*t14*t2-4.0*t39*t41+16.0*t44*t46-4.0*t39*t11*t2*
                     Pi-8.0*t27*t53+32.0*t56*t58-8.0*t27*t11*t22*
                     Pi)+(t68+2.0*t71-2.0*t74)*(-2.0*t19*t3*t9*t11*t13+2.0*t83-8.0
                                                *t89)+16.0*(-t94+t71-t74)*(-t89-t26*t33*t82+t83)+(-2.0*t94-t68)*(-4.0*t44*t41+
                                                                                                                 4.0*t93*t46-8.0*t56*t53+8.0*t28*t58)-4.0*t114;
  }
}
