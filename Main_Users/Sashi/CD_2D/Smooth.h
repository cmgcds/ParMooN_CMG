// ======================================================================
// smooth solution problem
// ======================================================================
#include <ConvDiff2D.h>
#define __SMOOTH__

void ExampleFile()
{
  OutPut("Example: Smooth.h, convection(P1,P2) (" << TDatabase::ParamDB->P1
	 <<","<<TDatabase::ParamDB->P2<<"), reaction(P3) " << TDatabase::ParamDB->P3
	 << endl) ;
  TDatabase::ParamDB->INTERNAL_QUAD_RULE = 99;
}

// exact solution
void Exact(double x, double y, double *values)
{
  double t1,t2,t3,t4,t5,t6,t7,t9,t10,t12,t13;

  t1 = 1.0-x;
  t2 = t1*t1;
  t3 = x*x;
  t10 = 100.0*t2*t3*y*(1.0-2.0*y)*(1.0-y);
  values[0]=t10;

  t1 = 1.0-x;
  t2 = x*x;
  t7 = y*(1.0-2.0*y)*(1.0-y);
  t9 = t1*t1;
  t12 = -200.0*t1*t2*t7+200.0*t9*x*t7;
  values[1]=t12;

  t2 = pow(1.0-x,2.0);
  t3 = x*x;
  t4 = t2*t3;
  t5 = 1.0-2.0*y;
  t6 = 1.0-y;
  t13 = 100.0*t4*t5*t6-200.0*y*t6*t4-100.0*t4*y*t5;
  values[2]=t13;

  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;
  double t1,t3,t4,t5,t7,t9,t10,t12,t15,t30;
  double x, y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps;
    coeff[1] = TDatabase::ParamDB->P1;
    coeff[2] = TDatabase::ParamDB->P2;
    coeff[3] = TDatabase::ParamDB->P3;

    x = X[i];
    y = Y[i];
    t1 = x*x;
    t3 = 1.0-2.0*y;
    t4 = 1.0-y;
    t5 = t3*t4;
    t7 = 1.0-x;
    t9 = y*t3;
    t10 = t9*t4; 
    t12 = t7*t7; 
    t15 = t12*t1; 
    t30 = -eps*(200.0*t1*y*t5-800.0*t7*x*t10+200.0*t12*y*t5
                -400.0*t15*t4-200.0*t15*t3+400.0*t15*y)
	+coeff[1]*(-200.0*t7*t1*t10+200.0*t12*x*t10)
	+coeff[2]*(100.0*t15*t5-200.0*t15*y*t4-100.0*t15*t9)
	+coeff[3]*(100.0*t15*t10);

    coeff[4] = t30;
  }
}

void ErrorInGivenPoints(TFEFunction2D *ufct)
{
    double x,y,val_u[4],val_uh[4], step;
    
    step = 1./64.;
    if (strcmp(TDatabase::ParamDB->GEOFILE,"TwoTriangles")==0)
	step /= 2.0;
    x = 0.5 + step;
    if (strcmp(TDatabase::ParamDB->GEOFILE,"TwoTriangles")==0)
	step = -step;
    y = 0.25 + step;
    ufct->FindGradient(x,y,val_uh);
    Exact(x,y,val_u);
    OutPut("error in ("<<x<<","<<y<<"): " << fabs(val_uh[0]-val_u[0]) <<endl);
    y = 0.5   + step;
    ufct->FindGradient(x,y,val_uh);
    Exact(x,y,val_u);
    OutPut("error in ("<<x<<","<<y<<"): " << fabs(val_uh[0]-val_u[0]) <<endl);
    y = 0.75  + step;
    ufct->FindGradient(x,y,val_uh);
    Exact(x,y,val_u);
    OutPut("error in ("<<x<<","<<y<<"): " << fabs(val_uh[0]-val_u[0]) <<endl);
}
