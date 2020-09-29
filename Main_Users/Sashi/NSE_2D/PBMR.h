// Navier-Stokes problem, PBMR
// 
// u(x,y) velosity
// p(x,y) pressure

void ExampleFile()
{
  OutPut("Example: PBMR.h" << endl) ;
}

//#include <NavierStokes.h>
#include <Database.h>

// start value for the nonlinear iteration
void InitialU1(double x, double y, double *values)
{
  static double v=TDatabase::ParamDB->REACTOR_P1;
  static double d=TDatabase::ParamDB->REACTOR_P2;
  
  if(TDatabase::ParamDB->REACTOR_P3==0)
  { 
    values[0]=v;
    values[2]=0;
    values[3]=0;
  }
  else 
  {
    values[0]=6*v*d*y*(d-d*y)/(d*d);
    values[2]=6*v-3*v*y;
    values[3] = -3*v;
  } 
  values[1] = 0;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// =======================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=0;
  if(i==1)
    cond = NEUMANN;  else 
    cond = DIRICHLET;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  static double v=TDatabase::ParamDB->REACTOR_P2;
  static double v_y=TDatabase::ParamDB->REACTOR_P5;
  static double d=TDatabase::ParamDB->REACTOR_P3;

  switch(BdComp)
  {
    case 0: if(Param<0.125 || Param>0.875)
               value=0;
	     else
	       value=v_y;
            break;
    case 1: value=0;
            break;
    case 2: if(Param<0.125 || Param>0.875)
               value=0;
	     else
	       value=-v_y;
            break;
     case 3:  value = 0;
             break;
	      default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  static double v=TDatabase::ParamDB->REACTOR_P2;

  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3:  if(Param<0.0001 || Param>0.9999) 
              value = 0;
	     else
	       value=v;
              break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  static double f_friction1 = TDatabase::ParamDB->REACTOR_P0;
  static double f_friction2 = TDatabase::ParamDB->REACTOR_P1;
  int i;
  double *coeff;
  double *param;
  double v;


  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];


    coeff[0] = eps;
    coeff[1] = 0;
    //    coeff[1] = f_friction*16*y[i]*y[i]*(1-y[i])*(1-y[i]); // f1
    coeff[2] = 0; // f2
    coeff[3] = f_friction1;
    coeff[4] = f_friction2;
    // cout << "Friction parameter:  " << coeff[3];
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================

void NonLinCoeffs(int n_points, double *x, double *y,
                  double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  static double f_friction = TDatabase::ParamDB->REACTOR_P0;
  int i;
  double *coeff, *param;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = eps;
    coeff[1] = param[0];
    coeff[2] = param[1];
    coeff[3] = f_friction;
  }
}
