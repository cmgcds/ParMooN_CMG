// ==========================================================================
// coupled linear system
// ==========================================================================
#include <TimeConvDiff2D_2.h>

void ExampleFile()
{
  OutPut("Exzmple: coupled1.h" << endl);
}

// exact solution
void Exact_U(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  double t1 = exp(-50*t);
  
  values[0] = t1*(x+2*y);
  values[1] = t1;
  values[2] = 2*t1;
  values[3] = 0;
}

void Exact_V(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  
  values[0] = t*(2*x+y);
  values[1] = 2*t;
  values[2] = t;
  values[3] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double Param, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue_U(int BdComp, double Param, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  double t3=exp(-50*t);
  
  switch(BdComp)
  {
    case 0: value=Param*t3;
            break;
    case 1: value=(1+2*Param)*t3;
            break;
    case 2: value= (3-Param)*t3;
            break;
    case 3: value=(2-2*Param)*t3;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  } // endswitch
}

void BoundValue_V(int BdComp, double Param, double &value)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  
  switch(BdComp)
  {
    case 0: value=2*Param*t;
            break;
    case 1: value=(2+Param)*t;
            break;
    case 2: value= (3-2*Param)*t;
            break;
    case 3: value=(1-Param)*t;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  } // endswitch
}

void BilinearCoeffs_U(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{ static double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;
  double x, y;
  double t=TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps; // Laplace u
    coeff[1] = 0;   // u_x
    coeff[2] = 0;   // u_y
    coeff[3] = 50;  // u
    coeff[4] = eps; // Laplace v
    coeff[5] = 1;   // v_x
    coeff[6] = 0;   // v_y
    coeff[7] = 0;   // v
    coeff[8] = t;   // f
  }
}

void BilinearCoeffs_V(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{ static double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;
  double x, y;
  double t=TDatabase::TimeDB->CURRENTTIME;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    x = X[i];
    y = Y[i];

    coeff[0] = eps; // Laplace u
    coeff[1] = 0;   // u_x
    coeff[2] = 1;   // u_y
    coeff[3] = 0;   // u
    coeff[4] = eps; // Laplace v
    coeff[5] = 0;   // v_x
    coeff[6] = 0;   // v_y
    coeff[7] = 1;   // v
    coeff[8] = (2*x+y)*(t+1)+2*exp(-50*t);   // f
  }
}


// exact solution
void Initial_U(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  double t1 = exp(-50*t);
  
  values[0] = t1*(x+2*y);
  values[1] = t1;
  values[2] = 2*t1;
  values[3] = 0;
}

void Initial_V(double x, double y, double *values)
{
  double t=TDatabase::TimeDB->CURRENTTIME;
  
  values[0] = t*(2*x+y);
  values[1] = 2*t;
  values[2] = t;
  values[3] = 0;
}
