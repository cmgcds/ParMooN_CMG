// Navier-Stokes problem, rotating flow 
// for rotating system
// 
// u(x,y) = 2*Pi*(y,-x)
// p(x,y) = 2*Pi*Pi*(x^2+y^2) - Pi*Pi 
//
// one rotation per second, omega = -2*Pi  

void ExampleFile()
{
  OutPut("Example: RotatingFlow.h" << endl) ;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
    values[0] = 2*Pi*y;
    //values[0] = 0;
}

void InitialU2(double x, double y, double *values)
{

    values[0] = -2*Pi*x;
    //values[0] = 0;
}

void InitialP(double x, double y, double *values)
{
    values[0] =  2*Pi*Pi*(x*x+y*y)-Pi*Pi;
    //values[0] = 2*Pi*Pi*(x*x+y*y)-Pi*Pi;
    //values[0] = 0;
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
    values[0] = 2*Pi*y;
    values[1] = 0;
    values[2] = 2*Pi;
    values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
    values[0] = -2*Pi*x;
    values[1] = -2*Pi;
    values[2] = 0;
    values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
    values[0] = 2*Pi*Pi*(x*x+y*y)-Pi*Pi;
    values[1] = 4*Pi*Pi*x;
    values[2] = 4*Pi*Pi*y;
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
    // for Circle.PRM, parametrization of circle  t \in [0,2pi]
    // radius r = 1;
    double omega = -TDatabase::ParamDB->SSMUM_ROT_PER_SECOND*2*Pi;
    double angle = TDatabase::ParamDB->SSMUM_ANGLE;
    omega = 0;
    value = (2*Pi + omega) * sin(2*Pi*Param-angle);
    //value = (2*Pi + omega) * sin(2*Pi*Param);
    //value = 0;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
    double omega = -TDatabase::ParamDB->SSMUM_ROT_PER_SECOND*2*Pi;
    double angle = TDatabase::ParamDB->SSMUM_ANGLE;

    omega = 0;
    value = -(2*Pi + omega)*cos(2*Pi*Param-angle);
    //value = -(2*Pi + omega)*cos(2*Pi*Param);
    //value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  static double nu = 1/TDatabase::ParamDB->RE_NR;
   int i;
  double *coeff, x, y; 
 
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = nu;
    // rhs 
    coeff[1] = 0;
    coeff[2] = 0;
  }
}


