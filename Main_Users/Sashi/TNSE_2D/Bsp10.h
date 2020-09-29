// Stokes problem

// 

// u(x,y) = (t^3 y^2, t^2 x)

// p(x,y) = t*x+y-1/2*t-1/2



void ExampleFile()

{

  OutPut("Example: Bsp10.h" << endl);

}



// ========================================================================

// initial solution

// ========================================================================

void InitialU1(double x, double y, double *values)

{

  double t=TDatabase::TimeDB->CURRENTTIME;



  values[0] = t*t*t*y*y;

}



void InitialU2(double x, double y, double *values)

{

  double t=TDatabase::TimeDB->CURRENTTIME;



  values[0] = t*t*x;

}



void InitialP(double x, double y, double *values)

{

  double t=TDatabase::TimeDB->CURRENTTIME;



  values[0] = t*x+y-0.5*(1+t);

}





// ========================================================================

// exact solution

// ========================================================================

void ExactU1(double x, double y, double *values)

{

  double t=TDatabase::TimeDB->CURRENTTIME;



  values[0] = t*t*t*y*y;

  values[1] = 0;

  values[2] = 2*t*t*t*y;

  values[3] = 2*t*t*t;

}



void ExactU2(double x, double y, double *values)

{

  double t=TDatabase::TimeDB->CURRENTTIME;



  values[0] = t*t*x;

  values[1] = t*t;

  values[2] = 0;

  values[3] = 0;

}



void ExactP(double x, double y, double *values)

{

  double t=TDatabase::TimeDB->CURRENTTIME;



  values[0] = t*x+y-0.5*(1+t);

  values[1] = t;

  values[2] = 1;

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

  double t=TDatabase::TimeDB->CURRENTTIME;



  switch(BdComp)

  {

    case 0:

      value = 0;

    break;

    case 1:

      value = t*t*t*Param*Param;

    break;

    case 2:

      value = t*t*t;

    break;

    case 3:

      value = t*t*t*(1-Param)*(1-Param);

    break;

  }

}



void U2BoundValue(int BdComp, double Param, double &value)

{

  double t = TDatabase::TimeDB->CURRENTTIME;



  switch(BdComp)

  {

    case 0:

      value = t*t*Param;

    break;

    case 1:

      value = t*t;

    break;

    case 2:

      value = t*t*(1-Param);

    break;

    case 3:

      value = 0;

    break;

  }

}



void U1BoundValue_diff(int BdComp, double Param, double &value)

{

  double t=TDatabase::TimeDB->CURRENTTIME;



  switch(BdComp)

  {

    case 0:

      value = 0;

    break;

    case 1:

      value = 3*t*t*Param*Param;

    break;

    case 2:

      value = 3*t*t;

    break;

    case 3:

      value = 3*t*t*(1-Param)*(1-Param);

    break;

  }

}



void U2BoundValue_diff(int BdComp, double Param, double &value)

{

  double t = TDatabase::TimeDB->CURRENTTIME;



  switch(BdComp)

  {

    case 0:

      value = 2*t*Param;

    break;

    case 1:

      value = 2*t;

    break;

    case 2:

      value = 2*t*(1-Param);

    break;

    case 3:

      value = 0;

    break;

  }

}



// ========================================================================

// coefficients for Stokes form: A, B1, B2, f1, f2

// ========================================================================

void LinCoeffs(int n_points, double *X, double *Y,

               double **parameters, double **coeffs)

{

  static double nu = 1/TDatabase::ParamDB->RE_NR;

  double t = TDatabase::TimeDB->CURRENTTIME;

  double t1, t2, t3, t5, t6, t8, t9, t10, t14, t16;

  double t17, t22, t23, t35, t38, t40;

  int i;

  double *coeff, x, y;



  for(i=0;i<n_points;i++)

  {

    coeff = coeffs[i];

    x = X[i];

    y = Y[i];



    coeff[0] = nu;



    coeff[1] = 3*t*t*y*y-nu*2*t*t*t+2*t*t*t*t*t*x*y+t;

    coeff[2] = 2*t*x+t*t*t*t*t*y*y+1;

    coeff[3] = 6*t*y*y - nu*6*t*t + 10*t*t*t*t*x*y + 1;

    coeff[4] = 2*x + 5*t*t*t*t*y*y;

  }

}





