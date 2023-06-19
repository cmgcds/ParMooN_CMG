// Navier-Stokes problem, Driven Cavity
// 
// u(x,y) = ?
// p(x,y) = ?

void ExampleFile()
{
  OutPut("Example: Driven.h" << endl) ;
}
// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 0.0;
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
  cond = DIRICHLET;

  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;
  double eps = 1e-8;
  
  double y = Param * 0.41 ;
  
  switch(BdComp)
  {
    case 3: case 1:  
            value= pow(0.41,-2) * sin(Pi*t/8) * (6*y*(0.41-y));
            // if(BdComp == 1)
            //   cout <<"BDcomp : " << BdComp <<  " Param Values : " << Param  << " y : " << y << " Val " << value << endl;
            break;
    case 0: 
            value=0;
            break;

    case 2: 
            value=0;
            break;
    case 4: 
            value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff, x, y; 
  static double eps=1/TDatabase::ParamDB->RE_NR;
  double t = TDatabase::TimeDB->CURRENTTIME;

  double v = sin(Pi*t/8) ;
  double reNr_new = v/0.001;
  // cout << "Reynolds Number : " << reNr_new << endl;
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = 0.001;

    coeff[1] = 0;  // f1
    coeff[2] = 0;  // f2
  }
}


