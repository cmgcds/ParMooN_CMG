// Navier-Stokes problem, Driven cavity
//
// u(x,y) = unknown
// p(x,y) = unknown

#include <iostream>

void ExampleFile()
{
  OutPut("Example: flow_Cylinder.h" << endl);
}
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
  // ---------- GMSH  ------------------- //
  // INLET = 0
  // OUTLET = 1
  // walls = 2
  // cylinder = 3

  // -- PARMOON MESH ------------ //
  // inlet = 3 , walls - 0,2 , outlet 1, cylinder - 4
  if (i == 0) //if (i == 1)
    cond = DIRICHLET;
  else if (i == 2 || i == 3)
    cond = DIRICHLET;
  else if (i == 4)
    cond = DIRICHLET;
  else
    cond = NEUMANN;

  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  // INLET = 0 // OUTLET = 1 // walls = 2 // cylinder = 3
  // if (Param > 0 ) cout << "Param Exists" << endl;
  switch (BdComp)
  {
  case 0:
    value = 0; //value = -(Param)*(Param-1)*4 ;
    break;
  case 1:
    value = 0;
    break;
  case 2:
    value = 0;
    break;
  case 3:
    if (abs(Param - 0) < 1e-12 || abs(Param - 1.0) < 1e-12)
      value = 0;
    else
      value = Param * ( 1 - Param  ) * 1 * 4;
    break;

  case 4:
    value = 0;
    break;
  default:
    cout << "wrong boundary part number" << endl;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  static double eps = 1 / TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for (i = 0; i < n_points; i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
  }
}
