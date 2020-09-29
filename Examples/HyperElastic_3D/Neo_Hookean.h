//Hyperelastic problem 
// u(x,y,z) = displacement
 

void ExampleFile()
{
  OutPut("Example: Neo_Hookean.h" << endl) ;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int BdComp, double x, double y, double z, BoundCond &cond)
{
//     if(BdComp == 0)
//         cond = NEUMANN;
//     else if(BdComp == 1)
        cond = DIRICHLET;
}

void U1BoundValue(int BdComp,  double x, double y, double z, double &value)
{
  //cout<<"Boundary component 1 "<<BdComp<<" "<<x<<" "<<y<<" "<<z<<endl;
  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value = 100;
            break;
    case 2:  
            value = 0;
            break;
    case 3: value = 0;
            break;
    case 4: value = 0;
            break;
    case 5: value = 0;
            break;
    default: cout << "wrong boundary part number 1" << endl;
    //value = 0.5;
  }
}

void U2BoundValue(int BdComp,  double x, double y, double z, double &value)
{
  switch(BdComp)
  {
   case 0: value = 0;
            break;
    case 1: value = 100;
             break;
    case 2:  
            value = 0;
            break;
    case 3: value = 0;
            break;
    case 4: value = 0;
            break;
    case 5: value = 0;
            break;
    default: cout << "wrong boundary part number 1" << endl;
 }
}

void U3BoundValue(int BdComp,  double x, double y, double z, double &value)
{
  //cout<<"Boundary component 3 "<<BdComp<<" "<<x<<" "<<y<<" "<<z<<endl;
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
  static double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    /*coeff[1] = 0.01; // f1
    coeff[2] = 0.01; // f2
    coeff[3] = 0.01; // f3   */ 
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
    coeff[3] = 0; // f3   
  }
}

