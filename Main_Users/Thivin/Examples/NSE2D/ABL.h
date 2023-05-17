// Navier-Stokes problem, Driven cavity
// 
// u(x,y) = unknown
// p(x,y) = unknown

void ExampleFile()
{
  OutPut("Example: flow_Cylinder.h" << endl) ;
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
   // INLET = 0
   // OUTLET = 1
   // walls = 2
   // cylinder = 3   
    if (i==0) //if (i == 1)
        cond = DIRICHLET;
    else if ( i ==1)
        cond = NEUMANN;
    else if ( i == 2)
        cond = NEUMANN;
    else if ( i == 3)
        cond = DIRICHLET;
    else
    {
        cout << "wrong boundary part number" << endl;
        exit(0);
    }
}

void U1BoundValue(int BdComp, double Param, double &value)
{
       // INLET = 0 // OUTLET = 1 // walls = 2 // cylinder = 3
    // if (Param > 0 ) cout << "Param Exists" << endl;
  switch(BdComp)
  {
    case 0: 
            // Convert the param value, which is fom 0 to 1 to actual x value which is from 0 to 500
            value = 0;
            break;
    case 1: value = 0;
            break;
    case 2: value =0;
            break;
    case 3: 
        double actual_y = Param * 500;
        
        double uStar = 0.3;
        double kappa = 0.41;
        double z0 = 4;
        // Compute natural log
        // q: syntax for natural log in C++
        // a: 

        double u = (uStar/kappa) * log( (actual_y + z0)/z0);
        cout <<  " Param : " << Param << " ;; actual_y: " << actual_y << " ;; u : " << u << endl;
        value = u;
        break;  
    
       
    default: cout << "wrong boundary part number" << endl;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
 switch(BdComp)
  {
    case 0: value = 0 ; //value = -(Param)*(Param-1)*4 ;
            break;
    case 1: value = 0;
            break;
    case 2: value =0;
            break;
    case 3: value = 0;
            break;
    case 4: value = 0;
            break;  
    default: cout << "wrong boundary part number" << endl;
  }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
  }
}

